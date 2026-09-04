import os
import json5
from pathlib import Path
import xarray as xr
import numpy as np
import pandas as pd
from datetime import datetime
import multiprocessing
from functools import partial
import xesmf as xe  # Required for regridding
import warnings

warnings.filterwarnings('ignore')

def build_target_grid(lat_min, lat_max, res=0.5):
    """Creates a target grid for xESMF with 0.5-degree resolution."""
    # Centers of the bins
    new_lats = np.arange(lat_min + res/2, lat_max, res)
    new_lons = np.arange(0 + res/2, 360, res)
    
    ds_target = xr.Dataset(
        coords={
            'lat': (['lat'], new_lats, {'units': 'degrees_north'}),
            'lon': (['lon'], new_lons, {'units': 'degrees_east'}),
        }
    )
    return ds_target

def coordinates_processors(data):
    if 'longitude' in data.coords:
        data = data.rename({'longitude': 'lon', 'latitude': 'lat'})
    if data.lat[1] - data.lat[0] < 0:
        data = data.reindex(lat=data.lat[::-1])
    if data.lon.min() < 0:
        data.coords['lon'] = (data.coords['lon'] + 360) % 360
        data = data.sortby('lon')
    return data

def process_single_year(year, stats_dir, mcs_id_dir, lat_bounds, p_threshold, regrid_res):
    lat_min, lat_max = lat_bounds
    month_list = np.arange(1, 13)
    
    mcsfiles = sorted(list((mcs_id_dir / f'{year}').glob(f'cloudid_PyFLEXTRKR_mcs_{year}*.nc')))
    if not mcsfiles:
        return None
    
    # 1. Setup Target Grid and Regridder
    ds_target = build_target_grid(lat_min, lat_max, res=regrid_res)
    
    # Open first file to initialize regridder (weights)
    with xr.open_dataset(mcsfiles[0]) as ds_sample:
        ds_sample = coordinates_processors(ds_sample).sel(lat=slice(lat_min, lat_max))
        # xESMF conservative regridding requires bounds. 
        # cf_xarray or xESMF's utility can add them if missing.
        ds_sample = ds_sample.cf.add_bounds(['lat', 'lon'])
        ds_target = ds_target.cf.add_bounds(['lat', 'lon'])
        
        regridder = xe.Regridder(ds_sample, ds_target, method='conservative', periodic=True)

    # Initialize storage for regridded data (now on target grid dims)
    n_lat, n_lon = len(ds_target.lat), len(ds_target.lon)
    pr_freq = np.zeros((12, n_lat, n_lon), dtype=np.float32)
    pr_mcssum = np.zeros_like(pr_freq)
    pr_totsum = np.zeros_like(pr_freq)
    mcs_freq = np.zeros_like(pr_freq)
    tpr = np.zeros_like(pr_freq)

    ds_year = xr.open_mfdataset(mcsfiles, chunks={'time': 100}, compat='override', coords='minimal')
    ds_year = coordinates_processors(ds_year)
    ds_year = ds_year[['precipitation', 'tb', 'mcs_flag']].sel(lat=slice(lat_min, lat_max))

    for i, month in enumerate(month_list):
        strt = datetime(year, month, 1)
        end = datetime(year, month + 1, 1) if month < 12 else datetime(year + 1, 1, 1)
        
        ds_sub = ds_year.sel(time=slice(strt, end)).compute()
        if len(ds_sub.time) == 0: continue

        valid_mask = ds_sub.tb.notnull()
        pr = ds_sub.precipitation.where(valid_mask).astype(np.float32)
        mcs_flag = ds_sub.mcs_flag.where(valid_mask)
        mcsmask = xr.where(mcs_flag > 0, 1.0, 0.0).where(valid_mask).astype(np.float32)
        effective_times = valid_mask.sum('time').astype(np.float32)
        
        # Calculate statistics on high-res first
        s_pr_freq = (pr > p_threshold).sum('time') / effective_times
        s_pr_mcssum = pr.where(mcsmask == 1).sum('time')
        s_pr_totsum = pr.sum('time')
        s_mcs_freq = mcsmask.sum('time') / effective_times
        s_tpr = pr.mean('time')

        # 2. Apply Regridder to the monthly results
        pr_freq[i] = regridder(s_pr_freq).values
        pr_mcssum[i] = regridder(s_pr_mcssum).values
        pr_totsum[i] = regridder(s_pr_totsum).values
        mcs_freq[i] = regridder(s_mcs_freq).values
        tpr[i] = regridder(s_tpr).values

    # 3. Construct Output Dataset
    ds_out = xr.Dataset(
        data_vars=dict(
            pr_freq=(['year', 'month', 'lat', 'lon'], np.expand_dims(pr_freq, 0)),
            pr_mcssum=(['year', 'month', 'lat', 'lon'], np.expand_dims(pr_mcssum, 0)),
            pr_totsum=(['year', 'month', 'lat', 'lon'], np.expand_dims(pr_totsum, 0)),
            mcs_freq=(['year', 'month', 'lat', 'lon'], np.expand_dims(mcs_freq, 0)),
            tpr=(['year', 'month', 'lat', 'lon'], np.expand_dims(tpr, 0))
        ),
        coords=dict(
            year=[int(year)], 
            month=month_list.astype(np.int32), 
            lat=ds_target.lat.values, 
            lon=ds_target.lon.values
        ),
        attrs=dict(description=f'Regridded {regrid_res}deg PyFLEXTRKR MCS stats')
    ).astype(np.float32)
    
    ds_out.to_netcdf(stats_dir / f'precip_mcsstats_regridded_monthly.{year}.{regrid_res}deg.nc')
    print(f"output saved as : {stats_dir / f'precip_mcsstats_regridded_monthly.{year}.{regrid_res}deg.nc'}")

    return year

def main():

    main_dir = Path(os.getenv('CODE_ROOT', '.'))
    pod_dir = main_dir / 'diagnostics/MCS_precip_buoy_stats'
    settings_path = pod_dir / 'settings.jsonc'
    with open(settings_path, 'r') as f:
        runjob_settings = json5.load(f)
    opts = runjob_settings['pod_options']
    # read regridding resolution from settings.jsonc, default is using 1-degree 
    regrid_res = opts.get('regrid_res', 1) 

    # Path setup stays the same
    #work_dir = Path(os.environ["WORK_DIR"]) # Note!!!! This becomes /wkdir/MDTF_output.v{}/MCS_precip_buoy_stats ! not /wkdir
    work_dir = Path('/scratch/wmtsai/mdtf_miniforge/wkdir/MDTF_output/MCS_precip_buoy_stats')
    mcsmask_dir = work_dir / f'model/netCDF/MCS_identifiers/'    
    stats_dir = work_dir / f'model/netCDF/stats'
    os.makedirs(str(stats_dir), exist_ok=True)
  
   # get a list of year with processed MCS identification
    year_list = sorted([
        int(d.name) for d in mcsmask_dir.iterdir() 
        if d.is_dir() and d.name.isdigit()
    ])
    num_workers = min(len(year_list), 4) # Conservative regridding is memory intensive, lower worker count slightly
    worker_func = partial(process_single_year, stats_dir=stats_dir, mcs_id_dir=mcsmask_dir,
                          lat_bounds=(-30, 30), p_threshold=0.01, regrid_res=regrid_res)

    with multiprocessing.Pool(num_workers) as pool:
        pool.map(worker_func, year_list)

    pool.close()
    pool.join()

if __name__ == "__main__":

   main()