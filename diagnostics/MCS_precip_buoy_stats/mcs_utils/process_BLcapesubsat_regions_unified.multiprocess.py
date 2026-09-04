import os
import sys
import time
import json5
import xarray as xr
import numpy as np
from pathlib import Path
import multiprocessing
from functools import partial
import xesmf as xe
import warnings

warnings.filterwarnings('ignore')

def get_histogram_2d(c_vals, s_vals, p_vals, bins_cape, bins_subsat):
    valid = np.isfinite(c_vals) & np.isfinite(s_vals) & np.isfinite(p_vals)
    empty = np.zeros((len(bins_cape)-1, len(bins_subsat)-1), dtype=float)

    if not np.any(valid):
        return empty, empty

    h, _, _  = np.histogram2d(c_vals[valid], s_vals[valid],
                             bins=[bins_cape, bins_subsat])
    hp, _, _ = np.histogram2d(c_vals[valid], s_vals[valid],
                             bins=[bins_cape, bins_subsat], weights=p_vals[valid])
    return h, hp

def process_single_month(t_file, mcs_dir, landseamask, bins_cape, bins_subsat, ds_out_grid):
    """
    Worker function: Calculates physics at native resolution before regridding to 
    preserve non-linear convective signals.
    """
    # 1. Load Data at Native Resolution
    ds_t = xr.open_dataset(t_file).compute()
    time_str = ds_t.time.dt.strftime('%Y%m').values[0]
    year_str, month_str = time_str[:4], time_str[4:]
   
    mcs_month_files = sorted(list(mcs_dir.glob(f'cloudid_PyFLEXTRKR_mcs_{year_str}{month_str}*.nc')))
    if not mcs_month_files: 
            raise FileNotFoundError(f"No MCS files found for {year_str}-{month_str} in {mcs_dir}")
    
    ds_m = xr.open_mfdataset(mcs_month_files, combine='nested', concat_dim='time', 
                                coords='minimal', compat='override')
    ds_m = ds_m.compute()
    
    # Align high-res datasets
    ds_t, ds_m = xr.align(ds_t, ds_m, join='inner')
    tb_threshold = float(ds_m.attrs.get('tb_threshold_coldanvil', 241))

    # --- 2. Physics Calculation at NATIVE Resolution ---
    # This prevents the 'CAPE dilution' shift seen when averaging T/q before calculation
    delta_pl = 1000 - 100 - 500
    wb = (100 / delta_pl) * np.log((delta_pl + 100) / 100)
    wl = 1 - wb
    
    cape_native = (wb * ((ds_t.thetae_bl - ds_t.thetae_lt_sat) / ds_t.thetae_lt_sat) * 340)
    subsat_native = (wl * ((ds_t.thetae_lt_sat - ds_t.thetae_lt) / ds_t.thetae_lt_sat) * 340)

    # Create a container for native variables to regrid all at once
    ds_phys = xr.Dataset({
        'CAPE': cape_native,
        'SUBSAT': subsat_native,
        'tb': ds_m.tb,
        'mcs_flag': ds_m.mcs_flag,
        'precip': ds_m.precipitation
    })

    # --- 3. Regridding to Analysis Resolution ---
    regridder = xe.Regridder(ds_phys, ds_out_grid, method='conservative', periodic=True)
    ds_reg = regridder(ds_phys)
    
    regridder_landsea = xe.Regridder(landseamask, ds_out_grid, method='conservative', periodic=True)
    mask_reg = regridder_landsea(landseamask)
    ocean_bool = (mask_reg > 0.5).values 
    land_bool = ~ocean_bool

    # Binary mask based on area-averaged flag
    m_flag_binary = xr.where(ds_reg.mcs_flag >= 0.5, 1, 0)

    # 4. Initialize Results
    cats, surfaces = ['mcs', 'deep', 'other'], ['land', 'ocean']
    local_res = {cat: {srf: {'samples': 0, 'prec_sum': 0} 
                 for srf in surfaces} for cat in cats}

    # 5. Sampling Loop
    for t in range(len(ds_reg.time)):
        m_flag = m_flag_binary.isel(time=t).values
        tb = ds_reg.tb.isel(time=t).values
        pr = ds_reg.precip.isel(time=t).values
        CAPE_vals = ds_reg.CAPE.isel(time=t).values
        SUBSAT_vals = ds_reg.SUBSAT.isel(time=t).values
        
        valid_mask = ~np.isnan(m_flag) & ~np.isnan(tb) & ~np.isnan(pr)

        cat_logic = {
            'mcs': (m_flag == 1) & valid_mask,
            'deep': (m_flag == 0) & (tb < tb_threshold) & valid_mask,
            'other': (m_flag == 0) & (tb >= tb_threshold) & valid_mask
        }

        for cat, c_mask in cat_logic.items():
            srf_logic = {'land': c_mask & land_bool, 'ocean': c_mask & ocean_bool}
            for srf, final_mask in srf_logic.items():
                if np.any(final_mask):
                    h, hp = get_histogram_2d(CAPE_vals[final_mask], SUBSAT_vals[final_mask], 
                                                 pr[final_mask],
                                                 bins_cape, bins_subsat)
                    local_res[cat][srf]['samples'] += h
                    local_res[cat][srf]['prec_sum'] += hp
                
    return local_res

def main():

    # command-line MDTF (./mdtf -f ./template/runtime_config.jsonc) or ESNB Jupyter notebook
    running_mode = sys.argv[1] # "MDTF" or "ESNB"

    start_time = time.time()
 
    if running_mode == "MDTF":

        print('Reading settings.jsonc for MDTF')
        work_dir = Path(os.environ['WORK_DIR'])
        main_dir = Path(os.getenv('CODE_ROOT', '.'))
        pod_dir = main_dir / 'diagnostics/MCS_precip_buoy_stats'
        settings_path = pod_dir / 'settings.jsonc'
        with open(settings_path, 'r') as f:
            runjob_settings = json5.load(f)
        opts = runjob_settings['pod_options']
        # read regridding resolution from settings.jsonc, default is using 1-degree 
        regrid_res = opts.get('regrid_res', 1) # default == 1 degree
        dx_str = f"{regrid_res:.1f}"

    elif running_mode == "ESNB":

        print('Reading settings.jsonc for ESNB')
        work_dir = Path(os.environ['WORK_DIR'])
        settings_path = work_dir / 'settings.jsonc'
        with open(settings_path, 'r') as f:
            runjob_settings = json5.load(f)
        opts = runjob_settings['settings']['pod_options']
        # read regridding resolution from settings.jsonc, default is using 1-degree 
        regrid_res = opts.get('regrid_res', 1) # default == 1 degree
        dx_str = f"{regrid_res:.1f}"       

    else:
        raise ValueError("Unknown running mode argument (either MDTF or ESNB)")    

    obsdata_dir = Path(os.environ['CODE_ROOT']).parent / 'inputdata/obs_data' # get obs_data directory
#    work_dir = Path('/pscratch/sd/w/wmtsai/mdtf/wkdir/MDTF_output/MCS_precip_buoy_stats')
    thetae_dir = work_dir / 'model/netCDF/layer_averaged_thetae'
    mcs_dir = work_dir / 'model/netCDF/MCS_identifiers'
    stats_dir = work_dir / 'model/netCDF/stats'
    stats_dir.mkdir(parents=True, exist_ok=True)
    
    year_list = sorted([int(d.name) for d in mcs_dir.iterdir() if d.is_dir() and d.name.isdigit()])
    bins_cape = np.arange(-15, 10.5, 0.5)
    bins_subsat = np.arange(-5, 25.5, 0.5)
    
    ds_out_grid = xr.Dataset({"lat": (["lat"], np.arange(-30, 30+regrid_res, regrid_res)),
                              "lon": (["lon"], np.arange(0, 360, regrid_res))})

    for year in year_list:
        print('proceesing year: ', year)
        thetae_files = sorted(list((thetae_dir / f'{year}').glob('layer_averaged_thetae_buoyancy*.nc')))
        
        # Load and Fix LandSea Mask coordinates
        ds_landsea = xr.open_dataset(obsdata_dir / 'LandSeaMask.0.25deg.nc').sel(lat=slice(-30,30))
        # Standardize longitude to 0-360 if it is -180/180
        if ds_landsea.lon.min() < 0:
            ds_landsea = ds_landsea.assign_coords(lon=(ds_landsea.lon % 360)).sortby('lon')
        landseamask = ds_landsea.landseamask / 100.0
        
        with multiprocessing.Pool(processes=min(len(thetae_files), 6)) as pool:
            worker_func = partial(process_single_month, mcs_dir=mcs_dir / f'{year}', 
                                landseamask=landseamask, bins_cape=bins_cape, 
                                bins_subsat=bins_subsat, ds_out_grid=ds_out_grid)
            results = pool.map(worker_func, thetae_files)

        # --- Aggregate and Save ---
        cats, surfaces = ['mcs', 'deep', 'other'], ['land', 'ocean']
        
        for cat in cats:
            agg_samples = {srf: 0 for srf in surfaces}
            agg_prec    = {srf: 0 for srf in surfaces}

            for res in results:
                if res is None: continue
                for srf in surfaces:
                    agg_samples[srf] += res[cat][srf]['samples']
                    agg_prec[srf]    += res[cat][srf]['prec_sum']

            ds_out = xr.Dataset(
                data_vars=dict(
                    samples=(['surface_type', 'bins_cape', 'bins_subsat'], np.stack([agg_samples['land'], agg_samples['ocean']])),
                    prec_sum=(['surface_type', 'bins_cape', 'bins_subsat'], np.stack([agg_prec['land'], agg_prec['ocean']]))
                    
                ),
                coords=dict(surface_type=surfaces, bins_cape=bins_cape[:-1], bins_subsat=bins_subsat[:-1]),
                attrs=dict(description=f"2D Histograms for {cat}", year=year, resolution=f"{dx_str}deg")
            )
            
            stats_dir.mkdir(parents=True, exist_ok=True)
            outfile = stats_dir / f'BLprec_2Dhist_{cat}_combined.{year}.{dx_str}deg.nc'
            ds_out.to_netcdf(outfile)
            print(f'output saved as : {outfile}')

        pool.close()
        pool.join()

    print(f"Finished. Total Time: {time.time() - start_time:.2f}s")

if __name__ == "__main__":

    main()
