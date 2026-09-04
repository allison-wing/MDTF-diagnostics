import os
import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path
from datetime import datetime
from scipy.ndimage import label, sum as ndi_sum
from dask import delayed, compute

# ---------------------------------------------------------
# 1. Fast Vectorized Worker Function (Pure Numpy)
# ---------------------------------------------------------
def process_timestep_vectorized(tb_2d, pr_2d, lat, lon, time_val, out_dir):
    """
    Processing logic using Scipy to avoid loops.
    Expects pure numpy arrays (not xarray objects) for speed.
    """
    # --- Constants ---
    dlon = (lon[1] - lon[0])
    dlat = (lat[1] - lat[0])
    dx = dlon * dlat * 100 * 100  # Area approx
    
    tb_threshold_core = 225.
    tb_threshold_coldanvil = 246.
    pf_rr_thres = 2 / (dlon*dlon/(0.1*0.1)) 
    pf_link_area_threshold = 1000
    mcs_tb_area_thresh = 40000

    # --- Masks ---
    # Prepare binary masks
    mask_cold = (tb_2d <= tb_threshold_coldanvil) # Boolean
    mask_core = (tb_2d <= tb_threshold_core)
    mask_pf   = (pr_2d >= pf_rr_thres)

    # Label features (clouds)
    # feature_map is an integer array where 0=background, 1=cloud1, 2=cloud2...
    feature_map, num_features = label(mask_cold)
    
    # Early exit if no clouds
    if num_features == 0:
        save_result(out_dir, time_val, lat, lon, 
                    np.zeros_like(tb_2d, dtype='int8'), 
                    np.zeros_like(tb_2d, dtype='int8'), 
                    np.zeros_like(tb_2d, dtype='int8'), 
                    tb_2d, pr_2d)
        return

    # --- Vectorized MCS Logic (Replacing the loop) ---
    labels_index = np.arange(1, num_features + 1)

    # 1. Calculate Area for all clouds at once
    # ndi_sum(input, labels, index) sums the 'input' values over the regions defined by 'labels'
    area_counts = ndi_sum(np.ones_like(feature_map), feature_map, index=labels_index)
    cloud_areas = area_counts * dx

    # 2. Check Overlaps (Sum masks inside cloud regions)
    # Core overlap
    has_core = ndi_sum(mask_core, feature_map, index=labels_index) > 0
    # PF overlap count and area
    pf_counts = ndi_sum(mask_pf, feature_map, index=labels_index)
    has_pf = pf_counts > 0
    pf_areas = pf_counts * dx

    # 3. Apply Criteria (Boolean Indexing)
    # Valid clouds must meet ALL criteria
    valid_mask = (
        (cloud_areas >= mcs_tb_area_thresh) & 
        (has_pf) & 
        (pf_areas > pf_link_area_threshold) & 
        (has_core)
    )
    
    # Get the ID numbers of valid clouds
    valid_ids = labels_index[valid_mask]

    # --- Reconstruct Output Maps ---
    if len(valid_ids) == 0:
        # No valid MCS found after filtering
        save_result(out_dir, time_val, lat, lon, 
                    np.zeros_like(tb_2d, dtype='int8'), 
                    np.zeros_like(tb_2d, dtype='int8'), 
                    np.zeros_like(tb_2d, dtype='int8'), 
                    tb_2d, pr_2d)
        return

    # Create Final Masks using np.isin (Fast look up)
    final_mcs_mask = np.isin(feature_map, valid_ids)
    
    # Filter PF and Core to only exist INSIDE the valid MCS shields
    final_pf_mask = mask_pf & final_mcs_mask
    final_core_mask = mask_core & final_mcs_mask

    # Save
    save_result(out_dir, time_val, lat, lon, 
                final_mcs_mask.astype('int8'), 
                final_pf_mask.astype('int8'), 
                final_core_mask.astype('int8'), 
                tb_2d, pr_2d)

# ---------------------------------------------------------
# 2. Helper: Fast Writer
# ---------------------------------------------------------
def save_result(out_dir, time_val, lat, lon, mcs, pf, core, tb, pr):
    # Same as before, just creates the netcdf
    timestamp = str(time_val)
    # Basic parsing
    try:
        # Fast string slicing for standard ISO format
        yr, month, day, hour, minute = timestamp[:4], timestamp[5:7], timestamp[8:10], timestamp[11:13], timestamp[14:16]
    except:
        ts = pd.to_datetime(timestamp)
        yr, month, day, hour, minute = ts.year, ts.month, ts.day, ts.hour, ts.minute

    fname = f'cloudid_PyFLEXTRKR_mcs_{yr}{month}{day}.{hour}{minute}.nc'
    
# Create dataset only at the last moment
    ds = xr.Dataset(
        data_vars={
            'mcs_flag': (('time', 'lat', 'lon'), mcs[None, ...]),
            'pf_flag':  (('time', 'lat', 'lon'), pf[None, ...]),
            'core_flag':(('time', 'lat', 'lon'), core[None, ...]),
            'tb':        (('time', 'lat', 'lon'), tb[None, ...]),
            'precipitation': (('time', 'lat', 'lon'), pr[None, ...]),
        },
        coords={'time': [time_val], 'lat': lat, 'lon': lon}
    )

    # --- Add Units and Metadata ---
    ds.precipitation.attrs['units'] = 'mm/hr'
    ds.precipitation.attrs['long_name'] = 'Precipitation Rate'

    ds.tb.attrs['units'] = 'K'
    ds.tb.attrs['long_name'] = 'Brightness Temperature'

    # Adding metadata for the masks is also good practice
    ds.mcs_flag.attrs['units'] = 'unitless'
    ds.mcs_flag.attrs['long_name'] = 'MCS Binary Mask'

    # Compression helps write speed for masks (lots of zeros)
    encoding = {v: {'zlib': True, 'complevel': 1} for v in ds.data_vars}
    ds.to_netcdf(out_dir / fname, encoding=encoding)

# ---------------------------------------------------------
# 3. Helper: OLR Calculation (Numpy)
# ---------------------------------------------------------
def olr_to_tb_numpy(OLR):
    a = 1.228
    b = -1.106e-3
    sigma = 5.67e-8 
    tf = (OLR/sigma)**0.25
    term = a**2 + 4*b*tf
    term = np.where(term < 0, 0, term) # Safety
    tb = (-a + np.sqrt(term))/(2*b)
    return tb

# ---------------------------------------------------------
# 4. Main Driver (Uses your original structure)
# ---------------------------------------------------------
def dask_write_cloudid_PyFLEXTRKR(input_data, year):
    
    # A. LOAD DATA ONCE (Original Strategy)
    print("Loading data into memory...")
    input_data = input_data.sel(time=slice(datetime(year,1,1,0), datetime(year+1,12,31,23))).compute() 
    
    # Pre-calculate simple arrays to pass to workers
    # (Avoid passing heavy Xarray objects if possible, pass numpy)
    print("Pre-calculating TB and Pr...")
   
    # We convert entire arrays to numpy at once (Very fast vectorized)
    all_rlut = input_data['rlut'].values
    all_pr   = input_data['pr'].values * 3600. # convert to mm/hr
    all_tb   = olr_to_tb_numpy(all_rlut) # Calculate TB here once
    
    lat = input_data.lat.values
    lon = input_data.lon.values
    times = input_data.time.values
  
    # B. SETUP OUTPUT
    work_dir = Path(os.environ["WORK_DIR"])
    mcsmask_dir = work_dir / f'model/netCDF/MCS_identifiers/{year}'
    os.makedirs(mcsmask_dir, exist_ok=True)
    
    print(f"Processing {len(times)} timesteps...")

    # C. PARALLEL EXECUTION
    # Use dask.delayed on the numpy arrays.
    # Since data is in memory, this is purely CPU bound now.
    tasks = []
    
    for i, t in enumerate(times):
        # Slice numpy arrays (Cheap pointers)
        tb_slice = all_tb[i, :, :]
        pr_slice = all_pr[i, :, :]
        
        # Create delayed task
        # We wrap the function with @delayed or call delayed(func)
        task = delayed(process_timestep_vectorized)(
            tb_slice, pr_slice, lat, lon, t, mcsmask_dir
        )
        tasks.append(task)
    
    # D. COMPUTE
    compute(tasks, scheduler='processes', num_workers=8)
    
    print("Done.")
  
if __name__ == "__main__":
  
    year = 2005
    work_dir = Path('/pscratch/sd/w/wmtsai/mdtf/wkdir/MDTF_output/MCS_precip_buoy_stats')
    mcsmask_dir = work_dir / f'model/netCDF/MCS_identifiers/{year}'    
    stats_dir = work_dir / f'model/netCDF/stats'
    os.makedirs(str(stats_dir), exist_ok=True)
 
    data_dir = work_dir.parent / 'MDTF_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_20050101000000_20060101000000/6hr'
    pr_path = data_dir / 'MPI-ESM1-2-HR_historical_r1i1p1f1_gn.pr.6hr.nc'
    rlut_path = data_dir / 'MPI-ESM1-2-HR_historical_r1i1p1f1_gn.rlut.6hr.nc'
    input_data = xr.merge([xr.open_dataset(pr_path), xr.open_dataset(rlut_path)])
    #input_data = input_data.rename({'latitude':'lat', 'longitude':'lon'})
    dask_write_cloudid_PyFLEXTRKR(input_data, year)