import os
import time
from pathlib import Path
import xarray as xr
import numpy as np
import pandas as pd
from numba import jit, prange
import multiprocessing
from functools import partial
import warnings

warnings.filterwarnings('ignore')

# --- Physics Engine (Compiled with Numba) ---

@jit(nopython=True)
def es_calc_full(temp):
    return 6.112 * np.exp(17.67 * (temp - 273.15) / (temp - 273.15 + 243.5))

@jit(nopython=True)
def theta_e_calc(temp, q, p_hpa):
    r = q / (1.0 - q)
    ev_hPa = max(p_hpa * r / (0.622 + r), 1e-8)
    TL = (2840. / ((3.5 * np.log(temp)) - (np.log(ev_hPa)) - 4.805)) + 55.
    chi_e = 0.2854 * (1. - (0.28 * r))
    return temp * ((1000.0 / p_hpa)**chi_e) * np.exp(((3.376 / TL) - 0.00254) * r * 1000.0 * (1. + (0.81 * r)))

@jit(nopython=True)
def trapz_1d(var, p):
    if len(var) < 2:
        return np.nan
    dp_total = p[-1] - p[0]
    if abs(dp_total) < 1e-6:
        return np.nan
    res = 0.0
    for z in range(1, len(var)):
        res += 0.5 * (var[z - 1] + var[z]) * (p[z] - p[z - 1])
    return res / dp_total

@jit(nopython=True, parallel=True)
def calculate_grid_metrics(T_3d, q_3d, p_level, sp_2d, T2m_2d, q2m_2d):
    nz, ny, nx = T_3d.shape
    out_bl = np.full((ny, nx), np.nan, dtype=np.float32)
    out_lt = np.full((ny, nx), np.nan, dtype=np.float32)
    out_sat_lt = np.full((ny, nx), np.nan, dtype=np.float32)

    for j in prange(ny):
        for i in range(nx):
            sfc_p = sp_2d[j, i]
            pbl_p = sfc_p - 100.0
            if pbl_p < 500.0:
                continue

            p_1d = np.zeros(nz + 1)
            t_1d = np.zeros(nz + 1)
            q_1d = np.zeros(nz + 1)

            p_1d[0] = sfc_p
            t_1d[0] = T2m_2d[j, i]
            q_1d[0] = q2m_2d[j, i]
            for k in range(nz):
                p_1d[k + 1] = p_level[k]
                t_1d[k + 1] = T_3d[k, j, i]
                q_1d[k + 1] = q_3d[k, j, i]

            th_e = np.zeros(nz + 1)
            th_e_s = np.zeros(nz + 1)
            for k in range(nz + 1):
                th_e[k] = theta_e_calc(t_1d[k], q_1d[k], p_1d[k])
                es = es_calc_full(t_1d[k])
                qs = (0.622 * es) / (p_1d[k] + ((0.622 - 1.0) * es))
                th_e_s[k] = theta_e_calc(t_1d[k], qs, p_1d[k])

            m_bl = (p_1d <= sfc_p + 0.1) & (p_1d >= pbl_p - 0.1)
            m_lt = (p_1d <= pbl_p + 0.1) & (p_1d >= 500.0 - 0.1)
            out_bl[j, i] = trapz_1d(th_e[m_bl], p_1d[m_bl])
            out_lt[j, i] = trapz_1d(th_e[m_lt], p_1d[m_lt])
            out_sat_lt[j, i] = trapz_1d(th_e_s[m_lt], p_1d[m_lt])

    return out_bl, out_lt, out_sat_lt

# --- CFTime-safe timestamp formatting ---

def time_to_stamp(t):
    """
    Convert a time-like object to YYYYMMDDHHMM string.
    Works for cftime objects and falls back to pandas for numpy datetime64 / strings.
    """
    try:
        return f"{t.year:04d}{t.month:02d}{t.day:02d}{t.hour:02d}{t.minute:02d}"
    except Exception:
        return pd.to_datetime(t).strftime("%Y%m%d%H%M")

# --- WORKER: Named by Timestamp ---

def process_worker(task, ds_ta, ds_hus, temp_dir):
    """Processes a single timestep and saves it to a timestamped file."""
    t_idx, timestamp_str = task
    try:
        temp_file = temp_dir / f"{timestamp_str}.nc"
        if temp_file.exists():
            return temp_file

        # NOTE: Keeping your original logic; only cftime fix requested.
        with ds_ta.isel(time=t_idx) as dta, \
             ds_hus.isel(time=t_idx) as dhus:

            p_lev = dta.plev.values.copy()
            if dta.plev.units != "hPa":
                p_lev /= 100.0

            # Use surface logic
            T = dta.values
            q = dhus.values
            sp = np.full(T.shape[1:], 1000.0)

            bl, lt, sat_lt = calculate_grid_metrics(T, q, p_lev, sp, T[0, :, :], q[0, :, :])

            ds = xr.Dataset(
                data_vars={
                    'thetae_bl': (['lat', 'lon'], bl.astype(np.float32)),
                    'thetae_lt': (['lat', 'lon'], lt.astype(np.float32)),
                    'thetae_lt_sat': (['lat', 'lon'], sat_lt.astype(np.float32)),
                },
                coords={'time': dta.time.values, 'lat': dta.lat.values, 'lon': dta.lon.values}
            )
            ds.to_netcdf(temp_file)
            return temp_file

    except Exception as e:
        print(f"Error at index {t_idx} ({timestamp_str}): {e}")
        return None

# --- Main Logic ---

def process_thetae_layers(ds_var3d, out_dir, ntime=None, num_workers=12):

    start_time = time.time()
    out_dir = Path(out_dir)
    temp_dir = out_dir / "temp_steps"
    temp_dir.mkdir(parents=True, exist_ok=True)

    # Pre-extract timestamps for worker naming
    times = ds_var3d.time.values
    if ntime is not None:
        times = times[:ntime]

    # Create list of (index, timestamp_string) -- CFTime safe
    tasks = [(i, time_to_stamp(t)) for i, t in enumerate(times)]

    print(f"Starting parallel processing for {len(tasks)} steps...")

    # Step 1: Compute and Save Individual files
    ctx = multiprocessing.get_context('spawn')
    with ctx.Pool(processes=num_workers) as pool:
        worker_func = partial(
            process_worker,
            ds_ta=ds_var3d['ta'],
            ds_hus=ds_var3d['hus'],
            temp_dir=temp_dir
        )
        pool.map(worker_func, tasks)

    # Step 2: Monthly Merge using glob pattern
    # Find all unique YYYYMM from the temp files
    all_temp_files = list(temp_dir.glob("*.nc"))
    unique_months = sorted(list(set([f.name[:6] for f in all_temp_files])))

    print(f"Merging into {len(unique_months)} monthly files...")

    for ym in unique_months:
        temp_dir.mkdir(parents=True, exist_ok=True)
        year, month = ym[:4], ym[4:]
        year_dir = out_dir / f'{year}'
        year_dir.mkdir(parents=True, exist_ok=True)
        monthly_output = year_dir / f"layer_averaged_thetae_buoyancy.{year}.{month}.nc"

        # This is the speed trick: open_mfdataset only looks for files of ONE month
        ds_month = xr.open_mfdataset(
            str(temp_dir / f"{ym}*.nc"),
            combine='nested',
            concat_dim='time'
        )

        # Add metadata/depths
        ds_month['depth_pb'] = 100
        ds_month['depth_lt'] = 400

        ds_month.to_netcdf(monthly_output)
        ds_month.close()
        print(f"Completed: {monthly_output.name}")

    print(f"Total time elapsed: {time.time() - start_time:.2f}s")
    # remove temporary files
    os.system(f'rm -r {str(temp_dir)}')

if __name__ == "__main__":

    # Test paths - Replace with your actual files
    data_dir = Path('/pscratch/sd/w/wmtsai/MDTF_POD_diags/inputdata/model/MPI-ESM1-2-HR_historical_r1i1p1f1_gn_200501010130-200912312230/6hr')
    ta_f = data_dir / 'MPI-ESM1-2-HR_historical_r1i1p1f1_gn_200501010600-201001010000.ta.6hr.nc'
    hus_f = data_dir / 'MPI-ESM1-2-HR_historical_r1i1p1f1_gn_200501010600-201001010000.hus.6hr.nc'

    out_path = "/pscratch/sd/w/wmtsai/mdtf/wkdir/MDTF_output/MCS_precip_buoy_stats/model/netCDF/layer_averaged_thetae"

    # Open dataset to pass into the function
    ds_test = xr.open_dataset(ta_f, chunks={'time': 1})
    ds_hus = xr.open_dataset(hus_f, chunks={'time': 1})
    ds_input = xr.Dataset({'ta': ds_test.ta, 'hus': ds_hus.hus})

    # Run test
    process_thetae_layers(ds_input, out_path, ntime=100, num_workers=12)