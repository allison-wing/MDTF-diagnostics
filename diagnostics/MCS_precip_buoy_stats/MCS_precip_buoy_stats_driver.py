import os
import sys
import subprocess
import time
import warnings
from pathlib import Path
import json5
import xarray as xr

# Suppress unnecessary warnings
warnings.filterwarnings('ignore')

def print_banner(message):
    print('\n' + '-'*50)
    print(f"Step: {message}")
    print('-'*50)

def read_settings(opts):
    def to_bool(val):
        if isinstance(val, bool):
            return val
        if val is None:
            return False
        return str(val).strip().lower() == "true"

    settings = {
        "run_mcs_id":     to_bool(opts.get("run_mcs_identification", False)),
        "run_precip_st":  to_bool(opts.get("run_precip_statistics", False)),
        "plot_precip_st": to_bool(opts.get("plot_precip_statistics", False)),
        "run_buoy_calc":  to_bool(opts.get("run_buoyancy_calculation", False)),
        "run_buoy_st":    to_bool(opts.get("run_buoyancy_statistics", False)),
        "plot_buoy_st":   to_bool(opts.get("plot_buoyancy_statistics", False)),
        "num_workers": (opts.get("num_workers", 8)) # default uses 8 workers to accelerate the POD 
    }
    # Derived bounds (expects keys to exist; will raise KeyError if missing)
    settings["lat_bounds"] = slice(opts["latitude_min"], opts["latitude_max"])
    return settings

if __name__ == "__main__":

    pod_name = "MCS_precip_buoy_stats"

    start_time_execute = time.time()
    #========================================================================
    #               STEP: Environment and Path Configuration 
    #=========================================================================
    pod_dir = Path(os.getenv('CODE_ROOT')) / f'diagnostics/{pod_name}'
    utils_dir = pod_dir / 'mcs_utils'
    work_dir = Path(os.environ["WORK_DIR"]) # i.e., /wkdir/MDTF_output/MCS_precip_buoy_stats
    #work_dir = Path('scratch/wmtsai/mdtf_miniforge/wkdir/MDTF_output/MCS_precip_buoy_stats')
    # Add utils to the system path to import modules saved under "/mcs_utils" 
    sys.path.append(str(utils_dir))
    from process_PyFLEXTRKR_MCSmask_writeout_parallel import dask_write_cloudid_PyFLEXTRKR
    from process_layer_thetae_writeout import process_thetae_layers

    #========================================================================
    #               STEP: Load Settings and Preprocessed Variables
    #========================================================================
    # Read the optional settings for step controls described in "settings.jsonc"
    # The default for each step is "True" to process every step below. Any change
    # is not required unless you want to test something step by step.
    settings_path = pod_dir / 'settings.jsonc'
    with open(settings_path, 'r') as f:
        runjob_settings = json5.load(f)    
    opts = runjob_settings['pod_options']
    settings = read_settings(opts) 

    # Read the MDTF-preprocessed variables under "wkdir/MDTF_output/casename_dir"
    # The files are saved as environment variables and read by os.environ['VAR_NAME']
    # All varialbes are then read using xarray and merged into 2-D and 3-D datasets
    ds_var2d = xr.merge([
            xr.open_dataset(os.environ['PR_FILE']), 
            xr.open_dataset(os.environ['RLUT_FILE'])
        ]).sel(lat=settings['lat_bounds'])
    
    ds_var3d = xr.merge([
            xr.open_dataset(os.environ["TA_FILE"]), 
            xr.open_dataset(os.environ["HUS_FILE"])
        ]).sel(lat=settings['lat_bounds'])

    #========================================================================   
    #                        STEP 1: MCS Identification 
    #========================================================================
    # Identify 2-D MCS masks based on a "snapshot-based" approach using thresholds
    # similar to PyFLEXTRKR (see ./doc/MCS_precip_buoy_stats.rst). 
    # 
    # Users can skip this step if they have alternative data (e.g., tracked MCSs from trackers).
    #
    # 1. In this case, set "run_mcs_identification" = False and make sure that you 
    #    create the corresponding directory "MDTF_output/model/netCDF/MCS_identifiers/" ahead
    #    and save the data following the directory structure described in "MCS_precip_buoy_stats.rst"
    #    before running "./mdtf -f /tempelate/runtime_config.yml"
    #
    # 2. Three variables are required in a single MCS data file: 
    #    - "mcs_flag": 2-D MCS binary mask
    #    - "precipitation": precipitation rate (mm/hr)
    #    - "tb": brightness temperature (K) 

    if settings['run_mcs_id']:
        print_banner("Running snapshot-based MCS identification")
        start_year = int(str(opts['start_time'])[:4])
        end_year = int(str(opts['end_time'])[:4])
        # Create "MCS_identifiers" directory and start identifying MCSs
        # Subdaily outputs will be saved under "MCS_identifiers/$year"
        for year in range(start_year, end_year + 1):
            dask_write_cloudid_PyFLEXTRKR(ds_var2d, year)
    else:
        # Check for existing data when skipping the identification step
        mcsmask_dir = work_dir / 'model/netCDF/MCS_identifiers/'
        if not (mcsmask_dir.exists() and any(mcsmask_dir.iterdir())):
            raise FileNotFoundError(
                f"No MCS mask files found in {mcsmask_dir}. Check settings.jsonc or see if alternative Masks are provided."
            )
        print("Existing MCS masks found. Skipping identification.")

    #========================================================================   
    #      STEP 2: Calculate MCS-associated Precipitation Statistics
    #========================================================================
    # This step writes out diagnostics outputs under "MDTF_output/MCS_precip_buoy_stats/model/netCDF/stats"
    # "precip_mcsstats_regridded_monthly.$year.$regrid.nc": MCS frequency, mean precip, accumulative precip,
    # accumulative MCS precip, ... etc. 
    # - Requried steps being run: STEP 1

    if settings['run_precip_st']:
        print_banner("Start MCS-precipitation diagnostics")
        try:
            subprocess.run(['python', str(utils_dir / 'process_prmcs_stats_writeout.py')], check=True)
        except subprocess.CalledProcessError:
            print("Error: process_prmcs_stats_writeout.py failed. If 'run_mcs_identification == False', check the required \
                   variables, 'mcs_flag' and 'precipitation', are provided if alternative data is used." )
            sys.exit(1)

    #========================================================================   
    #      STEP 3: Plot precipitaiton and MCS-associated maps 
    #========================================================================
    # This step generates "MDTF_output/MCS_precip_buoy_stats/fig" and saves diagnostics figures
    # - Requried steps being run: STEP 2 

    if settings['plot_precip_st']:
        print_banner("Plot precipitation and MCS-associated maps")
        try:
            subprocess.run(['python', str(utils_dir / 'plot_precip_stats.py')], check=True)
        except subprocess.CalledProcessError:
            print("Error: plot_precip_stats.py failed. Check 'MDTF_output/MCS_precip_buoy_stats/model/netCDF/stats'")
            sys.exit(1)
    #========================================================================   
    #      STEP 4: Calculate the Lower-Tropospheric Buoyancy 
    #========================================================================
    # This step calculate the two-layer buoyancy components (see ./doc/MCS_precip_buoy_stats.rst)
    # and save the intermediate data under "MDTF_output/MCS_precip_buoy_stats/model/netCDF/layer_averaged_thetae"
    # ***Note: Generated files ~ 500 M per year given 1-deg. and 6-hourly resolutions.
    # Users should be aware of available storage before running the POD. 

    if settings['run_buoy_calc']:
        print_banner("Calculate low-tropospheric buoyancy")
        thetae_dir = work_dir / 'model/netCDF/layer_averaged_thetae/'
        thetae_dir.mkdir(parents=True, exist_ok=True)
        try:
            process_thetae_layers(ds_var3d, thetae_dir, ntime=None, num_workers=settings['num_workers']
                                  )
        except Exception as e:
            print(f"Error during buoyancy calculation: {e}")
            sys.exit(1)

    #========================================================================   
    #      STEP 5: Calculate the Joint Histograms of Buoyancy Components  
    #========================================================================
    # This step generates the diagnostics outputs under "MDTF_output/MCS_precip_buoy_stats/model/netCDF/stats"
    # - Joint histograms of BL,CAPE and BL,SUBSAT and conditional precipitation associated with different 
    #   conditions: MCS, Non-MCS Deep, and others (see ./doc/MCS_precip_buoy_stats.rst)
    # - Requried steps being run: STEP 4 and STEP 1

    if settings['run_buoy_st']:
        print_banner("Calculate histogram statistics")
        try:
            subprocess.run(['python', str(utils_dir / 'process_BLcapesubsat_regions_unified.multiprocess.py'), 'MDTF'], check=True)
        except subprocess.CalledProcessError:
            print("Error: process_BLcapesubsat_regions_unified.multiprocess.py failed.")
            sys.exit(1)

    #========================================================================   
    #      STEP 6: Plot Joint Histograms and Conditional Precipitaiton 
    #========================================================================
    # This step generates "MDTF_output/MCS_precip_buoy_stats/fig" and saves diagnostics figures
    # - Requried steps being run: STEP 5

    if settings['plot_buoy_st']:
        print_banner("Plot 2-D buoyancy-precip statistics")
        try:
            subprocess.run(['python', str(utils_dir / 'plot_BLprecip_relations_2D.py')], check=True)
        except subprocess.CalledProcessError:
            print("Error: plot_BLprecip_relations_2D.py failed. Check 'MDTF_output/MCS_precip_buoy_stats/model/netCDF/stats'")
            sys.exit(1)
            
    # --- Finalize ---
    execution_time = time.time() - start_time_execute
    print(f"\n--- Total Execution Time: {execution_time:.2f} seconds ---")
    print("POD execution finished successfully!")
    sys.exit(0)