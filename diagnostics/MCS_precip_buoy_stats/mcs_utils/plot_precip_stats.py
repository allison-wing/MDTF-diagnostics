import os
import json5
from pathlib import Path
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from PIL import Image

import warnings
warnings.filterwarnings('ignore')

def plot_prec_clim_seasonal_variability(ds_all):
    """
    Input: ds_all (xarray Dataset) with 'month' coordinate as integers 1-12.
    """
    # --- 1. Internal Statistics Calculation ---
    # Annual Mean per year (average across the 12 numerical months)
    yearly_mean_map = 24 * ds_all.tpr.mean('month') 
    map_data = yearly_mean_map.mean('year')
    
    # Annual Zonal Mean & STD
    zonal_by_year = yearly_mean_map.mean('lon')
    zonal_mean_ann = zonal_by_year.mean('year')
    zonal_std_ann = zonal_by_year.std('year')

    # --- Seasonal Zonal Means (Manual Mapping for 1-12) ---
    # Define seasons based on month integers
    seasons = {
        'DJF': [12, 1, 2],
        'MAM': [3, 4, 5],
        'JJA': [6, 7, 8],
        'SON': [9, 10, 11]
    }
    
    # Calculate daily rate for zonal mean calculation
    ds_daily = ds_all.tpr * 24
    zonal_seasonal = {}
    
    for name, m_list in seasons.items():
        # Select months, then average over month, year, and longitude
        zonal_seasonal[name] = ds_daily.sel(month=m_list).mean(dim=['month', 'year', 'lon'])

    # --- 2. Figure Setup ---
    fig = plt.figure(figsize=(12, 4.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[4, 1], wspace=0.08, top=0.82)

    # Left Panel: Annual Mean Map
    ax = fig.add_subplot(gs[0], projection=ccrs.PlateCarree(180))
    ax.coastlines(linewidth=0.5, color='k')
    ax.add_feature(cfeat.LAND, facecolor='none', edgecolor='grey', linewidth=0.3)
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.5, color='grey', alpha=0.4, linestyle=':')
    gl.top_labels = False; gl.right_labels = False
    gl.xlabel_style = {'size': 8}; gl.ylabel_style = {'size': 8}

    cp = ax.contourf(ds_all.lon, ds_all.lat, map_data,
                     cmap='jet', levels=np.arange(0, 18, 2), extend='max',
                     transform=ccrs.PlateCarree())
 
    ax.contour(ds_all.lon, ds_all.lat, map_data,
                     colors=['k'], levels=np.arange(0, 18, 4), linewidths=1,
                     transform=ccrs.PlateCarree())
    ax.set_ylim([-30, 30])
    ax.set_title(rf'Annual Precipitation ({ds_all.year[0].values}-{ds_all.year[-1].values})', 
                 loc='left', fontweight='bold', fontsize=11, pad=12)

    # Right Panel: Zonal Statistics
    ax_zonal = fig.add_subplot(gs[1])
    
    # Force layout draw to align heights with Cartopy
    fig.canvas.draw()
    pos_map = ax.get_position()
    pos_zonal = ax_zonal.get_position()
    ax_zonal.set_position([pos_zonal.x0, pos_map.y0, pos_zonal.width, pos_map.height])
    
    # 1. Plot Shaded Annual Variability (1-sigma envelope)
    ax_zonal.fill_betweenx(zonal_mean_ann.lat, 
                           zonal_mean_ann - zonal_std_ann, 
                           zonal_mean_ann + zonal_std_ann, 
                           color='gray', alpha=0.25, label=r'Ann 1$\sigma$ IAV')
    
    # 2. Plot Annual Mean Line
    ax_zonal.plot(zonal_mean_ann, zonal_mean_ann.lat, color='black', 
                  linewidth=2, label='Annual Mean', zorder=5)

    # 3. Plot Seasonal Mean Lines
    season_colors = {'DJF': 'tab:blue', 'MAM': 'tab:green', 
                     'JJA': 'tab:red', 'SON': 'tab:orange'}
    
    for sea in ['DJF', 'MAM', 'JJA', 'SON']:
        ax_zonal.plot(zonal_seasonal[sea], zonal_seasonal[sea].lat, 
                      color=season_colors[sea], linewidth=1.2, 
                      label=sea, alpha=0.9)
    
    ax_zonal.set_ylim([-30, 30])
    ax_zonal.set_xlim([0, 10]) # Seasonal peaks are often higher than annual mean
    ax_zonal.set_xlabel('mm/day', fontsize=8)
    ax_zonal.set_yticklabels([]) 
    ax_zonal.grid(True, linestyle=':', alpha=0.5)
    
    # Position legend slightly outside to avoid clutter
    ax_zonal.legend(fontsize=7, loc='upper left', frameon=False, bbox_to_anchor=(1.02, 1))
    ax_zonal.set_title('Zonal Mean & Seasons', fontsize=9, fontweight='bold')
    ax_zonal.tick_params(labelsize=8)

    # --- 3. Colorbar ---
    cbar_width, cbar_height = 0.25, 0.025
    cbar_ax = fig.add_axes([pos_map.x1 - cbar_width, pos_map.y1 + 0.04, cbar_width, cbar_height])
    cbar = plt.colorbar(cp, cax=cbar_ax, orientation='horizontal')
    cbar.ax.tick_params(width=1, length=2, labelsize=7)
    cbar_ax.text(-0.05, 0.5, '(mm/day)', transform=cbar_ax.transAxes, 
                 fontsize=8, va='center', ha='right')

    return fig, ax, ax_zonal

def plot_mcsfreq_clim_seasonal(ds_all):
    """
    Input: ds_all (xarray Dataset) with 'month' coordinate as integers 1-12.
    Plots Annual MCS Frequency Map + Zonal Mean (Annual, IAV, and 4 Seasons).
    """
    # --- 1. Internal Statistics Calculation ---
    # Annual Mean per year (converted to %)
    yearly_mean_map = 100 * ds_all.mcs_freq.mean('month') 
    map_data = yearly_mean_map.mean('year')
    
    # Annual Zonal Mean & STD
    zonal_by_year = yearly_mean_map.mean('lon')
    zonal_mean_ann = zonal_by_year.mean('year')
    zonal_std_ann = zonal_by_year.std('year')

    # Seasonal Zonal Means (Manual Mapping for 1-12)
    seasons = {'DJF': [12, 1, 2], 'MAM': [3, 4, 5], 'JJA': [6, 7, 8], 'SON': [9, 10, 11]}
    zonal_seasonal = {}
    for name, m_list in seasons.items():
        zonal_seasonal[name] = 100 * ds_all.mcs_freq.sel(month=m_list).mean(dim=['month', 'year', 'lon'])

    # --- 2. Figure Setup ---
    fig = plt.figure(figsize=(12, 4.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[4, 1], wspace=0.08, top=0.82)

    ax = fig.add_subplot(gs[0], projection=ccrs.PlateCarree(180))
    ax.coastlines(linewidth=0.5); ax.add_feature(cfeat.LAND, facecolor='none', edgecolor='grey', linewidth=0.3)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.4, linestyle=':')
    gl.top_labels = False; gl.right_labels = False
    gl.xlabel_style = {'size': 8}; gl.ylabel_style = {'size': 8}

    cp = ax.contourf(ds_all.lon, ds_all.lat, map_data, cmap='magma_r', 
                     levels=np.arange(0, 30, 2.5), extend='max', transform=ccrs.PlateCarree())
    ax.set_ylim([-30, 30])
    ax.set_title(rf'MCS Frequency ({ds_all.year[0].values}-{ds_all.year[-1].values})'
                                         , loc='left', fontweight='bold', fontsize=11, pad=12)

    # Right Panel: Zonal Statistics
    ax_zonal = fig.add_subplot(gs[1])
    fig.canvas.draw()
    pos_map = ax.get_position(); pos_zonal = ax_zonal.get_position()
    ax_zonal.set_position([pos_zonal.x0, pos_map.y0, pos_zonal.width, pos_map.height])
    
    # Plot Annual IAV Shading and Mean
    ax_zonal.fill_betweenx(zonal_mean_ann.lat, zonal_mean_ann - zonal_std_ann, 
                           zonal_mean_ann + zonal_std_ann, color='gray', alpha=0.25, label=r'Ann 1$\sigma$ IAV')
    ax_zonal.plot(zonal_mean_ann, zonal_mean_ann.lat, color='black', linewidth=2, label='Annual Mean', zorder=5)

    # Plot Seasonal Lines
    colors = {'DJF': 'tab:blue', 'MAM': 'tab:green', 'JJA': 'tab:red', 'SON': 'tab:orange'}
    for sea in ['DJF', 'MAM', 'JJA', 'SON']:
        ax_zonal.plot(zonal_seasonal[sea], zonal_seasonal[sea].lat, color=colors[sea], linewidth=1.2, label=sea)
    
    ax_zonal.set_ylim([-30, 30]); ax_zonal.set_xlim([0, 15]); ax_zonal.set_xlabel('(%)', fontsize=8)
    ax_zonal.set_yticklabels([]); ax_zonal.grid(True, linestyle=':', alpha=0.5)
    ax_zonal.legend(fontsize=7, loc='upper left', frameon=False, bbox_to_anchor=(1.02, 1))
    ax_zonal.set_title('Zonal Mean & Seasons', fontsize=9, fontweight='bold')
    ax_zonal.tick_params(labelsize=8)

    # Colorbar
    cbar_ax = fig.add_axes([pos_map.x1 - 0.25, pos_map.y1 + 0.04, 0.25, 0.025])
    cbar = plt.colorbar(cp, cax=cbar_ax, orientation='horizontal')
    cbar.ax.tick_params(width=1, length=2, labelsize=7)
    cbar_ax.text(-0.05, 0.5, '(%)', transform=cbar_ax.transAxes, fontsize=8, va='center', ha='right')

    return fig, ax, ax_zonal

def plot_mcsprecip_contr_seasonal(ds_all):
    """
    Input: ds_all (xarray Dataset) with 'month' coordinate as integers 1-12.
    Plots MCS Contribution to Total Precip (%) + Seasonal Zonal Mean.
    """
    colors2 = plt.cm.jet(np.linspace(0.3, 0.9, 8))
    colors1 = plt.cm.binary_r(np.linspace(1, 0.8, 2))
    colors = np.vstack((colors1, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    # --- 1. Internal Statistics Calculation ---
    # Annual Contribution per year: (Sum of MCS Precip / Sum of Total Precip)
    # Calculated per year to get IAV
    yearly_mcs = ds_all.pr_mcssum.sum('month')
    yearly_tot = ds_all.pr_totsum.sum('month')
    yearly_contr = 100 * (yearly_mcs / yearly_tot)
    
    map_data = yearly_contr.mean('year')
    zonal_by_year = yearly_contr.mean('lon')
    zonal_mean_ann = zonal_by_year.mean('year')
    zonal_std_ann = zonal_by_year.std('year')

    # Seasonal Zonal Means
    seasons = {'DJF': [12, 1, 2], 'MAM': [3, 4, 5], 'JJA': [6, 7, 8], 'SON': [9, 10, 11]}
    zonal_seasonal = {}
    for name, m_list in seasons.items():
        s_mcs = ds_all.pr_mcssum.sel(month=m_list).sum(dim=['month', 'year'])
        s_tot = ds_all.pr_totsum.sel(month=m_list).sum(dim=['month', 'year'])
        zonal_seasonal[name] = 100 * (s_mcs / s_tot).mean('lon')

    # --- 2. Figure Setup ---
    fig = plt.figure(figsize=(12, 4.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[4, 1], wspace=0.08, top=0.82)

    ax = fig.add_subplot(gs[0], projection=ccrs.PlateCarree(180))
    ax.coastlines(linewidth=0.5); ax.add_feature(cfeat.LAND, facecolor='none', edgecolor='grey', linewidth=0.3)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.4, linestyle=':')
    gl.top_labels = False; gl.right_labels = False
    gl.xlabel_style = {'size': 8}; gl.ylabel_style = {'size': 8}

    cp = ax.contourf(ds_all.lon, ds_all.lat, map_data, cmap=mymap, levels=np.arange(10,100,10), transform=ccrs.PlateCarree(),extend='both')
    
    # Overlay 10% Frequency Contour (Annual Mean)
    mcs_freq_pct = 100 * ds_all.mcs_freq.mean(['month', 'year'])
    cf = ax.contour(ds_all.lon, ds_all.lat, mcs_freq_pct, levels=[5,10], colors=['k'], linewidths=1, transform=ccrs.PlateCarree())
    ax.clabel(cf, inline=True, fontsize=7, fmt='%1.0f%%')
    ax.set_ylim([-30, 30])
    ax.set_title(rf'MCS Precip Contribution ({ds_all.year[0].values}-{ds_all.year[-1].values})'
                                         , loc='left', fontweight='bold', fontsize=11, pad=12)

    # Right Panel: Zonal Statistics
    ax_zonal = fig.add_subplot(gs[1])
    fig.canvas.draw()
    pos_map = ax.get_position(); pos_zonal = ax_zonal.get_position()
    ax_zonal.set_position([pos_zonal.x0, pos_map.y0, pos_zonal.width, pos_map.height])
    
    ax_zonal.fill_betweenx(zonal_mean_ann.lat, zonal_mean_ann - zonal_std_ann, 
                           zonal_mean_ann + zonal_std_ann, color='gray', alpha=0.25, label=r'Ann 1$\sigma$ IAV')
    ax_zonal.plot(zonal_mean_ann, zonal_mean_ann.lat, color='black', linewidth=2, label='Annual Mean', zorder=5)

    colors = {'DJF': 'tab:blue', 'MAM': 'tab:green', 'JJA': 'tab:red', 'SON': 'tab:orange'}
    for sea in ['DJF', 'MAM', 'JJA', 'SON']:
        ax_zonal.plot(zonal_seasonal[sea], zonal_seasonal[sea].lat, color=colors[sea], linewidth=1.2, label=sea)
    
    ax_zonal.set_ylim([-30, 30]); ax_zonal.set_xlim([0, 70]); ax_zonal.set_xlabel('(%)', fontsize=8)
    ax_zonal.set_yticklabels([]); ax_zonal.grid(True, linestyle=':', alpha=0.5)
    ax_zonal.legend(fontsize=7, loc='upper left', frameon=False, bbox_to_anchor=(1.02, 1))
    ax_zonal.set_title('Zonal Mean & Seasons', fontsize=9, fontweight='bold')
    ax_zonal.tick_params(labelsize=8)

    # Colorbar
    cbar_ax = fig.add_axes([pos_map.x1 - 0.25, pos_map.y1 + 0.04, 0.25, 0.025])
    cbar = plt.colorbar(cp, cax=cbar_ax, orientation='horizontal', extend='both')
    cbar.ax.tick_params(width=1, length=1, labelsize=6)
    cbar_ax.text(-0.05, 0.5, 'Contr (%)', transform=cbar_ax.transAxes, fontsize=8, va='center', ha='right')

    return fig, ax, ax_zonal

def combine_model_obs_pngs(model_path, obs_path, model_info, obs_id, output_name='combined_plot.png', dpi=1000):
    """
    Combines two PNG images into a vertical 2-panel figure.
    """
    # Load images using PIL to handle various color modes/transparency
    img_model = Image.open(model_path)
    img_obs = Image.open(obs_path)

    # Create the figure
    # We use a taller figsize for the 2x1 vertical stack
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 3.5), layout='compressed')

    # Display Model Image
    ax1.imshow(img_model)
    ax1.set_title(model_info, fontsize=8, fontweight='bold')
    ax1.axis('off')

    # Display Observation Image
    ax2.imshow(img_obs)
    ax2.set_title(obs_id, fontsize=8, fontweight='bold')
    ax2.axis('off')

    # Adjust layout to minimize whitespace
    #plt.tight_layout()
    
    # Save the result
    fig.savefig(output_name, dpi=dpi, bbox_inches='tight')
    print(f"Combined figure saved as: {output_name}")
    
    return fig

def main(ds_file, model_id):

    plt.ioff() # turn off interactive plotting 

    colors2 = plt.cm.jet(np.linspace(0.3, 0.9, 8))
    colors1 = plt.cm.binary_r(np.linspace(1, 0.8, 2))
    colors = np.vstack((colors1, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
 
    main_dir = Path(os.getenv('CODE_ROOT', '.'))
    pod_dir = main_dir / 'diagnostics/MCS_precip_buoy_stats'
    settings_path = pod_dir / 'settings.jsonc'
    with open(settings_path, 'r') as f:
        runjob_settings = json5.load(f)
    opts = runjob_settings['pod_options']
    # read regridding resolution from settings.jsonc, default is using 1-degree 
    regrid_res = opts.get('regrid_res', 1) 

    work_dir = Path('/scratch/wmtsai/mdtf_miniforge/wkdir/MDTF_output/MCS_precip_buoy_stats')
#    work_dir = Path(os.environ["WORK_DIR"])
    parent_dir = work_dir.parent
    mcsmask_dir = work_dir / f'model/netCDF/MCS_identifiers'    
    stats_dir = work_dir / f'model/netCDF/stats'
    fig_dir = work_dir / 'fig'
    obs_dir = Path(os.environ["CODE_ROOT"]) / 'diagnostics/MCS_precip_buoy_stats/mcs_utils/obs_ref'
    os.makedirs(str(fig_dir), exist_ok=True)

    parent_dir = work_dir.parent
    # Get case name and model resolution info 
    print(ds_file)
    print(model_id)

    if (ds_file is None) and (model_id is None):
        df = pd.read_csv(parent_dir / "MDTF_postprocessed_data.csv")
        model_id = df["source_id"].iloc[0]
        ds_file = df["path"].iloc[0]
        lon_res = abs(np.diff(xr.open_dataset(ds_file).lon)[0])
        lat_res = abs(np.diff(xr.open_dataset(ds_file).lat)[0])
        model_info = (model_id, np.round(lon_res, 2), np.round(lat_res, 2))

    elif (ds_file is not None) and (model_id is not None):
        try:
            lon_res = abs(np.diff(ds_file.lon)[0])
            lat_res = abs(np.diff(ds_file.lat)[0])
        except: # lazy fix 
            lon_res = abs(np.diff(ds_file.longitude)[0])
            lat_res = abs(np.diff(ds_file.latitude)[0])           
        model_info = (model_id, np.round(lon_res, 2), np.round(lat_res, 2))

    else:
        raise ValueError("Provide both ds_file and model_id, or neither (both None).")    

    #################################################
    # 2. plotting POD figures 
    #################################################    
    files = sorted(list(stats_dir.glob(f'precip_mcsstats_regridded_monthly.*{regrid_res}deg.nc')))
    # load multi-year stats files 
    if len(files) >= 1:
        ds_merged = xr.open_mfdataset(files, combine='nested', concat_dim='year')
    else:
        raise IOError("An I/O error occurred: No 'precip_mcsstats' files found under /wkdir/stats/")

    # 1. Annual precipitation mean
    (fig, ax, ax_zonal) = plot_prec_clim_seasonal_variability(ds_merged)
    fig.savefig(fig_dir / 'precipitation_climate.png', transparent=False, bbox_inches='tight', dpi=1000)
    combine_model_obs_pngs(fig_dir / 'precipitation_climate.png', obs_dir / 'precipitation_climate.png', 
                           model_info=rf'{model_info[0]} ($\Delta$x:{model_info[1]}$^o$, $\Delta$y:{model_info[1]}$^o$)'+rf', Regrid: {regrid_res}$^o$',
                           obs_id= r'GPM-IMERG V06 Final ($\Delta$x:0.25$^{o}$, $\Delta$y:0.25$^{o}$)'+rf' Regrid: {regrid_res}$^o$',
                           output_name= fig_dir / 'precipitation_climate.MODvsOBS.png', dpi=1000)
    print('Combined figure saved as: ' + str(fig_dir / 'precipitation_climate.MODvsOBS.png'))
  
    # 2. MCS frequency
    (fig, ax, ax_zonal) = plot_mcsfreq_clim_seasonal(ds_merged)
    fig.savefig(fig_dir / 'MCS_frequency.png', transparent=False, bbox_inches='tight', dpi=1000)
    combine_model_obs_pngs(fig_dir / 'MCS_frequency.png', obs_dir / 'MCS_frequency.png',
                           model_info=rf'{model_info[0]} ($\Delta$x:{model_info[1]}$^o$, $\Delta$y:{model_info[1]}$^o$)'+rf', Regrid: {regrid_res}$^o$',
                           obs_id= r'GPM-IMERG V06 Final ($\Delta$x:0.25$^{o}$, $\Delta$y:0.25$^{o}$)'+rf' Regrid: {regrid_res}$^o$',
                           output_name= fig_dir / 'MCS_frequency.MODvsOBS.png', dpi=1000)
    print('Combined figure saved as: ' + str(fig_dir / 'MCS_frequency.MODvsOBS.png'))

    # 3. MCS precipitation contribution
    (fig, ax, ax_zonal) = plot_mcsprecip_contr_seasonal(ds_merged)
    fig.savefig(fig_dir / 'precipitation_MCScontribution.png', transparent=False, bbox_inches='tight', dpi=1000)
    combine_model_obs_pngs(fig_dir / 'precipitation_MCScontribution.png', obs_dir / 'precipitation_MCScontribution.png',
                           model_info=rf'{model_info[0]} ($\Delta$x:{model_info[1]}$^o$, $\Delta$y:{model_info[1]}$^o$)'+rf', Regrid: {regrid_res}$^o$',
                           obs_id= r'GPM-IMERG V06 Final ($\Delta$x:0.25$^{o}$, $\Delta$y:0.25$^{o}$)'+rf' Regrid: {regrid_res}$^o$',
                           output_name= fig_dir / 'precipitation_MCScontribution.MODvsOBS.png', dpi=1000)
    print('Combined figure saved as: ' + str(fig_dir / 'precipitation_MCScontribution.MODvsOBS.png'))

    del ds_merged
    plt.close()

if __name__ == "__main__":

    # call by mdtf
    main(ds_file=None, model_id=None)
     
