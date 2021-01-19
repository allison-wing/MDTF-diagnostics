#!/usr/bin/env python
# coding: utf-8

# # Code to calculate column-integrated moist static energy from 3D data and write out to 2D netcdf file

import sys
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy import config
import cartopy.feature as cfeature
from cartopy.vector_transform import vector_scalar_to_grid
from matplotlib.axes import Axes
import scipy as sp
import xarray as xr
import math
import datetime
import Nio

def read_data_box(latbox,lonbox):
    """
    Read in the temperature, humidity, and geopotential data (+ dimensions) from a 3D file at a single time, over a given region
    e.g., tabs,qv,gz,plev,lat,lon = read_data(latbox,lonbox)
    """
    #Read in dimensions
    plev = ds['lv_ISBL0']
    lat = ds['g4_lat_1'].sel(g4_lat_1=latbox)
    lon = ds['g4_lon_2'].sel(g4_lon_2=lonbox)
    
    #Read in variables
    tabs = ds['T_GDS4_ISBL'].sel(g4_lat_1=latbox,g4_lon_2=lonbox) #temperature (K)
    qv = ds['Q_GDS4_ISBL'].sel(g4_lat_1=latbox,g4_lon_2=lonbox) #specific humidity (kg/kg)
    gz = ds['Z_GDS4_ISBL'].sel(g4_lat_1=latbox,g4_lon_2=lonbox) #geopotential (m^2/s^2)
    ds.close()
    
    return tabs,qv,gz,plev,lat,lon

def read_data():
    """
    Read in the temperature, humidity, and geopotential data (+ dimensions) from a 3D file at a single time, over the whole domain
    e.g., tabs,qv,gz,plev,lat,lon = read_data
    """
    #Read in dimensions
    plev = ds['lv_ISBL0']
    lat = ds['g4_lat_1']
    lon = ds['g4_lon_2']
    
    #Read in variables
    tabs = ds['T_GDS4_ISBL'] #temperature (K)
    qv = ds['Q_GDS4_ISBL'] #specific humidity (kg/kg)
    gz = ds['Z_GDS4_ISBL'] #geopotential (m^2/s^2)
    ds.close()
    
    return tabs,qv,gz,plev,lat,lon

def compute_mse(tabs,qv,gz,plev,pbot):
    """
    Compute column-integrated moist static energy
    e.g., h = compute_mse(tabs,qv,gz,plev,pbot)
    """
    #Define constants
    cp = 1.00464e3
    g=9.8
    Lv=2.501e6

    #Compute moist static energy
    mse = cp*tabs + gz + Lv*qv

    #Select the range we are integrating over and define dp. This selects betweeh 1 hPa and the bottom pressure defined as pbot
    #indexing is from top to bottom so do this way for positive dp
    sz = mse.shape
    dp = np.diff(plev.sel(lv_ISBL0=slice(1,pbot)))
    dptile = np.transpose(np.tile(dp,(mse.shape[1],mse.shape[2],1)),(2, 0, 1)) #make plev x lat x lon
    dptile = dptile*100 #convert to Pa

    #Do vertical integral
    ibot = np.int(np.where(plev==pbot)[0]) #find index corresponding to pbot
    h = np.sum(mse[0:ibot,:,:]*dptile,axis=0)/g #sum over zeroth (plev) dimension
    
    return h

def write_to_file(h,filedir,year,mm,dd,hour):
    """
    Write column-integrated moist static energy (already defined as a data array) to a netcdf file
    """
    h.attrs['units']='J/m^2'
    h.attrs['long_name']='column-integrated moist static energy'
    h.attrs['_FillValue']=-9999
    h.attrs['GridType']='Gaussian Latitude/Longitude Grid'
    h=h.rename({'g4_lat_1':'latitude','g4_lon_2':'longitude'})

    hds = xr.Dataset({'h':h}, attrs={'note':'column integral from 950 hPa'})

    hds.to_netcdf(filedir+'erai.h.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+hour+'.nc', format='NETCDF4')

#Declare subdomain
latbox = slice(25,15)
lonbox = slice(260,270)

#What level to integrate from (hPa)
pbot = 950 

# ## Read 3D data

hour = ['00','06','12','18']
filebase = '/huracan/tank2/columbia/reanalysis/erai/3D/' #where to read files from
filedir = '/huracan/tank2/columbia/reanalysis/erai/2D/h/' #where to output files to
years = np.arange(1985,2019,1) #define range of years
for yy,year in enumerate(years): #loop over years
    for mm in range(1,13): #month 1 through 12
        if mm==9 or mm==4 or mm==6 or mm==11: #if September, April, June, or November
            for dd in range(1,31): #30 days in month
                for hh in range(len(hour)):
                    print(mm,dd,hour[hh])
                    #Open dataset
                    ds = xr.open_dataset(filebase+str(year)+'/'+str(year)+"{0:0=2d}".format(mm)+'/ei.oper.an.pl.regn128sc.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+hour[hh]+'.wing',engine='pynio')
                    #Read in data
                    #tabs,qv,gz,plev,lat,lon = read_data_box(latbox,lonbox)
                    tabs,qv,gz,plev,lat,lon = read_data()
                    #Compute column-integrated moist static energy
                    h = compute_mse(tabs,qv,gz,plev,pbot)
                    #Write h out to netcdf file
                    write_to_file(h,filedir,year,mm,dd,hour[hh])               
        elif mm==2 and year % 4 ==0: #leap year February
            for dd in range(1,30): #29 days
                for hh in range(len(hour)):
                    print(mm,dd,hour[hh])
                    #Open dataset
                    ds = xr.open_dataset(filebase+str(year)+'/'+str(year)+"{0:0=2d}".format(mm)+'/ei.oper.an.pl.regn128sc.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+hour[hh]+'.wing',engine='pynio')
                    #Read in data
                    tabs,qv,gz,plev,lat,lon = read_data_box(latbox,lonbox)
                    tabs,qv,gz,plev,lat,lon = read_data()
                    #Compute column-integrated moist static energy
                    h = compute_mse(tabs,qv,gz,plev,pbot)
                    #Write h out to netcdf file
                    write_to_file(h,filedir,year,mm,dd,hour[hh])
        elif mm==2 and year % 4 !=0: #non-leap year February
            for dd in range(1,29): #28 days
                for hh in range(len(hour)):
                    print(mm,dd,hour[hh])
                    #Open dataset
                    ds = xr.open_dataset(filebase+str(year)+'/'+str(year)+"{0:0=2d}".format(mm)+'/ei.oper.an.pl.regn128sc.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+hour[hh]+'.wing',engine='pynio')
                    #Read in data
                    #tabs,qv,gz,plev,lat,lon = read_data_box(latbox,lonbox)
                    tabs,qv,gz,plev,lat,lon = read_data()
                    #Compute column-integrated moist static energy
                    h = compute_mse(tabs,qv,gz,plev,pbot)
                    #Write h out to netcdf file
                    write_to_file(h,filedir,year,mm,dd,hour[hh])
#        elif mm==10:
#            for dd in range(20,32): #31 days (1,32)
#                for hh in range(len(hour)):
#                    print(mm,dd,hour[hh])
#                    #Open dataset
#                    ds = xr.open_dataset(filebase+str(year)+'/'+str(year)+"{0:0=2d}".format(mm)+'/ei.oper.an.pl.regn128sc.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+hour[hh]+'.wing',engine='pynio')
#                    #Read in data
#                    #tabs,qv,gz,plev,lat,lon = read_data_box(latbox,lonbox)
#                    tabs,qv,gz,plev,lat,lon = read_data()
#                    #Compute column-integrated moist static energy
#                    h = compute_mse(tabs,qv,gz,plev,pbot)
#                    #Write h out to netcdf file
#                    write_to_file(h,filedir,year,mm,dd,hour[hh])
        else:
            for dd in range(1,32): #31 days
                for hh in range(len(hour)):
                    print(mm,dd,hour[hh])
                    #Open dataset
                    ds = xr.open_dataset(filebase+str(year)+'/'+str(year)+"{0:0=2d}".format(mm)+'/ei.oper.an.pl.regn128sc.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+hour[hh]+'.wing',engine='pynio')
                    #Read in data
                    #tabs,qv,gz,plev,lat,lon = read_data_box(latbox,lonbox)
                    tabs,qv,gz,plev,lat,lon = read_data()
                    #Compute column-integrated moist static energy
                    h = compute_mse(tabs,qv,gz,plev,pbot)
                    #Write h out to netcdf file
                    write_to_file(h,filedir,year,mm,dd,hour[hh])

##Plot to check that it looks reasonable-ish
#fig = plt.figure(figsize=(10,10))
#ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
#im=ax.contourf(lon,lat,h,30,transform=ccrs.PlateCarree())
#
#ax.coastlines()
#ax.coastlines('50m', edgecolor='white')
#gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
#
#fig.colorbar(im,ax=ax)
#plt.title('Column Moist Static Energy [J/m^2]\n'+str(year)+'-'+"{0:0=2d}".format(mm)+'-'+"{0:0=2d}".format(dd)+'-'+str(hour[hh])+'\n')
#plt.show()


# ## How to compute column-integrated moist static energy

# In matlab, I would do  
# sz=size(mse);  
# ibot = find(plev==pbot);  
# p = plev(1:ibot);  
# dp = diff(p);  
# dp = repmat(dp,sz(2),sz(3),1);  
# dp = permute(dp,[3 1 2]);  
# for i = 1:sz(1)  
#     for j = 1:sz(2)  
#         h(i,j) = sum(mse(i,j,1:ibot-1).* dp(i,j,:),3)/g;  
#     end  
# end
#         
