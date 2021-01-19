#!/usr/bin/env python
# coding: utf-8

# # Code to calculate column-integrated moist static energy from 3D data and write out to 2D netcdf file

# In[2]:


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


# In[3]:


def read_data_box(latbox,lonbox):
    """
    Read in the temperature, humidity, and geopotential data (+ dimensions) from a 3D file at a single time, over a given region
    e.g., tabs,qv,gz,plev,lat,lon = read_data(latbox,lonbox)
    """
    #Read in dimensions
    plev = ds['lv_ISBL0']#pressure (Pa)
    lat = ds['lat_0'].sel(lat_0=latbox)
    lon = ds['lon_0'].sel(lon_0=lonbox)
    
    #Read in variables
    tabs = ds['TMP_P0_L100_GLL0'].sel(lat_0=latbox,lon_0=lonbox) #temperature (K)
    qv = ds['SPFH_P0_L100_GLL0'].sel(lat_0=latbox,lon_0=lonbox) #specific humidity (kg/kg)
    z = ds['HGT_P0_L100_GLL0'].sel(lat_0=latbox,lon_0=lonbox) #geopotential height (m)
    ds.close()
    
    return tabs,qv,z,plev,lat,lon


# In[4]:


def read_data():
    """
    Read in the temperature, humidity, and geopotential data (+ dimensions) from a 3D file at a single time, over the whole domain
    e.g., tabs,qv,gz,plev,lat,lon = read_data
    """
    #Read in dimensions
    plev = ds['lv_ISBL0']#pressure (Pa)
    lat = ds['lat_0']
    lon = ds['lon_0']
    
    #Read in variables
    tabs = ds['TMP_P0_L100_GLL0'] #temperature (K)
    qv = ds['SPFH_P0_L100_GLL0'] #specific humidity (kg/kg)
    z = ds['HGT_P0_L100_GLL0'] #geopotential height (m)
    ds.close()
    
    return tabs,qv,z,plev,lat,lon


# In[5]:


def compute_mse(tabs,qv,z,plev,pbot):
    """
    Compute column-integrated moist static energy
    e.g., h = compute_mse(tabs,qv,z,plev,pbot) [plev in Pa, pbot in hPa]
    """
    #Define constants
    cp = 1.00464e3
    g=9.8
    Lv=2.501e6

    #Compute moist static energy
    mse = cp*tabs + g*z + Lv*qv

    #Select the range we are integrating over and define dp. This selects betweeh 1 hPa and the bottom pressure defined as pbot
    #indexing is from top to bottom so do this way for positive dp
    sz = mse.shape
    dp = np.diff(plev.sel(lv_ISBL0=slice(1*100,pbot*100)))
    dptile = np.transpose(np.tile(dp,(sz[1],sz[2],1)),(2, 0, 1)) #make plev x lat x lon
    #dptile = dptile*100 #convert to Pa (already in Pa)

    #Do vertical integral
    ibot = np.int(np.where(plev==pbot*100)[0]) #find index corresponding to pbot
    h = np.sum(mse[0:ibot,:,:]*dptile,axis=0)/g #sum over zeroth (plev) dimension
    
    return h


# In[6]:


def write_to_file(h,filedir,year,mm,dd,hour):
    """
    Write column-integrated moist static energy (already defined as a data array) to a netcdf file
    """
    h.attrs['units']='J/m^2'
    h.attrs['long_name']='column-integrated moist static energy'
    h.attrs['_FillValue']=-9999
    h.attrs['GridType']='Latitude/Longitude Grid'
    h=h.rename({'lat_0':'latitude','lon_0':'longitude'})

    hds = xr.Dataset({'h':h}, attrs={'note':'column integral from 950 hPa to 1hPa'})

    hds.to_netcdf(filedir+'cfsr.h.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+hour+'.nc', format='NETCDF4')


# In[7]:


#Declare subdomain
latbox = slice(25,15)
lonbox = slice(260,270)

#What level to integrate from (hPa)
pbot = 950 


# In[8]:


#ds = xr.open_dataset('/huracan/tank2/columbia/reanalysis/cfsr/3D/1980/198001/pgbh00.gdas.1980013118.grb2.wing',engine='pynio')
#tabs = ds['TMP_P0_L100_GLL0']


# ## Read 3D data

# In[9]:


hour = ['00','06','12','18']
filebase = '/huracan/tank2/columbia/reanalysis/cfsr/3D/' #where to read files from
filedir = '/huracan/tank2/columbia/reanalysis/cfsr/2D/h/' #where to output files to
years = np.arange(2019,2020,1) #define range of years
for yy,year in enumerate(years): #loop over years
    for mm in range(1,13): #month 1 through 12
        if mm==9 or mm==4 or mm==6 or mm==11: #if September, April, June, or November
            for dd in range(1,31): #30 days in month
                for hh in range(len(hour)):
                    print(mm,dd,hour[hh])
                    #Open dataset
                    if year<2011:
                        ds = xr.open_dataset(filebase+str(year)+'/'+str(year)+"{0:0=2d}".format(mm)+'/pgbh00.gdas.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+hour[hh]+'.grb2.wing',engine='pynio')
                    elif year==2011 and mm<4:
                        ds = xr.open_dataset(filebase+str(year)+'/'+str(year)+"{0:0=2d}".format(mm)+'/pgbh00.gdas.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+hour[hh]+'.grb2.wing',engine='pynio')
                    else:
                        ds = xr.open_dataset(filebase+str(year)+'/'+str(year)+"{0:0=2d}".format(mm)+'/cdas1.pgrbh00.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+hour[hh]+'.grib2.wing',engine='pynio')                        
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
                    if year<2011:
                        ds = xr.open_dataset(filebase+str(year)+'/'+str(year)+"{0:0=2d}".format(mm)+'/pgbh00.gdas.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+hour[hh]+'.grb2.wing',engine='pynio')
                    elif year==2011 and mm<4:
                        ds = xr.open_dataset(filebase+str(year)+'/'+str(year)+"{0:0=2d}".format(mm)+'/pgbh00.gdas.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+hour[hh]+'.grb2.wing',engine='pynio')
                    else:
                        ds = xr.open_dataset(filebase+str(year)+'/'+str(year)+"{0:0=2d}".format(mm)+'/cdas1.pgrbh00.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+hour[hh]+'.grib2.wing',engine='pynio')                    #Read in data
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
                    if year<2011:
                        ds = xr.open_dataset(filebase+str(year)+'/'+str(year)+"{0:0=2d}".format(mm)+'/pgbh00.gdas.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+hour[hh]+'.grb2.wing',engine='pynio')
                    elif year==2011 and mm<4:
                        ds = xr.open_dataset(filebase+str(year)+'/'+str(year)+"{0:0=2d}".format(mm)+'/pgbh00.gdas.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+hour[hh]+'.grb2.wing',engine='pynio')
                    else:
                        ds = xr.open_dataset(filebase+str(year)+'/'+str(year)+"{0:0=2d}".format(mm)+'/cdas1.pgrbh00.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+hour[hh]+'.grib2.wing',engine='pynio')                  #Read in data
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
#                    if year<2011:
#                        ds = xr.open_dataset(filebase+str(year)+'/'+str(year)+"{0:0=2d}".format(mm)+'/pgbh00.gdas.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+hour[hh]+'.grb2.wing',engine='pynio')
#                    elif year==2011 and mm<4:
#                        ds = xr.open_dataset(filebase+str(year)+'/'+str(year)+"{0:0=2d}".format(mm)+'/pgbh00.gdas.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+hour[hh]+'.grb2.wing',engine='pynio')
#                    else:
#                        ds = xr.open_dataset(filebase+str(year)+'/'+str(year)+"{0:0=2d}".format(mm)+'/cdas1.pgrbh00.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+hour[hh]+'.grib2.wing',engine='pynio')                    #Read in data
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
                    if year<2011:
                        ds = xr.open_dataset(filebase+str(year)+'/'+str(year)+"{0:0=2d}".format(mm)+'/pgbh00.gdas.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+hour[hh]+'.grb2.wing',engine='pynio')
                    elif year==2011 and mm<4:
                        ds = xr.open_dataset(filebase+str(year)+'/'+str(year)+"{0:0=2d}".format(mm)+'/pgbh00.gdas.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+hour[hh]+'.grb2.wing',engine='pynio')
                    else:
                        ds = xr.open_dataset(filebase+str(year)+'/'+str(year)+"{0:0=2d}".format(mm)+'/cdas1.pgrbh00.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+hour[hh]+'.grib2.wing',engine='pynio')                #Read in data
                    #tabs,qv,gz,plev,lat,lon = read_data_box(latbox,lonbox)
                    tabs,qv,gz,plev,lat,lon = read_data()
                    #Compute column-integrated moist static energy
                    h = compute_mse(tabs,qv,gz,plev,pbot)
                    #Write h out to netcdf file
                    write_to_file(h,filedir,year,mm,dd,hour[hh])


# In[9]:


#h


# In[10]:


#Plot to check that it looks reasonable-ish
#fig = plt.figure(figsize=(10,10))
#ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
#im=ax.contourf(lon,lat,h,30,transform=ccrs.PlateCarree())

#ax.coastlines()
#ax.coastlines('50m', edgecolor='white')
#gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                  linewidth=2, color='gray', alpha=0.5, linestyle='--')

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
