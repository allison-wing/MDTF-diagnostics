#!/usr/bin/env python
# coding: utf-8

# # Code to calculate column-integrated moist static energy from 3D data and write out to 2D netcdf file

# In[1]:


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


# In[13]:


def read_data_box(latbox,lonbox):
    """
    Read in the temperature, humidity, and geopotential data (+ dimensions) from a 3D file at a single time, over a given region
    e.g., tabs,qv,gz,plev,lat,lon = read_data(latbox,lonbox)
    """
    #Read in dimensions
    plev = ds['lv_ISBL1']
    lat = ds['g0_lat_2'].sel(g0_lat_2=latbox)
    lon = ds['g0_lon_3'].sel(g0_lon_3=lonbox)
    
    #Read in variables
    tabs = ds['TMP_GDS0_ISBL'].sel(g0_lat_2=latbox,g0_lon_3=lonbox,initial_time0_hours=targettime) #temperature (K)
    qv = ds['SPFH_GDS0_ISBL'].sel(g0_lat_2=latbox,g0_lon_3=lonbox,initial_time0_hours=targettime) #specific humidity (kg/kg)
    z = ds['HGT_GDS0_ISBL'].sel(g0_lat_2=latbox,g0_lon_3=lonbox,initial_time0_hours=targettime) #geopotential height (m)
    ds.close()
    
    #qv only goes up to 100 hPa (dimension lv_ISBL4) so pad with zeros up to 1 hPa, to match size of other vars
    qv=qv.pad(lv_ISBL4=(10,0), constant_values=0) #pad beginning with 10 zeros, this pads lv_ISBL4 with NaNs
    qv=qv.assign_coords(lv_ISBL1=("lv_ISBL4", plev)) #assign new coordinate bsed on plev (up to 1 hPa)
    qv=qv.swap_dims({"lv_ISBL4": "lv_ISBL1"}) #swap which is the dimension 
    qv=qv.drop_vars('lv_ISBL4') #delete old one
    
    return tabs,qv,z,plev,lat,lon


# In[14]:


def read_data():
    """
    Read in the temperature, humidity, and geopotential data (+ dimensions) from a 3D file at a single time, over the whole domain
    e.g., tabs,qv,gz,plev,lat,lon = read_data
    """
    #Read in dimensions
    plev = ds['lv_ISBL1']
    lat = ds['g0_lat_2']
    lon = ds['g0_lon_3']
    
    #Read in variables
    tabs = ds['TMP_GDS0_ISBL'].sel(initial_time0_hours=targettime) #temperature (K)
    qv = ds['SPFH_GDS0_ISBL'].sel(initial_time0_hours=targettime) #specific humidity (kg/kg)
    z = ds['HGT_GDS0_ISBL'].sel(initial_time0_hours=targettime) #geopotential height (m)
    ds.close()
    
    #qv only goes up to 100 hPa (dimension lv_ISBL4) so pad with zeros up to 1 hPa, to match size of other vars
    qv=qv.pad(lv_ISBL4=(10,0), constant_values=0) #pad beginning with 10 zeros, this pads lv_ISBL4 with NaNs
    qv=qv.assign_coords(lv_ISBL1=("lv_ISBL4", plev)) #assign new coordinate bsed on plev (up to 1 hPa)
    qv=qv.swap_dims({"lv_ISBL4": "lv_ISBL1"}) #swap which is the dimension 
    qv=qv.drop_vars('lv_ISBL4') #delete old one
    
    return tabs,qv,z,plev,lat,lon


# In[4]:


def compute_mse(tabs,qv,z,plev,pbot):
    """
    Compute column-integrated moist static energy
    e.g., h = compute_mse(tabs,qv,gz,plev,pbot)
    """
    #Define constants
    cp = 1.00464e3
    g=9.8
    Lv=2.501e6

    #Compute moist static energy
    mse = cp*tabs + g*z + Lv*qv

    #Select the range we are integrating over and define dp. This selects betweeh 1 hPa and the bottom pressure defined as pbot
    #indexing is from top to bottom so do this way for positive dp. Note that 1 hPa is the top of JRA-55
    sz = mse.shape
    dp = np.diff(plev.sel(lv_ISBL1=slice(1,pbot)))
    dptile = np.transpose(np.tile(dp,(mse.shape[1],mse.shape[2],1)),(2, 0, 1)) #make plev x lat x lon
    dptile = dptile*100 #convert to Pa

    #Do vertical integral
    ibot = np.int(np.where(plev==pbot)[0]) #find index corresponding to pbot
    h = np.sum(mse[0:ibot,:,:]*dptile,axis=0)/g #sum over zeroth (plev) dimension
    
    return h


# In[5]:


def write_to_file(h,filedir,year,mm,dd,hour):
    """
    Write column-integrated moist static energy (already defined as a data array) to a netcdf file
    """
    h.attrs['units']='J/m^2'
    h.attrs['long_name']='column-integrated moist static energy'
    h.attrs['_FillValue']=-9999
    h.attrs['GridType']='Latitude/Longitude Grid'
    h=h.rename({'g0_lat_2':'latitude','g0_lon_3':'longitude'})

    hds = xr.Dataset({'h':h}, attrs={'note':'column integral from 950 hPa'})

    hds.to_netcdf(filedir+'jra55.h.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+hour+'.nc', format='NETCDF4')


# In[6]:


#Declare subdomain
latbox = slice(25,15)
lonbox = slice(260,270)

#What level to integrate from (hPa)
pbot = 950 


# In[7]:


ds = xr.open_dataset('/huracan/tank2/columbia/reanalysis/jra55/3D/1990/anl_p125.199004.grb.wing',engine='pynio')
ds


# ## Read 3D data

# In[15]:


hour = ['00','06','12','18']
filebase = '/huracan/tank2/columbia/reanalysis/jra55/3D/' #where to read files from
filedir = '/huracan/tank2/columbia/reanalysis/jra55/2D/h/' #where to output files to
years = np.arange(2018,2020,1) #define range of years
for yy,year in enumerate(years): #loop over years
    for mm in range(1,13): #month 1 through 12
        if mm==9 or mm==4 or mm==6 or mm==11: #if September, April, June, or November
            for dd in range(1,31): #30 days in month
                for hh in range(len(hour)):
                    print(mm,dd,hour[hh])
                    #Open dataset
                    ds = xr.open_dataset(filebase+str(year)+'/anl_p125.'+str(year)+"{0:0=2d}".format(mm)+'.grb.wing',engine='pynio')
                    #Format the time we want
                    targettime = str(year)+'-'+"{0:0=2d}".format(mm)+'-'+"{0:0=2d}".format(dd)+'T'+hour[hh]+':00:00.000000000'                    
                    #Read in data   
                    #tabs,qv,z,plev,lat,lon = read_data_box(latbox,lonbox)
                    tabs,qv,z,plev,lat,lon = read_data()
                    #Compute column-integrated moist static energy
                    h = compute_mse(tabs,qv,z,plev,pbot)
                    #Write h out to netcdf file
                    write_to_file(h,filedir,year,mm,dd,hour[hh])               
        elif mm==2 and year % 4 ==0: #leap year February
            for dd in range(1,30): #29 days
                for hh in range(len(hour)):
                    print(mm,dd,hour[hh])
                    #Open dataset
                    ds = xr.open_dataset(filebase+str(year)+'/anl_p125.'+str(year)+"{0:0=2d}".format(mm)+'.grb.wing',engine='pynio')
                    #Format the time we want
                    targettime = str(year)+'-'+"{0:0=2d}".format(mm)+'-'+"{0:0=2d}".format(dd)+'T'+hour[hh]+':00:00.000000000'
                    #Read in data
                    #tabs,qv,z,plev,lat,lon = read_data_box(latbox,lonbox)
                    tabs,qv,z,plev,lat,lon = read_data()
                    #Compute column-integrated moist static energy
                    h = compute_mse(tabs,qv,z,plev,pbot)
                    #Write h out to netcdf file
                    write_to_file(h,filedir,year,mm,dd,hour[hh])
        elif mm==2 and year % 4 !=0: #non-leap year February
            for dd in range(1,29): #28 days
                for hh in range(len(hour)):
                    print(mm,dd,hour[hh])
                    #Open dataset
                    ds = xr.open_dataset(filebase+str(year)+'/anl_p125.'+str(year)+"{0:0=2d}".format(mm)+'.grb.wing',engine='pynio')
                    #Format the time we want
                    targettime = str(year)+'-'+"{0:0=2d}".format(mm)+'-'+"{0:0=2d}".format(dd)+'T'+hour[hh]+':00:00.000000000'
                    #Read in data
                    #tabs,qv,z,plev,lat,lon = read_data_box(latbox,lonbox)
                    tabs,qv,z,plev,lat,lon = read_data()
                    #Compute column-integrated moist static energy
                    h = compute_mse(tabs,qv,z,plev,pbot)
                    #Write h out to netcdf file
                    write_to_file(h,filedir,year,mm,dd,hour[hh])
#        elif mm==10:
#            for dd in range(20,32): #31 days (1,32)
#                for hh in range(len(hour)):
#                    print(mm,dd,hour[hh])
#                    #Open dataset
#                    ds = xr.open_dataset(filebase+str(year)+'/anl_p125.'+str(year)+"{0:0=2d}".format(mm)+'.grb.wing',engine='pynio')
#                    #Format the time we want
#                    targettime = str(year)+'-'+"{0:0=2d}".format(mm)+'-'+"{0:0=2d}".format(dd)+'T'+hour[hh]+':00:00.000000000'
#                    #Read in data
#                    #tabs,qv,z,plev,lat,lon = read_data_box(latbox,lonbox)
#                    tabs,qv,z,plev,lat,lon = read_data()
#                    #Compute column-integrated moist static energy
#                    h = compute_mse(tabs,qv,z,plev,pbot)
#                    #Write h out to netcdf file
#                    write_to_file(h,filedir,year,mm,dd,hour[hh])
        else:
            for dd in range(1,32): #31 days
                for hh in range(len(hour)):
                    print(mm,dd,hour[hh])
                    #Open dataset
                    ds = xr.open_dataset(filebase+str(year)+'/anl_p125.'+str(year)+"{0:0=2d}".format(mm)+'.grb.wing',engine='pynio')
                    #Format the time we want
                    targettime = str(year)+'-'+"{0:0=2d}".format(mm)+'-'+"{0:0=2d}".format(dd)+'T'+hour[hh]+':00:00.000000000'
                    #Read in data
                    #tabs,qv,z,plev,lat,lon = read_data_box(latbox,lonbox)
                    tabs,qv,z,plev,lat,lon = read_data()
                    #Compute column-integrated moist static energy
                    h = compute_mse(tabs,qv,z,plev,pbot)
                    #Write h out to netcdf file
                    write_to_file(h,filedir,year,mm,dd,hour[hh])


