#!/usr/bin/env python
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
import datetime as dt

##################################################################

def read_data_box(latbox,lonbox):
    """
    Read in the temperature, humidity, and geopotential data (+ dimensions) from a 3D file at a single time, over a given region
    e.g., tabs,qv,gz,plev,lat,lon = read_data(latbox,lonbox)
    """    
    #Read in dimensions
    plev = dst['level'] #hPa
    lat = dst['latitude'].sel(latitude=latbox)
    lon = dst['longitude'].sel(longitude=lonbox)
    
    #Read in variables
    tabs = dst['T'].sel(latitude=latbox,longitude=lonbox,time=targettime) #temperature (K)
    qv = dsq['Q'].sel(latitude=latbox,longitude=lonbox,time=targettime) #specific humidity (kg/kg)
    gz = dsz['Z'].sel(latitude=latbox,longitude=lonbox,time=targettime) #geopotential (m^2/s^2)
    dst.close()
    dsq.close()
    dsz.close()
    
    return tabs,qv,gz,plev,lat,lon

##################################################################

def read_data():
    """
    Read in the temperature, humidity, and geopotential data (+ dimensions) from a 3D file at a single time, over the whole domain
    e.g., tabs,qv,gz,plev,lat,lon = read_data
    """   
    #Read in dimensions
    plev = dst['level'] #hPa
    lat = dst['latitude']
    lon = dst['longitude']
    
    #Read in variables
    tabs = dst['T'].sel(time=targettime) #temperature (K)
    qv = dsq['Q'].sel(time=targettime) #specific humidity (kg/kg)
    gz = dsz['Z'].sel(time=targettime) #geopotential (m^2/s^2)
    dst.close()
    dsq.close()
    dsz.close()
    
    return tabs,qv,gz,plev,lat,lon

##################################################################

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
    #indexing is from top to bottom so do this way for positive dp. Note that 1 hPa is the top of ERA-5 so can integrate over index 0 to ibot.
    sz = mse.shape
    dp = np.diff(plev.sel(level=slice(1,pbot)))
    dptile = np.transpose(np.tile(dp,(mse.shape[1],mse.shape[2],1)),(2, 0, 1)) #make plev x lat x lon
    dptile = dptile*100 #convert to Pa

    #Do vertical integral
    ibot = np.int(np.where(plev==pbot)[0]) #find index corresponding to pbot
    h = np.sum(mse[0:ibot,:,:]*dptile,axis=0)/g #sum over zeroth (plev) dimension
    
    return h
    
##################################################################

def write_to_file(h,filedir,year,mm,dd,hour):
    """
    Write column-integrated moist static energy (already defined as a data array) to a netcdf file
    """
    h.attrs['units']='J/m^2'
    h.attrs['long_name']='column-integrated moist static energy'
    h.attrs['_FillValue']=-9999
    h.attrs['GridType']='Latitude/Longitude Grid'

    hds = xr.Dataset({'h':h}, attrs={'note':'column integral from 950 hPa'})

    hds.to_netcdf(filedir+'era5.h.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+hour+'.nc', format='NETCDF4')


#Declare subdomain
latbox = slice(25,15)
lonbox = slice(260,270)

#What level to integrate from (hPa)
pbot = 950 


##################################################################


#ds = xr.open_dataset('/gpfs/fs1/collections/rda/data/ds633.0/e5.oper.an.pl/198001/e5.oper.an.pl.128_130_t.ll025sc.1980012900_1980012923.nc')
#ds = xr.open_dataset('/gpfs/fs1/collections/rda/data/ds633.0/e5.oper.an.pl/198010/e5.oper.an.pl.128_130_t.ll025sc.1980102500_1980102523.nc')
#targettime = str(year)+'-'+"{0:0=2d}".format(mm)+'-'+"{0:0=2d}".format(dd)+'T'+hour[hh]+':00:00.000000000'
#tabs = ds['T'].sel(latitude=latbox,longitude=lonbox,time=targettime)

##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
# ## Read 3D data

hour = ['00','06','12','18']
filebase = '/gpfs/fs1/collections/rda/data/ds633.0/e5.oper.an.pl/' #where to read files from
filedir = '/glade/p/univ/ufsu0014/era5/2D/h/' #where to output files to
years = np.arange(1982,1983,1) #define range of years
for yy,year in enumerate(years): #loop over years
    for mm in range(1,13): #month 1 through 12
        if mm==9 or mm==4 or mm==6 or mm==11: #if September, April, June, or November
            for dd in range(1,31): #30 days in month
                for hh in range(len(hour)):
                    print(mm,dd,hour[hh])
                    #Open dataset
                    dst = xr.open_dataset(filebase+str(year)+"{0:0=2d}".format(mm)+'/e5.oper.an.pl.128_130_t.ll025sc.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'00_'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'23.nc')
                    dsz = xr.open_dataset(filebase+str(year)+"{0:0=2d}".format(mm)+'/e5.oper.an.pl.128_129_z.ll025sc.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'00_'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'23.nc')
                    dsq = xr.open_dataset(filebase+str(year)+"{0:0=2d}".format(mm)+'/e5.oper.an.pl.128_133_q.ll025sc.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'00_'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'23.nc')
                    #Format the time we want
                    targettime = str(year)+'-'+"{0:0=2d}".format(mm)+'-'+"{0:0=2d}".format(dd)+'T'+hour[hh]+':00:00.000000000'
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
                    dst = xr.open_dataset(filebase+str(year)+"{0:0=2d}".format(mm)+'/e5.oper.an.pl.128_130_t.ll025sc.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'00_'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'23.nc')
                    dsz = xr.open_dataset(filebase+str(year)+"{0:0=2d}".format(mm)+'/e5.oper.an.pl.128_129_z.ll025sc.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'00_'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'23.nc')
                    dsq = xr.open_dataset(filebase+str(year)+"{0:0=2d}".format(mm)+'/e5.oper.an.pl.128_133_q.ll025sc.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'00_'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'23.nc')
                    #Format the time we want
                    targettime = str(year)+'-'+"{0:0=2d}".format(mm)+'-'+"{0:0=2d}".format(dd)+'T'+hour[hh]+':00:00.000000000'
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
                    dst = xr.open_dataset(filebase+str(year)+"{0:0=2d}".format(mm)+'/e5.oper.an.pl.128_130_t.ll025sc.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'00_'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'23.nc')
                    dsz = xr.open_dataset(filebase+str(year)+"{0:0=2d}".format(mm)+'/e5.oper.an.pl.128_129_z.ll025sc.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'00_'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'23.nc')
                    dsq = xr.open_dataset(filebase+str(year)+"{0:0=2d}".format(mm)+'/e5.oper.an.pl.128_133_q.ll025sc.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'00_'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'23.nc')
                    #Format the time we want
                    targettime = str(year)+'-'+"{0:0=2d}".format(mm)+'-'+"{0:0=2d}".format(dd)+'T'+hour[hh]+':00:00.000000000'
                    #Read in data
                    #tabs,qv,gz,plev,lat,lon = read_data_box(latbox,lonbox)
                    tabs,qv,gz,plev,lat,lon = read_data()
                    #Compute column-integrated moist static energy
                    h = compute_mse(tabs,qv,gz,plev,pbot)
                    #Write h out to netcdf file
                    write_to_file(h,filedir,year,mm,dd,hour[hh])
#        elif mm==10:
#            for dd in range(30,32): #31 days (1,32)
#                for hh in range(len(hour)):
#                    print(mm,dd,hour[hh])
#                    #Open dataset
#                    dst = xr.open_dataset(filebase+str(year)+"{0:0=2d}".format(mm)+'/e5.oper.an.pl.128_130_t.ll025sc.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'00_'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'23.nc')
#                    dsz = xr.open_dataset(filebase+str(year)+"{0:0=2d}".format(mm)+'/e5.oper.an.pl.128_129_z.ll025sc.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'00_'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'23.nc')
#                    dsq = xr.open_dataset(filebase+str(year)+"{0:0=2d}".format(mm)+'/e5.oper.an.pl.128_133_q.ll025sc.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'00_'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'23.nc')
#                    #Format the time we want
#                    targettime = str(year)+'-'+"{0:0=2d}".format(mm)+'-'+"{0:0=2d}".format(dd)+'T'+hour[hh]+':00:00.000000000'
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
                    dst = xr.open_dataset(filebase+str(year)+"{0:0=2d}".format(mm)+'/e5.oper.an.pl.128_130_t.ll025sc.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'00_'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'23.nc')
                    dsz = xr.open_dataset(filebase+str(year)+"{0:0=2d}".format(mm)+'/e5.oper.an.pl.128_129_z.ll025sc.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'00_'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'23.nc')
                    dsq = xr.open_dataset(filebase+str(year)+"{0:0=2d}".format(mm)+'/e5.oper.an.pl.128_133_q.ll025sc.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'00_'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'23.nc')
                    #Format the time we want
                    targettime = str(year)+'-'+"{0:0=2d}".format(mm)+'-'+"{0:0=2d}".format(dd)+'T'+hour[hh]+':00:00.000000000'
                    #Read in data
                    #tabs,qv,gz,plev,lat,lon = read_data_box(latbox,lonbox)
                    tabs,qv,gz,plev,lat,lon = read_data()
                    #Compute column-integrated moist static energy
                    h = compute_mse(tabs,qv,gz,plev,pbot)
                    #Write h out to netcdf file
                    write_to_file(h,filedir,year,mm,dd,hour[hh])

