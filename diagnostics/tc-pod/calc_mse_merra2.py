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
import datetime



def read_data_box(latbox,lonbox):
    """
    Read in the temperature, humidity, and geopotential data (+ dimensions) from a 3D file over all times in that file, over a given region
    e.g., tabs,qv,gz,plev,lat,lon = read_data(latbox,lonbox)
    """
    #Read in dimensions
    plev = ds['lev'] #hPa
    lat = ds['lat'].sel(lat=latbox)
    lon = ds['lon'].sel(lon=lonbox)
    time = ds['time']/60 #convert to hours
    
    #Read in variables
    tabs = ds['T'].sel(lat=latbox,lon=lonbox) #temperature (K)
    tabs = tabs.where(tabs < 1e15) #set missing data (1e15) to NaN
    qv = ds['QV'].sel(lat=latbox,lon=lonbox) #specific humidity (kg/kg)
    qv = qv.where(qv < 1e15) #set missing data (1e15) to NaN
    z = ds['H'].sel(lat=latbox,lon=lonbox) #geopotential height (m)
    z = z.where(z < 1e15) #set missing data (1e15) to NaN
    ds.close()
    
    return tabs,qv,z,plev,lat,lon, time




def read_data():
    """
    Read in the temperature, humidity, and geopotential data (+ dimensions) from a 3D file over all times in that file, over the whole domain
    e.g., tabs,qv,gz,plev,lat,lon = read_data 
    For merrea-2, dimensions are time x lev x lat x lon 
    """
    #Read in dimensions
    plev = ds['lev'] #hPa
    lat = ds['lat']
    lon = ds['lon']
    time = ds['time']/60 #convert to hours
    
    #Read in variables
    tabs = ds['T'] #temperature (K) 
    tabs = tabs.where(tabs < 1e15) #set missing data (1e15) to NaN
    qv = ds['QV'] #specific humidity (kg/kg)
    qv = qv.where(qv < 1e15) #set missing data (1e15) to NaN
    z = ds['H'] #geopotential height (m)
    z = z.where(z < 1e15) #set missing data (1e15) to NaN
    ds.close()
    
    return tabs,qv,z,plev,lat,lon,time




def compute_mse(tabs,qv,z,plev,pbot):
    """
    Compute column-integrated moist static energy
    e.g., h = compute_mse(tabs,qv,z,plev,pbot)
    tabs, qv, and z should be restricted to a given time and be of shape lev x lat x lon
    """
    #Define constants
    cp = 1.00464e3
    g=9.8
    Lv=2.501e6

    #Compute moist static energy
    mse = cp*tabs + g*z + Lv*qv

    #Select the range we are integrating over and define dp. This selects betweeh 1 hPa and the bottom pressure defined as pbot
    #indexing is from bottom to top so do this way for positive dp [note that MERRA2 goes up to 0.1 hPa but I am hard coding 
    #to 1 hPa to be consistent with other datasets]
    sz = mse.shape
    dp = -1*np.diff(plev.sel(lev=slice(pbot,1)))
    dptile = np.transpose(np.tile(dp,(mse.shape[1],mse.shape[2],1)),(2, 0, 1)) #make plev x lat x lon
    dptile = dptile*100 #convert to Pa

    #Do vertical integral
    ibot = np.int(np.where(plev==pbot)[0]) #find index corresponding to pbot
    itop = np.int(np.where(plev==1)[0]) #find index corresponding to ptop (1 hPa)
    h = np.sum(mse[ibot+1:itop+1,:,:]*dptile,axis=0)/g #sum over zeroth (plev) dimension indexngs means to go from ibot+1 to the end               
    
    return h




def write_to_file(h,filedir,year,mm,dd,hour):
    """
    Write column-integrated moist static energy (already defined as a data array) to a netcdf file
    """
    h.attrs['units']='J/m^2'
    h.attrs['long_name']='column-integrated moist static energy'
    h.attrs['_FillValue']=-9999
    h.attrs['GridType']='Latitude/Longitude Grid'
    h=h.rename({'lat':'latitude','lon':'longitude'})

    hds = xr.Dataset({'h':h}, attrs={'note':'column integral from 950 hPa'})

    hds.to_netcdf(filedir+'merra2.h.'+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+"{0:0=2d}".format(hour)+'.nc', format='NETCDF4')




##Declare subdomain
#latbox = slice(5,25)
#lonbox = slice(-90,-40) #merra2 uses -180:180

#What level to integrate from (hPa)
pbot = 950 


# ## Read 3D data


filedir = '/huracan/tank2/columbia/reanalysis/merra2/2D/h/' #where to output files to
filebase = '/huracan/tank2/columbia/reanalysis/merra2/3D/' #where to read files from
years = np.arange(1990,2019,1) #define range of years
for yy,year in enumerate(years): #loop over years
    if year<1992:
        filepart = 'MERRA2_100.inst6_3d_ana_Np.'
    elif year>=1992 and year<2001:
        filepart = 'MERRA2_200.inst6_3d_ana_Np.'
    elif year>=2001 and year<2011:
        filepart = 'MERRA2_300.inst6_3d_ana_Np.'
    else: #after 2011
        filepart = 'MERRA2_400.inst6_3d_ana_Np.'
    for mm in range(1,13): #month 1 through 12
        if mm==9 or mm==4 or mm==6 or mm==11:
            for dd in range(1,31): #30 days in month
                print(mm,dd)
                #Open dataset
                ds = xr.open_dataset(filebase+str(year)+'/'+str(year)+"{0:0=2d}".format(mm)+'/'+filepart+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'.nc4.wing',decode_cf=False)
            	#Read in data
            	#tabs,qv,z,plev,lat,lon,time = read_data_box(latbox,lonbox)
                tabs,qv,z,plev,lat,lon,time = read_data()
                for hh in range(len(time)):
                    print(mm,dd,int(time[hh]))             
                    #Compute column-integrated moist static energy, restricted to that time (feed for specific time)
                    h = compute_mse(tabs[hh,:,:,:],qv[hh,:,:,:],z[hh,:,:,:],plev,pbot)
                    #Write h out to netcdf file
                    write_to_file(h,filedir,year,mm,dd,int(time[hh]))               
        elif mm==2 and int(year) % 4 ==0: #leap year February
            for dd in range(1,30): #29 days
            	print(mm,dd)
            	#Open dataset
            	ds = xr.open_dataset(filebase+str(year)+'/'+str(year)+"{0:0=2d}".format(mm)+'/'+filepart+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'.nc4.wing',decode_cf=False)
            	#Read in data
            	#tabs,qv,z,plev,lat,lon,time = read_data_box(latbox,lonbox)
            	tabs,qv,z,plev,lat,lon,time = read_data()
            	for hh in range(len(time)):
                    print(mm,dd,int(time[hh]))             
                    #Compute column-integrated moist static energy, restricted to that time (feed for specific time)
                    h = compute_mse(tabs[hh,:,:,:],qv[hh,:,:,:],z[hh,:,:,:],plev,pbot)
                    #Write h out to netcdf file
                    write_to_file(h,filedir,year,mm,dd,int(time[hh]))    
        elif mm==2 and int(year) % 4 !=0: #non-leap year February
            for dd in range(1,29): #28 days
            	print(mm,dd)
            	#Open dataset
            	ds = xr.open_dataset(filebase+str(year)+'/'+str(year)+"{0:0=2d}".format(mm)+'/'+filepart+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'.nc4.wing',decode_cf=False)
            	#Read in data
            	#tabs,qv,z,plev,lat,lon,time = read_data_box(latbox,lonbox)
            	tabs,qv,z,plev,lat,lon,time = read_data()
            	for hh in range(len(time)):
                    print(mm,dd,int(time[hh]))             
                    #Compute column-integrated moist static energy, restricted to that time (feed for specific time)
                    h = compute_mse(tabs[hh,:,:,:],qv[hh,:,:,:],z[hh,:,:,:],plev,pbot)
                    #Write h out to netcdf file
                    write_to_file(h,filedir,year,mm,dd,int(time[hh]))    
    	#elif mm==10:
    	#    for dd in range(1,32): #31 days (1,32)
    	#        #Open dataset
    	#        ds = xr.open_dataset(filebase+str(year)+'/'+str(year)+"{0:0=2d}".format(mm)+'/'+filepart+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'.nc4.wing',decode_cf=False)
    	#        #Read in data
    	#        #tabs,qv,z,plev,lat,lon,time = read_data_box(latbox,lonbox)
    	#        tabs,qv,z,plev,lat,lon,time = read_data()
    	#        for hh in range(len(time)):
    	#            print(mm,dd,int(time[hh]))             
    	#            #Compute column-integrated moist static energy, restricted to that time (feed for specific time)
    	#            h = compute_mse(tabs[hh,:,:,:],qv[hh,:,:,:],z[hh,:,:,:],plev,pbot)
    	#            #Write h out to netcdf file
    	#            write_to_file(h,filedir,year,mm,dd,int(time[hh])) 
        else:
            for dd in range(1,32): #31 days
            	print(mm,dd)
            	#Open dataset
            	ds = xr.open_dataset(filebase+str(year)+'/'+str(year)+"{0:0=2d}".format(mm)+'/'+filepart+str(year)+"{0:0=2d}".format(mm)+"{0:0=2d}".format(dd)+'.nc4.wing',decode_cf=False)
            	#Read in data
            	#tabs,qv,z,plev,lat,lon,time = read_data_box(latbox,lonbox)
            	tabs,qv,z,plev,lat,lon,time = read_data()
            	for hh in range(len(time)):
                    print(mm,dd,int(time[hh]))             
                    #Compute column-integrated moist static energy, restricted to that time (feed for specific time)
                    h = compute_mse(tabs[hh,:,:,:],qv[hh,:,:,:],z[hh,:,:,:],plev,pbot)
                    #Write h out to netcdf file
                    write_to_file(h,filedir,year,mm,dd,int(time[hh]))    





#Plot to check that it looks reasonable-ish
#fig = plt.figure(figsize=(10,10))
#ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
#im=ax.contourf(lon,lat,h,30,transform=ccrs.PlateCarree())

#ax.coastlines()
#ax.coastlines('50m', edgecolor='white')
#gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                  linewidth=2, color='gray', alpha=0.5, linestyle='--')

#fig.colorbar(im,ax=ax)
#plt.title('Column Moist Static Energy [J/m^2]\n'+str(year)+'-'+"{0:0=2d}".format(mm)+'-'+"{0:0=2d}".format(dd)+'-'+str(int(time[hh]))+'\n')
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
#         h(i,j) = sum(mse(i,j,ibot+1:end).*dp(i,j,:),3)/g;
#     end  
# end
#         
