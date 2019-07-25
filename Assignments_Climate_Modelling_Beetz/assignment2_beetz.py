import math
import numpy as np
import matplotlib.pyplot as pl
import netCDF4 as nc
from netCDF4 import Dataset
from netCDF4 import num2date
import datetime 
import pandas as pd


# Import of NetCDF Files with every variables
def import_nc(pathf,var):
    file = Dataset(pathf)
    lat = np.array(file.variables['latitude'])
    lon = np.array(file.variables['longitude'])
    data = np.array(file.variables[var])
    time = file.variables['time']
    dates = nc.num2date(time[:], units=time.units)
    return lat,lon,data,dates

lat,lon,data,dates = import_nc('./HadISST_sst.nc','sst') #call the NetCDF Import function

# Dataset correction (deleting NaNs)
data[data<=-1000.]=np.nan  #-1000,-1.0e+30 # where data is value x put in NaNs instead

lat1 = 5        # 5°N
lat2 = -5       # 5°S
lon1 = -170     # 170°W
lon2 = -120     # 120°W

data_area = data[:,(lat<=lat1) & (lat>=lat2),:][:,:,(lon>=lon1) & (lon<=lon2)]

#time period
begin = 1971
end = 1990
x = np.where((dates <= datetime.datetime(begin,12,31)) & (dates >= datetime.datetime(begin,1,1))); print(x)
x = x[0][0]
y = np.where((dates <= datetime.datetime(end,12,31)) & (dates >= datetime.datetime(end,1,1))); print(y)
y = y[0][11]

data_area = data[x:y+1,:,:][:,(lat<=lat1) & (lat>=lat2),:][:,:,(lon>=lon1) & (lon<=lon2)]

dates_area = dates[x:y+1]

del(begin,end)
del(x,y)
del(lat1,lat2,lon1,lon2)