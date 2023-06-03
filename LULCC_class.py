# -*- coding: utf-8 -*-
"""
Define global main LULCC classes and distributions
author: Liu Yingying
"""

import numpy as np
import netCDF4 as nc
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap,BoundaryNorm
import matplotlib.patches as mpatches
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.ticker as mticker
from matplotlib.axis import Axis

#%%read lat and lon data in CCI dataset
os.chdir("E:/paper/LULCC_impact/data")
nf       = nc.Dataset('LULC_data.nc','r')
varname  = nf.variables.keys()
varname  = list(varname)
lat_lcc  = np.array(nf.variables[varname[0]][:])
lon_lcc  = np.array(nf.variables[varname[1]][:])

# read proportion of PFT from CCI dataset
PFT     = np.array(nf.variables[varname[4]][:])
PFT28   = PFT-PFT[0,:,:,:]

#%%ranked in decreasing order of proportion of PFTs in each grid
def each_grid_rate(file_data):
  accum_change = np.zeros((360,720,29,8))
  sort_data = np.zeros((360,720,5))
  index_data = np.zeros((360,720,5))
  for lat in range(0,360):
    for lon in range(0,720):
      empty_array = list()
      for year in range(0,29):
        each_year_type = file_data[year,:,lat,lon]
        #calculate the proportion of each PFT in each grids
        sum_type = np.nansum(each_year_type)
        each_grid = each_year_type/sum_type
        empty_array.append(each_grid[np.newaxis,:])
        diff_year_type = np.concatenate(empty_array,axis = 0)
      year_change = diff_year_type-diff_year_type[0,:]
      accum_change[lat,lon,:,:]=year_change
      #rank PFT
      sort_type_rate = -(np.sort(-np.abs(year_change[28,:])))
      index_type_rate = np.argsort(-np.abs(year_change[28,:]))
      sort_data[lat,lon,:] = sort_type_rate[0:5]
      index_data[lat,lon,:] = index_type_rate[0:5]
  return sort_data,index_data,accum_change

# sort_data    : high to low proportion of PFT variations in each grid
# index_data   : PFT in order of rank from high to low according to PFT vairaitons
# accum_change : accumulation of changes in each PFT from 1992-2020
sort_data,index_data,accum_change = each_grid_rate(PFT)

# missing value
where_sort_zero = np.where(sort_data[:,:,0]==0)
index_data[where_sort_zero[0],where_sort_zero[1],:] = -999

#calculate the map amplification factor
def data_num_m():
  data_num     = np.array(np.where(index_data[:,:,0]!=-999))
  lat_data_num = lat_lcc[data_num[0,:]]
  m            = np.cos(np.deg2rad(lat_data_num))
  data_num     = np.nansum(m)
  return data_num
data_num       = data_num_m()

#define the dramatic changed grids
def check_magnitude():
  check_sum = np.where(sort_data[:,:,0]<5*10**(-2))
  index_data[check_sum[0],check_sum[1],:] = -999
  return check_sum,index_data

check_sum,index_data = check_magnitude()
index_data           = index_data.astype(int)

#%% find the groups with same PFT conversion
def combine_PFT():
  lat_x_data   = list()
  lon_y_data   = list()
  # the first main change PFT of grid
  max_one_type = index_data[:,:,0]
  # the second main change PFT of grid
  max_two_type = index_data[:,:,1]
  #find the all types of main change PFTs across the global
  unique_max_one,unique_where_max,where_max = np.unique(max_one_type,return_index = True,return_inverse = True)

  #find the two main change PFTs in each grid
  unique_group = list()
  for i in range(1,unique_where_max.size):
    # classify the grids with the same first main PFT
    grid_unique = list(np.where(max_one_type == unique_max_one[i]))
    lat_x       = grid_unique[0]
    lon_y       = grid_unique[1]
    unique_max_two,unique_where_two,where_two = np.unique(max_two_type[lat_x,lon_y],return_index = True,return_inverse = True)

    for j in range(unique_where_two.size):
      # classify the second main PFT
      grid_unique_two = list(np.where(max_two_type[lat_x,lon_y] == unique_max_two[j]))
      lat_x_two = lat_x[np.array(grid_unique_two)]
      lon_y_two = lon_y[np.array(grid_unique_two)]
      # combine the groups with same first main and second main PFT
      # such as combine the Tree>Shrub and Shrub>Tree group
      if unique_max_one[i]<unique_max_two[j]:
        unique_max   = [unique_max_one[i]*10 + unique_max_two[j]]
        unique_group = unique_group          + unique_max
        lat_x_data   = lat_x_data            + [lat_x_two]
        lon_y_data   = lon_y_data            + [lon_y_two]
      if unique_max_one[i]>unique_max_two[j]:
        where_small  = np.where(unique_group == unique_max_two[j]*10+unique_max_one[i])
        if where_small[0].size != 0:
          lat_x_data[where_small[0][0]]  = np.concatenate((lat_x_data[where_small[0][0]], lat_x_two),axis = 1)
          lon_y_data[where_small[0][0]]  = np.concatenate((lon_y_data[where_small[0][0]], lon_y_two),axis = 1)
        if where_small[0].size == 0:
          lat_x_data = lat_x_data + [lat_x_two]
          lon_y_data = lon_y_data + [lon_y_two]
  return lat_x_data,lon_y_data,unique_group,where_small
lat_x_data,lon_y_data,unique_group,where_small = combine_PFT()

# seperate the PFT with different change direction from the groups with same PFT conversion
def find_LULCC_type():
  #seperate the LULCC class with different change direction
  lat_x_group = list()
  lon_y_group = list()
  for i in range(len(lat_x_data)):
    group_lat = lat_x_data[i][0]
    group_lon = lon_y_data[i][0]
    PFT_one   = index_data[group_lat[0],group_lon[0],0]
    PFT_two   = index_data[group_lat[0],group_lon[0],1]
    lcc_dire  = np.zeros((len(group_lat),2))
    lat_x_grid_pos = list()
    lon_y_grid_pos = list()
    lat_x_grid_neg = list()
    lon_y_grid_neg = list()

    for j in range(len(group_lat)):
      grid_data_one      = PFT28[:,PFT_one,group_lat[j],group_lon[j]]
      grid_data_two      = PFT28[:,PFT_two,group_lat[j],group_lon[j]]
      grid_data_one_mean = np.nanmean(grid_data_one)
      grid_data_two_mean = np.nanmean(grid_data_two)

      # >0 increase ; <0 decrease
      if grid_data_one_mean>0:
        lcc_dire[j,0] = 1
      if grid_data_one_mean<0:
        lcc_dire[j,0] = -1
      if grid_data_two_mean>0:
        lcc_dire[j,1] = 1
      if grid_data_two_mean<0:
        lcc_dire[j,1] = -1
      # all increase or decrease
      if (lcc_dire[j,0] ==lcc_dire[j,1]):
        if lcc_dire[j,0]==-1:
          index_data[group_lat[j],group_lon[j],:]=-999
          continue
      if lcc_dire[j,0]==1:
        lat_x_grid_pos = lat_x_grid_pos + [group_lat[j]]
        lon_y_grid_pos = lon_y_grid_pos + [group_lon[j]]
      if lcc_dire[j,0]==-1:
        lat_x_grid_neg = lat_x_grid_neg + [group_lat[j]]
        lon_y_grid_neg = lon_y_grid_neg + [group_lon[j]]
    lat_x_group = lat_x_group + [lat_x_grid_pos] + [lat_x_grid_neg]
    lon_y_group = lon_y_group + [lon_y_grid_pos] + [lon_y_grid_neg]
  return lat_x_group,lon_y_group
lat_x_group,lon_y_group = find_LULCC_type()

#%% define the LULCC classes
def grid_data_search(accum_change,lat_x_data,lon_y_data):
  index_plot  = list()
  diff_plot   = list()
  lat_x       = list()
  lon_y       = list()
  lcc_num     = 0

  # find the grids in each LULCC classes
  for i in range(0,len(lat_x_group)):
    data_x    = lat_x_group[i]
    data_y    = lon_y_group[i]
    lat_m     = lat_lcc[data_x]
    m         = np.cos(np.deg2rad(lat_m)) #Map amplification factor
    group_num = np.nansum(m)
    diff_data = accum_change[data_x,data_y,:,:]

    #average change from 1992 to 2020 in each LULCC class
    diff_sum  = np.nansum(diff_data,axis = 0)
    diff_len  = diff_data.shape[0]
    diff_mean = diff_sum/diff_len
    sort_grid = sort_data[data_x,data_y,:]
    sort_num  = sort_grid.size/5
    index_grid= index_data[data_x,data_y,0:2]
    index_use = np.unique(index_grid)
    lcc_num   = lcc_num + group_num

    # reserve the LULCC classes with more than 400 grids
    if sort_num>400:
      index_plot = index_plot+ [index_use]
      diff_plot  = diff_plot + [diff_mean]
      lat_x      = lat_x     + [lat_x_group[i]]
      lon_y      = lon_y     + [lon_y_group[i]]

  data_Z = np.zeros((360,720))
  data_Z = np.where(data_Z == 0 ,np.nan,data_Z)
  for i in range(0,len(lat_x)):
    data_Z[lat_x[i],lon_y[i]] = i

  return index_plot,diff_plot,lat_x,lon_y,lcc_num,data_Z

# index_plot : Two main PFTs in global main LULCC classes
# diff_plot  : Average changes in each main LULCC classes from 1992-2020
# lat_x      : lat of grids in each LULCC classes
# lon_y      : lon of grids in each LULCC classes
# lcc_num    : total number of dramatically changed grids after considering map amplification factor
# data_Z     : array of LULCC classes
index_plot,diff_plot,lat_x,lon_y,lcc_num,data_Z = grid_data_search(accum_change,lat_x_group,lon_y_group)
