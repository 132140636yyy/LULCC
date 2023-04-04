# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 09:20:57 2022
@author: Liu Yingying
Define global main LULCC classes and distributions
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
filepath = "C:/Users/26308/Desktop/work/CCI/lcc/"
filename = 'ESACCI-LC-L4-LCCS-Map-300m-P1Y-aggregated-0.500000Deg-2005-v2.0.7.nc'
file = filepath+filename
nf = nc.Dataset(file,'r')
varname = nf.variables.keys()
varname = list(varname)
lat_lcc = np.array(nf.variables[varname[0]][:])
lon_lcc = np.array(nf.variables[varname[1]][:])
#%% read CCI data
os.chdir("C:/Users/26308/Desktop/work/CCI/lcc/")
file_new = np.load('file_new_data.npy')
file_new28 = file_new-file_new[0,:,:,:]
sort_data = np.load('sort_data.npy')
index_data = np.load('index_data.npy')
accum_change = np.load('accum_change.npy')
# missing value
where_sort_zero = np.where(sort_data[:,:,0]==0)
index_data[where_sort_zero[0],where_sort_zero[1],:] = -999
#calculate the map amplification factor
def data_num_m():
  data_num = np.array(np.where(index_data[:,:,0]!=-999))
  lat_data_num = lat_lcc[data_num[0,:]]
  m = np.cos(np.deg2rad(lat_data_num))
  data_num = np.nansum(m)
  return data_num
data_num = data_num_m()
#define the dramatic changed grids
def check_magnitude():
  check_sum = np.where(sort_data[:,:,0]<5*10**(-2))
  index_data[check_sum[0],check_sum[1],:] = -999
  return check_sum,index_data
check_sum,index_data = check_magnitude()
index_data = index_data.astype(int)
#%% define the main LULCC class
def combine_PFT():
  lat_x_data = list()
  lon_y_data = list()
  max_one_type = index_data[:,:,0]
  max_two_type = index_data[:,:,1]
  unique_max_one,unique_where_max,where_max = np.unique(max_one_type,return_index = True,return_inverse = True)
  #find the two main PFTs in each grid
  unique_group = list()
  for i in range(1,unique_where_max.size):
    grid_unique = list(np.where(max_one_type == unique_max_one[i]))
    lat_x = grid_unique[0]
    lon_y = grid_unique[1]
    unique_max_two,unique_where_two,where_two = np.unique(max_two_type[lat_x,lon_y],return_index = True,return_inverse = True)
    for j in range(unique_where_two.size):
      grid_unique_two = list(np.where(max_two_type[lat_x,lon_y] == unique_max_two[j]))
      lat_x_two = lat_x[np.array(grid_unique_two)]
      lon_y_two = lon_y[np.array(grid_unique_two)]
      if unique_max_one[i]<unique_max_two[j]:
        unique_max = [unique_max_one[i]*10+unique_max_two[j]]
        unique_group = unique_group + unique_max
        lat_x_data = lat_x_data + [lat_x_two]
        lon_y_data = lon_y_data + [lon_y_two]
      if unique_max_one[i]>unique_max_two[j]:
        where_small = np.where(unique_group == unique_max_two[j]*10+unique_max_one[i])
        if where_small[0].size != 0:
          lat_x_data[where_small[0][0]] = np.concatenate((lat_x_data[where_small[0][0]], lat_x_two),axis = 1)
          lon_y_data[where_small[0][0]] = np.concatenate((lon_y_data[where_small[0][0]], lon_y_two),axis = 1)
        if where_small[0].size == 0:
          lat_x_data = lat_x_data + [lat_x_two]
          lon_y_data = lon_y_data + [lon_y_two]
  return lat_x_data,lon_y_data,unique_group,where_small
lat_x_data,lon_y_data,unique_group,where_small = combine_PFT()

def find_LULCC_type():
  #seperate the LULCC class with different change direction
  lat_x_group = list()
  lon_y_group = list()
  for i in range(len(lat_x_data)):
    group_lat = lat_x_data[i][0]
    group_lon = lon_y_data[i][0]
    PFT_one = index_data[group_lat[0],group_lon[0],0]
    PFT_two = index_data[group_lat[0],group_lon[0],1]
    lcc_dire = np.zeros((len(group_lat),2))
    lat_x_grid_pos = list()
    lon_y_grid_pos = list()
    lat_x_grid_neg = list()
    lon_y_grid_neg = list()
    for j in range(len(group_lat)):
      grid_data_one = file_new28[:,PFT_one,group_lat[j],group_lon[j]]
      grid_data_two = file_new28[:,PFT_two,group_lat[j],group_lon[j]]
      grid_data_one_mean = np.nanmean(grid_data_one)
      grid_data_two_mean = np.nanmean(grid_data_two)
      if grid_data_one_mean>0:
        lcc_dire[j,0] = 1
      if grid_data_one_mean<0:
        lcc_dire[j,0] = -1
      if grid_data_two_mean>0:
        lcc_dire[j,1] = 1
      if grid_data_two_mean<0:
        lcc_dire[j,1] = -1
      if (lcc_dire[j,0]==lcc_dire[j,1]):
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
#%%
def grid_data_search(accum_change,lat_x_data,lon_y_data):
  index_plot = list()
  diff_plot = list()
  lat_x = list()
  lon_y = list()
  lcc_num = 0
  for i in range(0,len(lat_x_group)):
    data_x = lat_x_group[i]
    data_y = lon_y_group[i]
    lat_m = lat_lcc[data_x]
    m = np.cos(np.deg2rad(lat_m)) #Map amplification factor
    group_num = np.nansum(m)
    diff_data = accum_change[data_x,data_y,:,:]
    #average change from 1992 to 2020 in each LULCC class
    diff_sum = np.nansum(diff_data,axis = 0)
    diff_len = diff_data.shape[0]
    diff_mean = diff_sum/diff_len
    sort_grid = sort_data[data_x,data_y,:]
    sort_num = sort_grid.size/5
    index_grid = index_data[data_x,data_y,0:2]
    index_use = np.unique(index_grid)
    lcc_num = lcc_num + group_num
    if sort_num>400:
      index_plot = index_plot + [index_use]
      diff_plot = diff_plot + [diff_mean]
      lat_x = lat_x + [lat_x_group[i]]
      lon_y = lon_y + [lon_y_group[i]]
  return index_plot,diff_plot,lat_x,lon_y,lcc_num
index_plot,diff_plot,lat_x,lon_y,lcc_num = grid_data_search(accum_change,lat_x_group,lon_y_group)
zero = np.zeros((29,8),dtype=int)
zero_index = np.zeros(2,dtype=int)
diff_plot = diff_plot + [zero]
index_plot = index_plot + [zero_index]
#%% plot the accumulated change of main LULCC classes
font_tick = {'size':12,'weight':'heavy'}
font_title = {'size':14,'weight':'heavy'}
def to_percent(mean_pos_data):
  percent_mean_pos_data = list()
  for i in range(0,len(mean_pos_data)):
    percent_mean_pos_data.append(mean_pos_data[i]*100)
  return percent_mean_pos_data
diff_plot_percent = to_percent(diff_plot)
color_line=['#D8383A', '#FC8A61', '#2F7FC1', '#96C37D','#C497B2', '#A9B8C6', '#BB9727', '#6F6F6F']
text_word = ['S-T','T-S','T-G','C-T','T-C','B-T','T-B','S-C','B-G','B-C','C-U']
year = range(1992,2021,1)
def grid_data_plot(num_x,num_y):
  ax[num_x][num_y].set_prop_cycle(color=color_line)
  ax[num_x][num_y].set_xlim(1991.5, 2020.5)
  ax[num_x][num_y].set_xticks(range(1992,2021,2),fontdict=font_label)
  ax[num_x][num_y].get_xaxis().set_visible(False)
  ax[num_x][num_y].xaxis.set_major_formatter('{x:.0f}')
  labels = ax[num_x][num_y].get_xticklabels()
  plt.setp(labels, rotation=45, horizontalalignment='right')
  ax[num_x][num_y].axhline(y=0, color="black", linestyle="--",lw = 2)
  ax[num_x][num_y].tick_params(axis='both', which='both',labelsize=12,
                 bottom=False, top=False, labelbottom=True,
                 left=False, right=False, labelleft=True)
  labels = ax[num_x][num_y].get_xticklabels() + ax[num_x][num_y].get_yticklabels()
  [label.set_fontname('Times New Roman') for label in labels]
  [label.set_fontweight('medium') for label in labels]

  formatter = mpl.ticker.FormatStrFormatter('%d%%')
  Axis.set_major_formatter(ax[num_x][num_y].yaxis, formatter)
  order = index_plot[num_x*2+num_y]
  line_other = ax[num_x][num_y].plot(year,diff_plot_percent[num_x*2+num_y],lw = 1,alpha=0.8)
  color_main = [color_line[order[0]]]+[color_line[order[1]]]
  ax[num_x][num_y].set_prop_cycle(color=color_main)
  line_main = ax[num_x][num_y].plot(year,diff_plot_percent[num_x*2+num_y][:,order],lw = 2.2)
fig, ax = plt.subplots(int(len(diff_plot)/2),2,dpi=400,figsize=(10,9),sharex = True)
font_title = {'family':'Times New Roman','weight':'heavy','size':16}
font_label = {'family':'Times New Roman','weight':'heavy','size':12}
for num_x in range(0,int(len(diff_plot)/2)):
  for num_y in range(0,2):
    grid_data_plot(num_x,num_y)
    num_order = int(num_x*2+num_y+1)
    string = "C"+str(num_order)
    if num_order<10:
      ascii_num = ord(str(num_order))+48
      ax[num_x][num_y].set_title(chr(ascii_num),loc='left',fontdict=font_title)
    if num_order==10:
      ax[num_x][num_y].set_title('j',loc='left',fontdict=font_title)
    if num_order==11:
      ax[num_x][num_y].set_title('k',loc='left',fontdict=font_title)
    if ax[num_x][num_y].get_subplotspec().is_last_row():
      ax[num_x][num_y].get_xaxis().set_visible(True)
      ax[num_x][num_y].set_xticks(range(1992,2021,2),fontdict=font_label)
    if num_x == 4 and num_y == 1:
      ax[num_x][num_y].get_xaxis().set_visible(True)
      ax[num_x][num_y].set_xticks(range(1992,2021,2),fontdict=font_label)
      labels = ax[num_x][num_y].get_xticklabels()
      plt.setp(labels, rotation=45, horizontalalignment='right')
    if num_x == 5 and num_y == 1:
      ax[num_x][num_y].set_frame_on(False)
      ax[num_x][num_y].xaxis.set_visible(False)
      ax[num_x][num_y].yaxis.set_visible(False)
      ax[num_x][num_y].plot(year,diff_plot[num_x*2+num_y],color = "white",alpha=1,lw=5)
      ax[num_x][num_y].axhline(False)
    if (num_x==0 and num_y==1)|(num_x==3 and num_y==1):
      ax[num_x][num_y].set_ylim(-11,11)
    if (num_x==1 and num_y==0):
      ax[num_x][num_y].set_ylim(-9,9)
    if (num_x==3 and num_y==1):
      ax[num_x][num_y].set_ylim(-13,13)
ax[0][0].set_title('S-T',loc='center',fontdict=font_title)
ax[0][1].set_title('T-S',loc='center',fontdict=font_title)
ax[1][0].set_title('T-G',loc='center',fontdict=font_title)
ax[1][1].set_title('C-T',loc='center',fontdict=font_title)
ax[2][0].set_title('T-C',loc='center',fontdict=font_title)
ax[2][1].set_title('B-T',loc='center',fontdict=font_title)
ax[3][0].set_title('T-B',loc='center',fontdict=font_title)
ax[3][1].set_title('S-C',loc='center',fontdict=font_title)
ax[4][0].set_title('B-G',loc='center',fontdict=font_title)
ax[4][1].set_title('B-C',loc='center',fontdict=font_title)
ax[5][0].set_title('C-U',loc='center',fontdict=font_title)
fig.text(-0.02, 0.5,'Change of PFT Fractions (%)', va='center', rotation='vertical',fontdict = font_title)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                wspace=0.2, hspace=0.3)
fig.tight_layout()
plt.show()
#%% label
fig, ax = plt.subplots(dpi=600,figsize=(5, 4))
patch1 = mpatches.Patch(color='#D8383A', label='Tree')
patch2 = mpatches.Patch(color='#FC8A61', label='Shrub')
patch3 = mpatches.Patch(color='#2F7FC1', label='Nature Grass')
patch4 = mpatches.Patch(color='#96C37D', label='Managed Grass')
patch5 = mpatches.Patch(color='#C497B2', label='Bare Soil')
patch6 = mpatches.Patch(color='#A9B8C6', label='Water')
patch7 = mpatches.Patch(color='#BB9727', label='Snow/Ice')
patch8 = mpatches.Patch(color='#6F6F6F', label='Urban areas')
legend = ax.legend(handles=[patch1,patch2,patch3,patch4,patch5,patch6,patch7,patch8],ncol = 4)
ax.set_frame_on(False)
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
frame = legend.get_frame()
frame.set_alpha(0)
frame.set_facecolor('none')
plt.rc('font',family='Times New Roman',weight='semibold',size=14)
plt.show()
data_Z = np.zeros((360,720))
for i in range(0,len(lat_x)):
  data_Z[lat_x[i],lon_y[i]] = i+1.5
data_Z = np.where(data_Z == 0 ,np.nan,data_Z)
#%% spatial distribution and fractions of main LULCC classes
def group12_plot(x,tick_label):
  # cmap_use = plt.get_cmap("jet")
  cmap_use = cm.get_cmap('Paired_r',11)
  fig = plt.figure(dpi = 600,figsize=(10,7))
  # ax1 = fig.add_subplot(1,2,1,projection=ccrs.PlateCarree())
  ax1 = fig.add_axes([0.1,0.1,0.8,0.8],projection=ccrs.PlateCarree())
  gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0, linestyle='-.')
  ax1.coastlines(alpha=1.,linestyle='-',color = "black",lw=0.7)
  gl.xlabels_top = False
  gl.ylabels_right = False
  gl.xlocator = mticker.FixedLocator([-180, -120, -60, 0, 60, 120, 180])
  gl.ylocator = mticker.FixedLocator([-90,-60, -30, 0, 30, 60,90])
  gl.ypadding = 20
  gl.xpadding = 35
  gl.xlabel_style = {'size': 12, 'color': 'black'}
  gl.ylabel_style = {'size': 12, 'color': 'black'}
  lon_formatter = LongitudeFormatter(zero_direction_label=True)
  lat_formatter = LatitudeFormatter()
  ax1.xaxis.set_major_formatter(lon_formatter)
  ax1.yaxis.set_major_formatter(lat_formatter)
  lcc_per = lcc_num/data_num*100
  tick = np.linspace(1.5,11.5,11)
  ax1.set_extent([-179.75, 179.75, -89.75, 89.75], ccrs.PlateCarree())
  where_x1 = np.array(np.where(~np.isnan(data_Z)))[1]
  where_y1 = np.array(np.where(~np.isnan(data_Z)))[0]
  x1 = lon_lcc[where_x1]
  y1 = lat_lcc[where_y1]
  con_plot = ax1.scatter(x1,y1,c=data_Z[where_y1,where_x1],transform=ccrs.PlateCarree(),cmap = cmap_use,s=0.3)
  cb = plt.colorbar(con_plot,fraction = 0.05,pad = 0.08,label = 'LULCC Class',orientation="horizontal",boundaries=np.linspace(1,12,12), cmap=cmap_use,ticks=tick)
  cb.ax.set_xticklabels(tick_label,fontsize=12)

  tab20 = cm.get_cmap('Paired_r',11)
  colors = tab20(np.linspace(0, 1, 11))
  width = 0.3
  ax2 = fig.add_axes([0.58,0.26,0.2,0.1])
  ax2.spines['right'].set_visible(False)
  ax2.spines['top'].set_visible(False)
  ax2.bar(x,rate, width,color = colors,edgecolor="gray")
  ax2.set_xticks(x, labels,fontsize=9.4,rotation=305)
  ax2.set_ylim(0,max(rate)+2)
  # ax2.grid(linestyle='-.',alpha=.5)
  ax2.set_title('Fractions (%)',fontdict=font_label)
  formatter = mpl.ticker.FormatStrFormatter('%d%%')
  Axis.set_major_formatter(ax2.yaxis, formatter)
  plt.rcParams['font.size'] = 12
  plt.rc('font',family='Times New Roman')

  ax3 = fig.add_axes([0.15,0.22,0.15,0.15])
  ax3.pie(x=[lcc_per,100-lcc_per],colors=['#FF4500','#DCDCDC'])
  ax1.arrow(-85,-58,75,0, length_includes_head=False,head_width=2.5, fc='grey', ec='k')
  text_str = ' %.2f%%' %(lcc_per)
  text_lcc = ' %.2f%%' %(rate[4])
  font_small = {'family':'Times New Roman','weight':'heavy','size':10}
  ax1.text(-175,-70,'Unchanged',fontdict=font_small)
  ax1.text(-108,-40,'Changed',fontdict=font_small)
  ax1.text(-108,-33,text_str,fontdict=font_small)
  ax2.text(2.2,24,text_lcc,fontdict=font_small)
  return colors,con_plot

tick_label = text_word
labels = text_word
def nums_statis():
  nums = np.zeros(len(lat_x))
  for i in range(0,len(lat_x)):
    lat_each = lat_x[i]
    nums_each = len(lat_each)
    lat_each = lat_lcc[lat_x[i]]
    m = np.cos(np.deg2rad(lat_each))
    nums[i] = np.nansum(m)
  num_all = sum(nums)
  print(num_all)
  rate = nums/num_all*100
  return rate
rate = nums_statis()
x = np.arange(len(labels))/2  # the label locations
colors,con_plot = group12_plot(x,text_word)
plt.show()
# plt.savefig('E:/paper/LULCC_impact/plot/Figure2.pdf',dpi=600)
# plt.close()