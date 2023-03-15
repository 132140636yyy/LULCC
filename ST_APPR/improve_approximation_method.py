# -*- coding: utf-8 -*-
"""
Improved space-for-time approximation method
"""
#%%average of each seven years  1992-1999 1999-2006 2006-2013 2013-2020
year_str = np.zeros((360,720,11))
year_str = np.where(year_str==0,np.nan,year_str)
lcc_chag = np.zeros((360,720,11))
prep_rate = np.zeros((360,720,11))
prep_all_plot = np.zeros((360,720))
evap_rate = np.zeros((360,720,11))
evap_all_plot = np.zeros((360,720))
prep_lcc_rate = np.zeros((360,720,11))
prep_lcc_plot = np.zeros((360,720))
evap_lcc_rate = np.zeros((360,720,11))
evap_lcc_plot = np.zeros((360,720))
prep_chag = np.zeros((360,720,11))
evap_chag = np.zeros((360,720,11))
year_type_group = np.array([[0,1,0],
                     [1,2,0],
                     [2,3,0],
                     [3,2,1],
                     [4,3,1],
                     [5,3,2]],dtype=int)
year_devide = np.array([[1,8],
               [8,15],
               [15,22],
               [22,29]],dtype=int)
grid_main = np.array([0,1,2,0,3,0,4,3,2,3,7],dtype=int)

def LULCC_imp(i,j,prep_all28,prep_obs28,grid_type):
  prep_lcc_avgdata = np.zeros(4)
  prep_unchag_avgdata = np.zeros(4)
  prep_lcc_data = np.zeros(29)
  prep_obs_grid = prep_obs28[:,i,j]
  def other_change(i,j,year_choose,grid_type,prep_all28,prep_obs28,grid_size):
    #dramatic changed grids
    prep_area_obs = prep_obs28[year_devide[year_choose,0]:year_devide[year_choose,1],i-grid_size:i+grid_size+1,j-grid_size:j+grid_size+1]
    prep_area_all = prep_all28[year_devide[year_choose,0]:year_devide[year_choose,1],i-grid_size:i+grid_size+1,j-grid_size:j+grid_size+1]
    #unchanged grids
    prep_obs = prep_obs28[:,i-grid_size:i+grid_size+1,j-grid_size:j+grid_size+1]
    prep_area = prep_all28[:,i-grid_size:i+grid_size+1,j-grid_size:j+grid_size+1]
    prep_unchag_all = np.where(~np.isnan(prep_obs),np.nan,prep_area)
    prep_area_unchag = np.where(~np.isnan(prep_area_obs),np.nan,prep_area_all)
    prep_unchag = np.nanmean(prep_area_unchag)
    prep_all_range = prep_all28[year_devide[year_choose,0]:year_devide[year_choose,1],i,j]
    prep_lcc = np.nanmean(prep_all_range)
    prep_lcc_data[year_devide[year_choose,0]:year_devide[year_choose,1]] = prep_all_range-prep_unchag
    return prep_lcc,prep_unchag,prep_lcc_data
  for year_choose in range(0,4):
    prep_lcc_avgdata[year_choose],prep_unchag_avgdata[year_choose],prep_lcc_data = other_change(i,j,year_choose,grid_type,prep_all28,prep_obs28,3)
  for year_type in range(0,6):
    if year_str[i,j,int(grid_type)]==year_type_group[year_type,0]:
      #calculate the changes caused by LULCC between corresponding period
      time_range1 = year_type_group[year_type,1]
      time_range2 = year_type_group[year_type,2]
      prep_lcc_cause = (prep_lcc_avgdata[time_range1]-prep_lcc_avgdata[time_range2])-(prep_unchag_avgdata[time_range1]-prep_unchag_avgdata[time_range2])
      #calculate the total changes between corresponding period
      prep_total_cause = prep_lcc_avgdata[time_range1]-prep_lcc_avgdata[time_range2]
      #two series for ttest
      prep_lcc1 = prep_lcc_data[year_devide[time_range1,0]:year_devide[time_range1,1]]
      prep_lcc2 = prep_lcc_data[year_devide[time_range2,0]:year_devide[time_range2,1]]
      prep_all1 = prep_obs_grid[year_devide[time_range1,0]:year_devide[time_range1,1]]
      prep_all2 = prep_obs_grid[year_devide[time_range2,0]:year_devide[time_range2,1]]
  return prep_lcc_cause,prep_total_cause,prep_lcc1,prep_lcc2,prep_all1,prep_all2

def find_year(i,j,grid_type,lcc_data):
  for year_type in range(0,6):
    if year_str[i,j,int(grid_type)]==year_type_group[year_type,0]:
      lcc_mean = lcc_data[year_type_group[year_type,1]]-lcc_data[year_type_group[year_type,2]]
  return lcc_mean

def per_percent_imp(grid_size,prep_all28,prep_obs28):
  #calculate the LULCC impact in each dramatic changed grids
  prep_lcc_p = np.zeros((360,720))
  evap_lcc_p = np.zeros((360,720))
  wa_lcc_p = np.zeros((360,720))
  wa_all_p = np.zeros((360,720))
  prep_lcc_cause = np.zeros((360,720))
  evap_lcc_cause = np.zeros((360,720))
  prep_total_cause = np.zeros((360,720))
  evap_total_cause = np.zeros((360,720))
  WA_MSWEP_lcc_rate = np.zeros((360,720,11))
  WA_MSWEP_all_rate = np.zeros((360,720,11))
  wa_lcc_pval = np.zeros((360,720))
  wa_all_pval = np.zeros((360,720))
  n = 28
  t_test = stats.t.ppf(0.95,n-1)
  for i in range(0+grid_size,360-grid_size):
    for j in range(0+grid_size,720-grid_size):
      grid_type = data_Z[i,j]-1.5
      #PFT change in CCI data
      if abs(grid_type)<20:
        land_main_type = grid_main[int(grid_type)]
        prep_obs_grid = prep_obs28[:,i,j]
        evap_obs_grid = evap_obs28[:,i,j]
        year_mean = np.zeros(6)
        lcc_data = np.zeros(4)
        lcc_data[0] = np.nanmean(file_new28[1:8,int(land_main_type),i,j])
        lcc_data[1] = np.nanmean(file_new28[8:15,int(land_main_type),i,j])
        lcc_data[2] = np.nanmean(file_new28[15:22,int(land_main_type),i,j])
        lcc_data[3] = np.nanmean(file_new28[22:29,int(land_main_type),i,j])
        def find_lcc_year(lcc_year_type):
          year_mean[lcc_year_type] = lcc_data[year_type_group[lcc_year_type,1]]-lcc_data[year_type_group[lcc_year_type,2]]
          return year_mean
        for lcc_year_type in range(0,6):
          year_mean = find_lcc_year(lcc_year_type)
        year_mean_abs = list(np.abs(year_mean))
        #select the periods with max and min PFT
        year_str[i,j,int(grid_type)] = year_mean_abs.index(max(year_mean_abs))
        lcc_mean = find_year(i,j,grid_type,lcc_data)
        if abs(lcc_mean)<0.008:
          lcc_chag[i,j,:] = np.nan
          continue
        else:
          lcc_chag[i,j,int(grid_type)] = lcc_mean
        #changes of precipitation and evapotranspiration in corresponding periods
        prep_lcc_cause[i,j],prep_total_cause[i,j],prep_lcc_data1,prep_lcc_data2,prep_all_data1,prep_all_data2 = LULCC_imp(i,j,prep_all28,prep_obs28,grid_type)
        evap_lcc_cause[i,j],evap_total_cause[i,j],evap_lcc_data1,evap_lcc_data2,evap_all_data1,evap_all_data2 = LULCC_imp(i,j,evap_all28,evap_obs28,grid_type)
        #the total impact of per percent LULCC on precipitation and evapotranspiration
        prep_rate[i,j,int(grid_type)] = prep_total_cause[i,j]/lcc_chag[i,j,int(grid_type)]*0.01
        evap_rate[i,j,int(grid_type)] = evap_total_cause[i,j]/lcc_chag[i,j,int(grid_type)]*0.01
        prep_all_plot[i,j] = prep_total_cause[i,j]/lcc_chag[i,j,int(grid_type)]*0.01
        evap_all_plot[i,j] = evap_total_cause[i,j]/lcc_chag[i,j,int(grid_type)]*0.01
        WA_MSWEP_all_rate[i,j,int(grid_type)] = (prep_total_cause[i,j]-evap_total_cause[i,j])/lcc_chag[i,j,int(grid_type)]*0.01
        #the selected impact of per percent LULCC on precipitation and evapotranspiration
        prep_lcc_rate[i,j,int(grid_type)] = prep_lcc_cause[i,j]/lcc_chag[i,j,int(grid_type)]*0.01
        evap_lcc_rate[i,j,int(grid_type)] = evap_lcc_cause[i,j]/lcc_chag[i,j,int(grid_type)]*0.01
        prep_lcc_plot[i,j] = prep_lcc_cause[i,j]/lcc_chag[i,j,int(grid_type)]*0.01
        evap_lcc_plot[i,j] = evap_lcc_cause[i,j]/lcc_chag[i,j,int(grid_type)]*0.01
        WA_MSWEP_lcc_rate[i,j,int(grid_type)] = (prep_lcc_cause[i,j]-evap_lcc_cause[i,j])/lcc_chag[i,j,int(grid_type)]*0.01
        def grid_ttest(prep_all_data1,prep_all_data2):
          tstat, pval = stats.ttest_ind(a=prep_all_data1, b=prep_all_data2, alternative="two-sided")
          if pval>0.05:
            p_grid_test = 1
          else:
            p_grid_test = np.nan
          return p_grid_test,pval
        prep_lcc_p[i,j],pval = grid_ttest(prep_lcc_data1,prep_lcc_data2)
        evap_lcc_p[i,j],pval = grid_ttest(evap_lcc_data1,evap_lcc_data2)
        wa_lcc_p[i,j],wa_lcc_pval[i,j] = grid_ttest(prep_lcc_data1-evap_lcc_data1,prep_lcc_data2-evap_lcc_data2)
        wa_all_p[i,j],wa_all_pval[i,j] = grid_ttest(prep_all_data1-evap_all_data1,prep_all_data2-evap_all_data2)
  prep_lcc_p = np.where(prep_lcc_p==0,np.nan,prep_lcc_p)
  evap_lcc_p = np.where(evap_lcc_p==0,np.nan,evap_lcc_p)
  return year_str,lcc_chag,prep_lcc_cause,evap_lcc_cause,prep_rate,evap_rate,prep_all_plot,evap_all_plot,prep_lcc_rate,evap_lcc_rate,prep_lcc_plot,evap_lcc_plot,prep_lcc_p,evap_lcc_p,wa_lcc_p,wa_all_p,WA_MSWEP_all_rate,WA_MSWEP_lcc_rate
year_str,lcc_chag,prep_lcc_cause,evap_lcc_cause,prep_rate,evap_rate,prep_all_plot,evap_all_plot,prep_lcc_rate,evap_lcc_rate,prep_lcc_plot,evap_lcc_plot,prep_lcc_p,evap_lcc_p,wa_lcc_p,wa_all_p,WA_MSWEP_all_rate,WA_MSWEP_lcc_rate = per_percent_imp(3,prep_all28,prep_obs28)
lcc_chag = np.where(lcc_chag==0,np.nan,lcc_chag)
evap_rate = np.where(evap_rate==0,np.nan,evap_rate)
evap_all_plot = np.where(evap_all_plot==0,np.nan,evap_all_plot)
evap_lcc_rate = np.where(evap_lcc_rate==0,np.nan,evap_lcc_rate)
evap_lcc_plot = np.where(evap_lcc_plot==0,np.nan,evap_lcc_plot)
def np_where(prep_rate,prep_all_plot,prep_lcc_rate,prep_lcc_plot,WA_MSWEP_all_rate,WA_MSWEP_lcc_rate):
  prep_rate = np.where(prep_rate==0,np.nan,prep_rate)
  prep_all_plot = np.where(prep_all_plot==0,np.nan,prep_all_plot)
  prep_lcc_rate = np.where(prep_lcc_rate==0,np.nan,prep_lcc_rate)
  prep_lcc_plot = np.where(prep_lcc_plot==0,np.nan,prep_lcc_plot)
  WA_MSWEP_all_rate = np.where(WA_MSWEP_all_rate==0,np.nan,WA_MSWEP_all_rate)
  WA_MSWEP_lcc_rate = np.where(WA_MSWEP_lcc_rate==0,np.nan,WA_MSWEP_lcc_rate)
  return prep_rate,prep_all_plot,prep_lcc_rate,prep_lcc_plot,WA_MSWEP_all_rate,WA_MSWEP_lcc_rate
prep_rate,prep_all_plot,prep_lcc_rate,prep_lcc_plot,WA_MSWEP_all_rate,WA_MSWEP_lcc_rate = np_where(prep_rate,prep_all_plot,prep_lcc_rate,prep_lcc_plot,WA_MSWEP_all_rate,WA_MSWEP_lcc_rate)

#%% t-test
def m_value(prep_lcc_rate):
  m = np.cos(np.deg2rad(lat_lcc))
  prep_new = np.zeros((360,720))
  for j in range(0,360):
    if np.all(np.isnan(prep_lcc_rate[j,:]))==True:
      prep_new[j,:] = prep_lcc_rate[j,:]
    else:
      prep_new[j,:] = prep_lcc_rate[j,:]*m[j]
  return prep_new

def t_test(prep_lcc_rate):
  prep_plot_rate = np.zeros(11)
  t_test1 = np.zeros(11)
  t_test2 = np.zeros(11)
  t_test3 = np.zeros(11)
  t = np.zeros(11)
  p = np.zeros(11)
  prep_new = np.zeros((360,720,11))
  for i in range(0,11):
    prep_new[:,:,i] = m_value(prep_lcc_rate[:,:,i])
    prep_mean = np.nanmean(prep_new[:,:,i])
    prep_plot_rate[i] = np.nanmean(prep_new[:,:,i])
    prep_std = np.nanstd(prep_new[:,:,i])
    prep_where = np.unique(prep_new[:,:,i])
    n = int(np.array(prep_where.shape))
    t[i] = abs(prep_mean)/(prep_std/np.sqrt(n))
    t_test1[i] = stats.t.ppf(0.95,n-1)
    t_test2[i] = stats.t.ppf(0.99,n-1)
    t_test3[i] = stats.t.ppf(0.999,n-1)
    if t[i]>t_test1[i]:
      p[i] = 1
    if t[i]>t_test2[i]:
      p[i] = 2
    if t[i]>t_test3[i]:
      p[i] = 3
  p = np.where(p==0,np.nan,p)
  return prep_plot_rate,p
prep_plot_lcc,p_prep_lcc = t_test(prep_lcc_rate)
evap_plot_lcc,p_evap_lcc = t_test(evap_lcc_rate)
prep_plot_all,p_prep_all = t_test(prep_rate)
evap_plot_all,p_evap_all = t_test(evap_rate)
WA_plot_all,p_WA_all = t_test(WA_MSWEP_all_rate)
WA_plot_lcc,p_WA_lcc = t_test(WA_MSWEP_lcc_rate)