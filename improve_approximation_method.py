# -*- coding: utf-8 -*-
"""
Improved space-for-time approximation method
author: Liu Yingying
"""

#%% average of each seven years  1992-1999 1999-2006 2006-2013 2013-2020

year_str      = np.zeros((360,720,11)) # temporary variable for save the biggest change period, the third dimension is used to identify LULCC class
year_str      = np.where(year_str==0,np.nan,year_str)
lcc_chag      = np.zeros((360,720,11))
prep_rate     = np.zeros((360,720,11))
prep_tot_plot = np.zeros((360,720))
evap_rate     = np.zeros((360,720,11))
evap_tot_plot = np.zeros((360,720))
prep_lcc_rate = np.zeros((360,720,11))
prep_lcc_plot = np.zeros((360,720))
evap_lcc_rate = np.zeros((360,720,11))
evap_lcc_plot = np.zeros((360,720))

# temporary variable for calculate the differences between each post-division periods
year_type_group = np.array([[0,1,0],
                            [1,2,0],
                            [2,3,0],
                            [3,2,1],
                            [4,3,1],
                            [5,3,2]],dtype=int)

# split time period (1993-1999,2000-2006,2007-2013,2014-2020)
year_devide = np.array([[1,8],
                        [8,15],
                        [15,22],
                        [22,29]],dtype=int)

# Main PFT of each LULCC classes (0:Tree  1:Shrub  2:Grass  3:Cropland  4:Bare Soil  5:Water  6:Snow/Ice  7:Urban areas)
grid_main = np.array([0,1,2,0,3,0,4,3,2,3,7],dtype=int)

# exclude the large-scale impacts on water availability variations
def LULCC_imp(i,j,prep_tot28,prep_obs28,grid_type):
    # prep_tot28: changes in all grids of land
    # prep_obs28: changes in dramatically changed grids of land
    # grid_type : variable for identify the LULCC class
    prep_lcc_avgdata    = np.zeros(4)
    prep_unchag_avgdata = np.zeros(4)
    prep_lcc_data       = np.zeros(29)
    prep_obs_grid       = prep_obs28[:,i,j]

    def avg_change(i,j,year_choose,prep_tot28,prep_obs28,grid_size):
    # calculate average LULCC impacts
        #dramatic changed grids
        prep_area_obs    = prep_obs28[year_devide[year_choose,0]:year_devide[year_choose,1],i-grid_size:i+grid_size+1,j-grid_size:j+grid_size+1]
        prep_area_tot    = prep_tot28[year_devide[year_choose,0]:year_devide[year_choose,1],i-grid_size:i+grid_size+1,j-grid_size:j+grid_size+1]

        #unchanged grids
        prep_obs         = prep_obs28[:,i-grid_size:i+grid_size+1,j-grid_size:j+grid_size+1]
        prep_area        = prep_tot28[:,i-grid_size:i+grid_size+1,j-grid_size:j+grid_size+1]
        prep_area_unchag = np.where(~np.isnan(prep_area_obs),np.nan,prep_area_tot)
        prep_unchag      = np.nanmean(prep_area_unchag)
        prep_tot_range   = prep_tot28[year_devide[year_choose,0]:year_devide[year_choose,1],i,j]
        prep_lcc         = np.nanmean(prep_tot_range)

        # calculate difference between each two period
        prep_lcc_data[year_devide[year_choose,0]:year_devide[year_choose,1]] = prep_tot_range-prep_unchag

        return prep_lcc,prep_unchag,prep_lcc_data

    for year_choose in range(0,4):
    # calculate the average LULCC impacts in each periods
        prep_lcc_avgdata[year_choose],prep_unchag_avgdata[year_choose],prep_lcc_data = avg_change(i,j,year_choose,prep_tot28,prep_obs28,3)

    for year_type in range(0,6):
    # calculate precipitation and evapotranspiration changes between the two periods with biggest land cover change
        if year_str[i,j,int(grid_type)]==year_type_group[year_type,0]:

            # calculate the changes caused by LULCC between corresponding period
            time_range1    = year_type_group[year_type,1]
            time_range2    = year_type_group[year_type,2]
            prep_lcc_cause = (prep_lcc_avgdata[time_range1]-prep_lcc_avgdata[time_range2])-(prep_unchag_avgdata[time_range1]-prep_unchag_avgdata[time_range2])

            # calculate the total changes between corresponding period
            prep_total_cause = prep_lcc_avgdata[time_range1]-prep_lcc_avgdata[time_range2]

            # two series for ttest
            prep_lcc1 = prep_lcc_data[year_devide[time_range1,0]:year_devide[time_range1,1]]
            prep_lcc2 = prep_lcc_data[year_devide[time_range2,0]:year_devide[time_range2,1]]
            prep_tot1 = prep_obs_grid[year_devide[time_range1,0]:year_devide[time_range1,1]]
            prep_tot2 = prep_obs_grid[year_devide[time_range2,0]:year_devide[time_range2,1]]

    return prep_lcc_cause,prep_total_cause,prep_lcc1,prep_lcc2,prep_tot1,prep_tot2

def find_year(i,j,grid_type,lcc_data):
    for year_type in range(0,6):
    # find the corresponding periods with biggest land cover change
        if year_str[i,j,int(grid_type)]==year_type_group[year_type,0]:
            lcc_mean = lcc_data[year_type_group[year_type,1]]-lcc_data[year_type_group[year_type,2]]
    return lcc_mean

def per_percent_imp(grid_size,prep_tot28,prep_obs28):
    #calculate the LULCC impact in each dramatic changed grids
    prep_lcc_p        = np.zeros((360,720))
    evap_lcc_p        = np.zeros((360,720))
    wa_lcc_p          = np.zeros((360,720))
    wa_tot_p          = np.zeros((360,720))
    prep_lcc_cause    = np.zeros((360,720))
    evap_lcc_cause    = np.zeros((360,720))
    prep_total_cause  = np.zeros((360,720))
    evap_total_cause  = np.zeros((360,720))
    WA_MSWEP_lcc_rate = np.zeros((360,720,11))
    WA_MSWEP_tot_rate = np.zeros((360,720,11))
    wa_lcc_pval       = np.zeros((360,720))
    wa_tot_pval       = np.zeros((360,720))
    n                 = 28
    t_test            = stats.t.ppf(0.95,n-1)

    for i in range(0+grid_size,360-grid_size):
        for j in range(0+grid_size,720-grid_size):
            grid_type = data_Z[i,j]-1.5 # data_Z : array of LULCC classes

            #PFT change in CCI data
            if abs(grid_type)<20:
                land_main_type = grid_main[int(grid_type)]
                prep_obs_grid  = prep_obs28[:,i,j]
                evap_obs_grid  = evap_obs28[:,i,j]
                year_mean      = np.zeros(6)
                lcc_data       = np.zeros(4)
                lcc_data[0]    = np.nanmean(file_new28[ 1: 8,int(land_main_type),i,j])
                lcc_data[1]    = np.nanmean(file_new28[ 8:15,int(land_main_type),i,j])
                lcc_data[2]    = np.nanmean(file_new28[15:22,int(land_main_type),i,j])
                lcc_data[3]    = np.nanmean(file_new28[22:29,int(land_main_type),i,j])

                def find_lcc_year(lcc_year_type):
                # calculate land cover changes
                    year_mean[lcc_year_type] = lcc_data[year_type_group[lcc_year_type,1]]-lcc_data[year_type_group[lcc_year_type,2]]
                    return year_mean

                for lcc_year_type in range(0,6):
                    # calculate the land cover changes between each two period
                    year_mean = find_lcc_year(lcc_year_type)
                year_mean_abs = list(np.abs(year_mean))

                #select the periods with max and min PFT
                year_str[i,j,int(grid_type)] = year_mean_abs.index(max(year_mean_abs))
                lcc_mean = find_year(i,j,grid_type,lcc_data)
                if abs(lcc_mean)<0.01:
                    lcc_chag[i,j,:] = np.nan
                    continue
                else:
                    lcc_chag[i,j,int(grid_type)] = lcc_mean
                    #changes of precipitation and evapotranspiration in corresponding periods
                    prep_lcc_cause[i,j],prep_total_cause[i,j],prep_lcc_data1,prep_lcc_data2,prep_tot_data1,prep_tot_data2 = LULCC_imp(i,j,prep_tot28,prep_obs28,grid_type)
                    evap_lcc_cause[i,j],evap_total_cause[i,j],evap_lcc_data1,evap_lcc_data2,evap_tot_data1,evap_tot_data2 = LULCC_imp(i,j,evap_tot28,evap_obs28,grid_type)

                    #the total impact of per percent LULCC on precipitation and evapotranspiration
                    prep_rate[i,j,int(grid_type)] = prep_total_cause[i,j]/lcc_chag[i,j,int(grid_type)]*0.01
                    evap_rate[i,j,int(grid_type)] = evap_total_cause[i,j]/lcc_chag[i,j,int(grid_type)]*0.01
                    prep_tot_plot[i,j] = prep_total_cause[i,j]/lcc_chag[i,j,int(grid_type)]*0.01
                    evap_tot_plot[i,j] = evap_total_cause[i,j]/lcc_chag[i,j,int(grid_type)]*0.01
                    WA_MSWEP_tot_rate[i,j,int(grid_type)] = (prep_total_cause[i,j]-evap_total_cause[i,j])/lcc_chag[i,j,int(grid_type)]*0.01

                    #the selected impact of per percent LULCC on precipitation and evapotranspiration
                    prep_lcc_rate[i,j,int(grid_type)] = prep_lcc_cause[i,j]/lcc_chag[i,j,int(grid_type)]*0.01
                    evap_lcc_rate[i,j,int(grid_type)] = evap_lcc_cause[i,j]/lcc_chag[i,j,int(grid_type)]*0.01
                    prep_lcc_plot[i,j] = prep_lcc_cause[i,j]/lcc_chag[i,j,int(grid_type)]*0.01
                    evap_lcc_plot[i,j] = evap_lcc_cause[i,j]/lcc_chag[i,j,int(grid_type)]*0.01
                    WA_MSWEP_lcc_rate[i,j,int(grid_type)] = (prep_lcc_cause[i,j]-evap_lcc_cause[i,j])/lcc_chag[i,j,int(grid_type)]*0.01

                    def grid_ttest(prep_tot_data1,prep_tot_data2):
                    # t-test for the precipitation and evapotranspiration changes between the corresponding periods
                      stats_lev, p_value = stats.levene(prep_tot_data1,prep_tot_data2) # Homogeneity test for data
                      if p_value>0.05:
                        tstat, pval = stats.ttest_ind(a=prep_tot_data1, b=prep_tot_data2, equal_var = True)
                      else:
                        tstat, pval = stats.ttest_ind(a=prep_tot_data1, b=prep_tot_data2, equal_var = False)
                      if pval>0.05:
                          p_grid_test = 1
                      else:
                          p_grid_test = np.nan
                      return p_grid_test,pval

                    prep_lcc_p[i,j],pval = grid_ttest(prep_lcc_data1,prep_lcc_data2)
                    evap_lcc_p[i,j],pval = grid_ttest(evap_lcc_data1,evap_lcc_data2)
                    wa_lcc_p[i,j],wa_lcc_pval[i,j] = grid_ttest(prep_lcc_data1-evap_lcc_data1,prep_lcc_data2-evap_lcc_data2)
                    wa_tot_p[i,j],wa_tot_pval[i,j] = grid_ttest(prep_tot_data1-evap_tot_data1,prep_tot_data2-evap_tot_data2)

    prep_lcc_p = np.where(prep_lcc_p==0,np.nan,prep_lcc_p)
    evap_lcc_p = np.where(evap_lcc_p==0,np.nan,evap_lcc_p)
    return year_str,lcc_chag,prep_lcc_cause,evap_lcc_cause,prep_rate,evap_rate,prep_tot_plot,evap_tot_plot,prep_lcc_rate,evap_lcc_rate,prep_lcc_plot,evap_lcc_plot,prep_lcc_p,evap_lcc_p,wa_lcc_p,wa_tot_p,WA_MSWEP_tot_rate,WA_MSWEP_lcc_rate
year_str,lcc_chag,prep_lcc_cause,evap_lcc_cause,prep_rate,evap_rate,prep_tot_plot,evap_tot_plot,prep_lcc_rate,evap_lcc_rate,prep_lcc_plot,evap_lcc_plot,prep_lcc_p,evap_lcc_p,wa_lcc_p,wa_tot_p,WA_MSWEP_tot_rate,WA_MSWEP_lcc_rate = per_percent_imp(3,prep_tot28,prep_obs28)

lcc_chag      = np.where(lcc_chag     ==0,np.nan,lcc_chag)
evap_rate     = np.where(evap_rate    ==0,np.nan,evap_rate)
evap_tot_plot = np.where(evap_tot_plot==0,np.nan,evap_tot_plot)
evap_lcc_rate = np.where(evap_lcc_rate==0,np.nan,evap_lcc_rate)
evap_lcc_plot = np.where(evap_lcc_plot==0,np.nan,evap_lcc_plot)

# assign nan for the unchanged grids and ocean
def np_where(prep_rate,prep_tot_plot,prep_lcc_rate,prep_lcc_plot,WA_MSWEP_tot_rate,WA_MSWEP_lcc_rate):
    prep_rate     = np.where(prep_rate    ==0,np.nan,prep_rate)
    prep_tot_plot = np.where(prep_tot_plot==0,np.nan,prep_tot_plot)
    prep_lcc_rate = np.where(prep_lcc_rate==0,np.nan,prep_lcc_rate)
    prep_lcc_plot = np.where(prep_lcc_plot==0,np.nan,prep_lcc_plot)
    WA_MSWEP_tot_rate = np.where(WA_MSWEP_tot_rate==0,np.nan,WA_MSWEP_tot_rate)
    WA_MSWEP_lcc_rate = np.where(WA_MSWEP_lcc_rate==0,np.nan,WA_MSWEP_lcc_rate)
    return prep_rate,prep_tot_plot,prep_lcc_rate,prep_lcc_plot,WA_MSWEP_tot_rate,WA_MSWEP_lcc_rate
prep_rate,prep_tot_plot,prep_lcc_rate,prep_lcc_plot,WA_MSWEP_tot_rate,WA_MSWEP_lcc_rate = np_where(prep_rate,prep_tot_plot,prep_lcc_rate,prep_lcc_plot,WA_MSWEP_tot_rate,WA_MSWEP_lcc_rate)

#%% t-test of each LULCC classes
def m_value(prep_lcc_rate):
    # calculate the map magnification factor
    m = np.cos(np.deg2rad(lat_lcc))
    prep_new = np.zeros((360,720))
    n_new = np.zeros((360,720))
    for lat in range(0,360):
        if np.all(np.isnan(prep_lcc_rate[lat,:]))==True:
            prep_new[lat,:] = prep_lcc_rate[lat,:]
            n_new[lat,:]    = np.nan
        else:
            prep_new[lat,:] = prep_lcc_rate[lat,:]*m[lat]
            for lon in range(0,720):
              if np.isnan(prep_lcc_rate[lat,lon]) == True:
                n_new[lat,lon] = np.nan
              else:
                n_new[lat,lon] = m[lat]

    return prep_new,n_new

# t-test for the mean change in specific LULCC classes
def t_test(prep_lcc_rate):
    prep_plot_rate = np.zeros(11)
    t_test         = np.zeros(11)
    p_test         = np.zeros(11)
    t              = np.zeros(11)
    p              = np.zeros(11)
    prep_new       = np.zeros((360,720,11))
    n_new          = np.zeros((360,720,11))
    for i in range(0,11):
        prep_new[:,:,i],n_new[:,:,i] = m_value(prep_lcc_rate[:,:,i])
        prep_plot_rate[i] = np.nansum(prep_new[:,:,i])/np.nansum(n_new[:,:,i])
        prep_unique       = np.unique(prep_lcc_rate[:,:,i]) # extract values in dramatically changed grids of each LULCC class for t-test
        t_test[i],p_test[i] = stats.ttest_1samp(prep_unique[0:-1],0)
        if p_test[i]<0.05:
            p[i] = 1
        if p_test[i]<0.01:
            p[i] = 2
        if p_test[i]<0.001:
            p[i] = 3
    p = np.where(p==0,np.nan,p)
    return prep_plot_rate,p

prep_plot_lcc,p_prep_lcc = t_test(prep_lcc_rate)
evap_plot_lcc,p_evap_lcc = t_test(evap_lcc_rate)
prep_plot_tot,p_prep_tot = t_test(prep_rate)
evap_plot_tot,p_evap_tot = t_test(evap_rate)
WA_plot_tot,p_WA_tot     = t_test(WA_MSWEP_tot_rate)
WA_plot_lcc,p_WA_lcc     = t_test(WA_MSWEP_lcc_rate)