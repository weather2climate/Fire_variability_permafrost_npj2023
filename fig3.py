#--------------------------------------------------------------------------------------
#This python code is used for generating Figure 3 in 
#"Interannual fires as a source for subarctic summer decadal climate variability
# mediated by permafrost thawing" by Ji-Eun Kim et al.,
# npj Climate and Atmospheric Science, 2023
# Author email jieunkim(at)pusan(dot)ac(dot)kr
#--------------------------------------------------------------------------------------


import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 20})


# Fig. 3a-c ----------------------------------------------------------------------------
def scatter_polyfit_eachmonth(var1,var2):
    y1, y2 = 1997, 2014
    aerosol_region = '30N-90N'
    region = '50N-70N_land'
 
    # rd_bmb_sm loads aerosol data. Please see NCAR CESM2 website for reference.
    # rd_diff_tyx_tt loads difference data BBA_CMIP6-BBA_Smooth
    if var1 == 'BBA':
        smbb_tyx,cmip_tyx, smbb_tt,cmip_tt = rd_bmb_sum('BB',aerosol_region,y1,y2)  
        tt1 = (cmip_tt - smbb_tt)*1e-28
    if not var1 == 'BBA': tyx,tt1 = rd_diff_tyx_tt(var1,region,y1,y2)
    time = tt1.time
 
    tyx,tt2 = rd_diff_tyx_tt(var2,region,y1,y2)
    
    nyear = y2-y1+1
    year = np.arange(nyear)+y1
 
    tt1_chg6 = np.zeros(nyear)
    tt2_chg6 = np.zeros(nyear)
    for ii in range(nyear):
        if var1 == 'BBA':
            tt1_chg6[ii] = tt1.sel(time=slice(str(y1+ii)+'-06-01',str(y1+ii)+'-06-30')).sum().values
        if not var1 == 'BBA':
            tt1_chg6[ii] = tt1.sel(time=slice(str(y1+ii)+'-06-01',str(y1+ii)+'-06-30')).mean().values
        tt2_chg6[ii] =    tt2.sel(time=slice(str(y1+ii)+'-06-01',str(y1+ii)+'-06-30')).mean().values    
 
 
    tt1_chg7 = np.zeros(nyear)
    tt2_chg7 = np.zeros(nyear)
    for ii in range(nyear):
        if var1 == 'BBA':
            tt1_chg7[ii] = tt1.sel(time=slice(str(y1+ii)+'-07-01',str(y1+ii)+'-07-31')).sum().values
        if not var1 == 'BBA':
            tt1_chg7[ii] = tt1.sel(time=slice(str(y1+ii)+'-07-01',str(y1+ii)+'-07-31')).mean().values
        tt2_chg7[ii] =    tt2.sel(time=slice(str(y1+ii)+'-07-01',str(y1+ii)+'-07-31')).mean().values    


    tt1_chg8 = np.zeros(nyear)
    tt2_chg8 = np.zeros(nyear)
    for ii in range(nyear):
        if var1 == 'BBA':
            tt1_chg8[ii] = tt1.sel(time=slice(str(y1+ii)+'-08-01',str(y1+ii)+'-08-31')).sum().values
        if not var1 == 'BBA':
            tt1_chg8[ii] = tt1.sel(time=slice(str(y1+ii)+'-08-01',str(y1+ii)+'-08-31')).mean().values
        tt2_chg8[ii] =    tt2.sel(time=slice(str(y1+ii)+'-08-01',str(y1+ii)+'-08-31')).mean().values    


    tt1_chg678 = np.zeros(nyear)
    tt2_chg678 = np.zeros(nyear)
    for ii in range(nyear):
        if var1 == 'BBA':
            tt1_chg678[ii] = tt1.sel(time=slice(str(y1+ii)+'-06-01',str(y1+ii)+'-08-31')).sum().values
        if not var1 == 'BBA':
            tt1_chg678[ii] = tt1.sel(time=slice(str(y1+ii)+'-06-01',str(y1+ii)+'-08-31')).mean().values
        tt2_chg678[ii] =    tt2.sel(time=slice(str(y1+ii)+'-06-01',str(y1+ii)+'-08-31')).mean().values    
 
 
    if var1 == 'CDNUMC' or 'CCN' in var1 : tt1_chg6,tt1_chg7,tt1_chg8 = tt1_chg6*1e-10,tt1_chg7*1e-10,tt1_chg8*1e-10
    if var2 == 'CDNUMC' or 'CCN' in var2 : tt2_chg6,tt2_chg7,tt2_chg8 = tt2_chg6*1e-10,tt2_chg7*1e-10,tt2_chg8*1e-10
 
    fig = plt.figure(figsize=(6,6))
    plt.subplots_adjust(left = 0.2, bottom = 0.2, right = 0.95, top = 0.95)
 
    plt.scatter(tt1_chg6,tt2_chg6,label='Jun')
    plt.scatter(tt1_chg7,tt2_chg7,label='Jul')
    plt.scatter(tt1_chg8,tt2_chg8,label='Aug')
    plt.legend()
 
    # To remove outliers-----------------------------------
    tt1_chgS = np.concatenate([tt1_chg6,tt1_chg7,tt1_chg8])
    tt2_chgS = np.concatenate([tt2_chg6,tt2_chg7,tt2_chg8])
 
    tt1_sort = np.sort(tt1_chgS)
    tt2_sort = tt2_chgS[np.argsort(tt1_chgS)]
 
    mymodel = np.poly1d(np.polyfit(tt1_sort,tt2_sort,2))
 
    sigC = 3
    err = (tt2_sort - mymodel(tt1_sort))**2
 
    tt1_dead = tt1_sort[np.where(err > sigC*np.std(err))]
    tt2_dead = tt2_sort[np.where(err > sigC*np.std(err))]
    plt.scatter(tt1_dead,tt2_dead,c='white',s=5)
    tt1_sort = tt1_sort[np.where(err < sigC*np.std(err))]
    tt2_sort = tt2_sort[np.where(err < sigC*np.std(err))]
    mymodel = np.poly1d(np.polyfit(tt1_sort,tt2_sort,2))
 
    # Repeat removing outliers----------------------------
    err = (tt2_sort - mymodel(tt1_sort))**2
    tt1_dead = tt1_sort[np.where(err > sigC*np.std(err))]
    tt2_dead = tt2_sort[np.where(err > sigC*np.std(err))]
    plt.scatter(tt1_dead,tt2_dead,c='white',s=5)
    tt1_sort = tt1_sort[np.where(err < sigC*np.std(err))]
    tt2_sort = tt2_sort[np.where(err < sigC*np.std(err))]
    mymodel = np.poly1d(np.polyfit(tt1_sort,tt2_sort,2))
    
    myline = np.linspace(min(tt1_sort)*1. , max(tt1_sort)*1. , 100)
    plt.plot(myline, mymodel(myline),'k')
 
    plt.xlabel(renamevar(var1)+' ('+units(var1)+')')
    plt.ylabel(renamevar(var2)+' ('+units(var2)+')')
    plt.show()
    fig.savefig('polyfit.'+var1+'-'+var2+'.png', transparent=True, dpi=200)
    return None
scatter_polyfit_eachmonth('BBA','CCN3')
scatter_polyfit_eachmonth('CCN3','CDNUMC')
scatter_polyfit_eachmonth('CDNUMC','FSNS-FSNSC')
scatter_polyfit_eachmonth('BBA','FSNS')


# Fig. 3e ------------------------------------------------------------------------------
def scatter_polyfit_eachmonth_tendency(var1,var2):
    y1, y2 = 1997, 2014
    aerosol_region = '30N-90N'
    region = '50N-70N_land'
 
    if var1 == 'BBA':
        smbb_tyx,cmip_tyx, smbb_tt,cmip_tt = rd_bmb_sum('BB',aerosol_region,y1,y2)  
        tt1 = (cmip_tt - smbb_tt)*1e-28
    if not var1 == 'BBA': tyx,tt1 = rd_diff_tyx_tt(var1,region,y1,y2)
    time = tt1.time
 
    tyx,tt2 = rd_diff_tyx_tt(var2,region,y1,y2)
    
    nyear = y2-y1+1
    nyear = y2-y1+1 - 1
    year = np.arange(nyear)+y1
 
    tt1_chg6 = np.zeros(nyear)
    tt2_chg6 = np.zeros(nyear)
    for ii in range(nyear):
        if var1 == 'BBA':
            tt1_chg6[ii] = tt1.sel(time=slice(str(y1+ii)+'-06-01',str(y1+ii)+'-06-30')).sum().values
        if not var1 == 'BBA':
            tt1_chg6[ii] = tt1.sel(time=slice(str(y1+ii)+'-06-01',str(y1+ii)+'-06-30')).mean().values
        tt2_tmp6 =    tt2.sel(time=slice(str(y1+ii)+'-06-01',str(y1+ii)+'-06-30')).mean().values    
        tt2_tmp8 =    tt2.sel(time=slice(str(y1+ii)+'-08-01',str(y1+ii)+'-08-31')).mean().values    
        tt2_chg6[ii]= tt2_tmp8 - tt2_tmp6
 
    tt1_chg7 = np.zeros(nyear)
    tt2_chg7 = np.zeros(nyear)
    for ii in range(nyear):
        if var1 == 'BBA':
            tt1_chg7[ii] = tt1.sel(time=slice(str(y1+ii)+'-07-01',str(y1+ii)+'-07-31')).sum().values
        if not var1 == 'BBA':
            tt1_chg7[ii] = tt1.sel(time=slice(str(y1+ii)+'-07-01',str(y1+ii)+'-07-31')).mean().values
        tt2_tmp7 =    tt2.sel(time=slice(str(y1+ii)+'-07-01',str(y1+ii)+'-07-31')).mean().values    
        tt2_tmp9 =    tt2.sel(time=slice(str(y1+ii)+'-09-01',str(y1+ii)+'-09-30')).mean().values    
        tt2_chg7[ii]= tt2_tmp9 - tt2_tmp7
 
 
    tt1_chg8 = np.zeros(nyear)
    tt2_chg8 = np.zeros(nyear)
    for ii in range(nyear):
        if var1 == 'BBA':
            tt1_chg8[ii] = tt1.sel(time=slice(str(y1+ii)+'-08-01',str(y1+ii)+'-08-31')).sum().values
        if not var1 == 'BBA':
            tt1_chg8[ii] = tt1.sel(time=slice(str(y1+ii)+'-08-01',str(y1+ii)+'-08-31')).mean().values
        tt2_tmp8 =    tt2.sel(time=slice(str(y1+ii)+'-08-01',str(y1+ii)+'-08-31')).mean().values    
        tt2_tmp10 =    tt2.sel(time=slice(str(y1+ii)+'-10-01',str(y1+ii)+'-10-31')).mean().values    
        tt2_chg8[ii]= tt2_tmp10 - tt2_tmp8
 
 
    fig = plt.figure(figsize=(6,6))
  
    plt.subplots_adjust(left = 0.2, bottom = 0.2, right = 0.95, top = 0.95)
 
    plt.scatter(tt1_chg6,tt2_chg6,label='Jun->Aug')
    plt.scatter(tt1_chg7,tt2_chg7,label='Jul->Sep')
    plt.scatter(tt1_chg8,tt2_chg8,label='Aug-Oct')
    plt.legend()
 
    # To remove outliers-----------------------------------
    tt1_chgS = np.concatenate([tt1_chg6,tt1_chg7,tt1_chg8])
    tt2_chgS = np.concatenate([tt2_chg6,tt2_chg7,tt2_chg8])
 
    tt1_sort = np.sort(tt1_chgS)
    tt2_sort = tt2_chgS[np.argsort(tt1_chgS)]
 
    mymodel = np.poly1d(np.polyfit(tt1_sort,tt2_sort,2))
 
    sigC = 3
    err = (tt2_sort - mymodel(tt1_sort))**2
 
    tt1_dead = tt1_sort[np.where(err > sigC*np.std(err))]
    tt2_dead = tt2_sort[np.where(err > sigC*np.std(err))]
    plt.scatter(tt1_dead,tt2_dead,c='white',s=5)
    tt1_sort = tt1_sort[np.where(err < sigC*np.std(err))]
    tt2_sort = tt2_sort[np.where(err < sigC*np.std(err))]
    mymodel = np.poly1d(np.polyfit(tt1_sort,tt2_sort,2))
 
    # Repeat removing outliers----------------------------
    err = (tt2_sort - mymodel(tt1_sort))**2
    tt1_dead = tt1_sort[np.where(err > sigC*np.std(err))]
    tt2_dead = tt2_sort[np.where(err > sigC*np.std(err))]
    plt.scatter(tt1_dead,tt2_dead,c='white',s=5)
    tt1_sort = tt1_sort[np.where(err < sigC*np.std(err))]
    tt2_sort = tt2_sort[np.where(err < sigC*np.std(err))]
    mymodel = np.poly1d(np.polyfit(tt1_sort,tt2_sort,2))
    
    myline = np.linspace(min(tt1_sort)*1. , max(tt1_sort)*1. , 100)
    plt.plot(myline, mymodel(myline),'k')
 
    plt.xlabel(renamevar(var1)+' ('+units(var1)+')')
    plt.ylabel(renamevar(var2)+' ('+units(var2)+')')
    plt.show()
    fig.savefig('polyfit.'+var1+'-'+var2+'.png', transparent=True, dpi=200)
    return None
scatter_polyfit_eachmonth_tendency('FSNS','TOPd_SOILICE')
