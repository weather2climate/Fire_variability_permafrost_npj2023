#--------------------------------------------------------------------------------------
#This python code is used for generating Figure 1 in 
#"Interannual fires as a source for subarctic summer decadal climate variability
# mediated by permafrost thawing" by Ji-Eun Kim et al.,
# npj Climate and Atmospheric Science, 2023
# Author email jieunkim(at)pusan(dot)ac(dot)kr
#--------------------------------------------------------------------------------------

import numpy as np
import xarray as xr  
import glob
import matplotlib.pyplot as plt 
plt.rcParams.update({'font.size': 10})
import matplotlib.ticker as tic
from matplotlib.ticker import MultipleLocator
import matplotlib.lines as lines


# Composites -------------------------------------------------------------------------------------

def find_compo_years(Nneg,Npos,MMM='MJJAS'):
    if Nneg == 0 and Npos == 0:
        #Random sampling (random8mean - random4mean) for significance test
        yearpool = list(range(1997,2014))
        yearA8 = np.array(sorted(random.sample(yearpool,8)))
        for kill in yearA8: yearpool.remove(kill)
        yearA4 = np.array(sorted(random.sample(yearpool,4)))

        negyear,posyear = yearA8,yearA4
        bothyear = np.concatenate((yearA8,yearA4))

    if Nneg != 0 and Npos != 0:
        aerosol_region = '30N-90N'
        #loading prescribed aerosol amounts for the two large ensemble subgroups.
        #smbb_tt and cmip_tt are time series of them.
        #For original data, find NCAR CESM2-LE data page.
        smbb_tyx,cmip_tyx, smbb_tt,cmip_tt = rd_bmb_sum('BB',aerosol_region,y1,y2) 
        pom_tt = (cmip_tt - smbb_tt)*1e-28
    
        time = cmip_tt.time
    
        nyear = y2-y1+1
        year = np.arange(nyear)+y1

        pom_chg = np.zeros(nyear) # annual sum for Mar-Aug

        for ii in range(nyear):
            if ii > 0 : #and ii < nyear-1: 
                if MMM == 'MJJAS': pom_chg[ii] = pom_tt.sel(time=slice(str(y1+ii)+'-05-01',str(y1+ii)+'-09-30')).sum().values
  
        idx = np.argsort(pom_chg)
        negyear = year[idx[:Nneg]]
        posyear = year[idx[-Npos:]]
        bothyear = np.concatenate((year[idx[:Nneg]],year[idx[-Npos:]]))
    return negyear, posyear, bothyear

def cmip6BB_compo(var,region,MMM,Nneg=8,Npos=4):   # for Fig. 1c-f calculations
    negY, posY, bothyear = find_compo_years(Nneg,Npos,MMM=MMM) #from pdd.py
    if Nneg == 0 and Npos == 0:
        Nneg,Npos = 8,4 #Change after random sampling

    tyx,ttdiff = model_eachyear(var,region,1997,diff=True)
    lat = tyx.lat.values
    lon = tyx.lon.values
 
    tyx_cmip_neg = np.zeros((12,len(lat),len(lon)))
    tyx_diff_neg = np.zeros((12,len(lat),len(lon)))
    for yi in range(Nneg):
        #loading BBA_Smooth data and BBA_CMIP6 - BBA_Smooth from saved files.
        #For original data, find NCAR CESM2-LE data page.
        #model_eachyear calls data for each year.
        tyxsmbb,ttsmbb = model_eachyear(var,region,negY[yi],diff=None) #read SMBB
        tyxdiff,ttdiff = model_eachyear(var,region,negY[yi],diff=True) #read CMIP6-SMBB
        tyxcmip = tyxdiff + tyxsmbb
        tyx_cmip_neg = tyx_cmip_neg + tyxcmip.values
        tyx_diff_neg = tyx_diff_neg + tyxdiff.values
    tyx_cmip_neg = tyx_cmip_neg/Nneg
    tyx_diff_neg = tyx_diff_neg/Nneg

    tyx_cmip_pos = np.zeros((12,len(lat),len(lon)))
    tyx_diff_pos = np.zeros((12,len(lat),len(lon)))
    for yi in range(Npos):
        tyxsmbb,ttsmbb = model_eachyear(var,region,posY[yi],diff=None) #read SMBB
        tyxdiff,ttdiff = model_eachyear(var,region,posY[yi],diff=True) #read CMIP6-SMBB
        tyxcmip = tyxdiff + tyxsmbb
        tyx_cmip_pos = tyx_cmip_pos + tyxcmip.values
        tyx_diff_pos = tyx_diff_pos + tyxdiff.values
    tyx_cmip_pos = tyx_cmip_pos/Npos
    tyx_diff_pos = tyx_diff_pos/Npos

    if MMM == 'MJJAS': tyx = tyx_cmip_neg[4:9,:,:] - tyx_cmip_pos[4:9,:,:]

    yx = tyx.mean(axis=0)
    return lon,lat,yx


# Significance Calculations for Composites -----------------------------------------------------------------

def naming(var,MMM):
    if os.path.isfile('random001.'+var+'.8Yminus4Y.'+MMM+'.nc'):
        fns = glob.glob('random*.'+var+'.8Yminus4Y.'+MMM+'.nc')
        fns = sorted(fns)
        pp = fns[-1].find('random')
        num = fns[-1][pp+6:pp+9]
        nextnum = str(int(num)+1).zfill(3)

    if not os.path.isfile('random001.'+var+'.8Yminus4Y.'+MMM+'.nc'):
        nextnum = '001'
    return nextnum

def boot_cyx(nboot,var,MMM):
    for ii in np.arange(nboot):
        lon,lat,yx = cmip6BB_compo(var,'NH',MMM,Nneg=0,Npos=0)
        yx = xr.DataArray(yx, name=var, coords=[lat,lon], dims=['lat','lon'])
        if ii == 0:
            cyx = yx
        if ii > 0:
            cyx = xr.concat([cyx, yx], dim='sampling')

    fname = 'random'+naming(var,MMM)+'.'+var+'.8Yminus4Y.'+MMM+'.nc'
    cyx['sampling'] = np.arange(ii+1)
    cyx.to_netcdf(fname)
    return None

def get_std_boot_compo_cmip(var,MMM):
    if os.path.isfile('random001.'+var+'.8Yminus4Y.'+MMM+'.nc'):
        fns = glob.glob('random00*.'+var+'.8Yminus4Y.'+MMM+'.nc')
        ds = xr.open_mfdataset(fns, concat_dim='sampling',combine='nested')
        ds_sd = ds.std('sampling')
        ds_mean = ds.mean('sampling')
        std_yx = ds_sd[var].values
        mean_yx = ds_mean[var].values
    if not os.path.isfile('random001.'+var+'.8Yminus4Y.'+MMM+'.nc'):
        std_yx = np.zeros(192*288).reshape(192,288)
        mean_yx = np.zeros(192*288).reshape(192,288)
    return mean_yx,std_yx


# For plotting ---------------------------------------------------------------------------------------------

cc = ['darkblue',  'b','cornflowerblue','cyan','lime','green','hotpink',  'r','orange','saddlebrown','dimgray','darkviolet']

def subp_monthly(tt,bigmonth, annual=False):
    time = tt.time
    if not annual:
        yyyy = time.dt.year
        mm = time.dt.month
        yyyy = yyyy + (mm-1)/12.

        ax.fill_between(yyyy,y1=tt.where(tt>0),y2=np.zeros(len(yyyy)),alpha=0.2,color='r')
        ax.fill_between(yyyy,y1=tt.where(tt<0),y2=np.zeros(len(yyyy)),alpha=0.2,color='b')

        ax.text(2004.5, 3.9, 'Month:', fontsize=11, weight='bold')
        for mi in range(1,13):
            ttpik = tt.where(tt.time.dt.month == mi) 
            ms = 2.5
            ax.plot(yyyy,ttpik,'o',markersize=ms,color=cc[mi-1],label=str(mi))
            if ii == 0: ax.text(2005.5+mi*1.3, 3.9, str(mi), color=cc[mi-1], fontsize=11, weight='bold')
    return None

def setting(ii, ylabel,y1,y2):
    ax.set_frame_on(False)

    xmin, xmax = ax.get_xaxis().get_view_interval()
    xmin, xmax = y1,y2
    ymin, ymax = ax.get_yaxis().get_view_interval()
    rr = 0.92
    if (ii % 2) ==0:
        ax.get_yaxis().tick_left()
        ax.add_artist(lines.Line2D((xmin, xmin), (ymin*rr, ymax*rr), color='black', linewidth=1))
        ax.add_artist(lines.Line2D((xmin, xmax), (0,0), color='black', linewidth=1))
    if (ii % 2) !=0:
        ax.get_yaxis().tick_right()
        ax.yaxis.set_label_position("right")
        ax.add_artist(lines.Line2D((xmax, xmax), (ymin*rr, ymax*rr), color='black', linewidth=2))
        ax.add_artist(lines.Line2D((xmin, xmax), (0,0), color='black', linewidth=1))

    if ii != 0: ax.axes.get_xaxis().set_visible(False)
    ax.set_ylabel(ylabel)
    ax.set_xlim(y1,y2)
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    return None

def subp_map_with_value(var,lon,lat,yx):
    #### plotting####################################################
    cmap=plt.cm.RdBu_r
    if var == 'BBA' : 
        levels = np.arange(10)-4.5
        cbarlabel='10$^{11}$ mol cm$^{-2}$s$^{-1}$'
        extend = 'max' 
        ticks = levels[1::3]
    if 'AOD' in var or 'AEROD_v' in var: 
        levels = (np.arange(15)*1e-3 - 7e-3)
        cbarlabel = None
        extend = 'both' 
        ticks = levels[1::6]
    if var == 'FSNS': 
        levels = (np.arange(15) - 7)*1.5
        cbarlabel='W m$^{-2}$'
        extend = 'both' 
        ticks = levels[1::4]
    if 'CCN' in var : 
        levels = (np.arange(15) - 7)*1e10
        cbarlabel='cm$^{-3}$'
        extend = 'both' 
        ticks = levels[1::4]
    if 'TS' in var : 
        levels = (np.arange(19)-9)*0.02
        cbarlabel='\u2103'
        extend = 'both' 
        levels = (np.arange(19)-9)*0.06
        ticks = levels[1::8]


    transform=ccrs.PlateCarree()

    theta = np.linspace(0, 2*np.pi, 100) 
    center, radius = [0.5, 0.5], 0.5 
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)
    ax.set_extent([-180,180,35,90], crs=transform)

    cf = ax.contourf(lon,lat,yx,levels=levels,cmap=cmap,transform=transform,extend=extend)


    #Plotting Significance--------------------------------------------------
    mean,std = get_std_boot_compo_cmip(var,MMM)
    hatch = np.zeros_like(yxin) - 1 
    hatch[np.where(abs(yxin) < std)] = 1
    cf1 = ax.contourf(lon,lat,hatch,levels=[-2,0,2],colors='none',hatches=[' ',4*'.'],transform=transform)

    ax.coastlines()
    fig.canvas.draw()

    gl = ax.gridlines(transform, draw_labels=True,linewidth=1,color='gray',alpha=0.2)#,linestyle='--')
    gl.xlabels_top = False
    gl.xlabels_bottom = False
    gl.xlines = False
    gl.xlocator = mticker.FixedLocator([-999])
    gl.ylocator = mticker.FixedLocator([50,70])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    projx, projy = ax.projection.transform_point(342,50, transform)
    ax.text(projx,projy, ' {0} '.format(str(int(50))),ha='center',va='center',color='gray',size=9)
    projx, projy = ax.projection.transform_point(358,70, transform)
    ax.text(projx,projy, ' {0} '.format(str(int(70))),ha='center',va='center',color='gray',size=9)

    var_list = [' ','b. BBA Variability ','c. AOD','d. CCN','e. Surface SW ','f. Temperature']

    ax.set_title(var_list[ii])
    cbar = fig.colorbar(cf,orientation='horizontal',fraction=0.035,aspect=25,pad=0.06,ax=ax,ticks=ticks)
    if cbarlabel: cbar.set_label(cbarlabel) 
    return None


##--------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------
#                       MAIN
##--------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------

#Fig. 1a calculation
y1 = 1995 ; y2 = 2022
smbb_tyx,cmip_tyx, smbb_tt,cmip_tt = rd_bmb_sum('BB','30N-90N',y1,y2)  
bmb_tt = (cmip_tt - smbb_tt) 

#Fig. 1b calculation
smbb_tyx,cmip_tyx, smbb_tt,cmip_tt = rd_bmb_sum('BB','Global',1997,2014)
diff_yx = (cmip_tyx - smbb_tyx)*1e-11
bmb_sd_yx = diff_yx.resample(time='YS').mean().std('time')

#Fig. 1c-f calculation
MMM = 'MJJAS'
lon,lat,aod_yx = cmip6BB_compo('AEROD_v','NH',MMM,Nneg=8,Npos=4)
lon,lat,ccn_yx = cmip6BB_compo('CCN3','NH',   MMM,Nneg=8,Npos=4)
lon,lat,fsns_yx = cmip6BB_compo('FSNS','NH',  MMM,Nneg=8,Npos=4)
lon,lat,ts_yx = cmip6BB_compo('TS','NH',      MMM,Nneg=8,Npos=4)


#================================FIGURE=====================================
aa = 0.035 ; bb = 0.19
leftA,leftB,leftC,leftD,leftE,leftF = 0.11, aa, aa+bb, aa+bb*2, aa+bb*3, aa+bb*4
btmA,btmB,btmC = 0.64, 0.13, 0.13
widthA,widthB,heightA,heightB = 0.7625, 0.315, 0.28, 0.33
widthA,widthB,heightA,heightB = 0.82, 0.17, 0.28, 0.33

fig, axb = plt.subplots(1,1, figsize=(11,6))
plt.subplots_adjust(left=leftA, bottom=btmA, right=leftA+widthA, top=btmA+heightA)
plt.rcParams.update({'font.size': 12})

axb.xaxis.set_visible(False)
axb.yaxis.set_visible(False)
axb.set_xlim(1995,2022)
axb.axvspan(y1,2000,color='grey',alpha=0.05,lw=0)
axb.axvspan(2000,2005,color='grey',alpha=0.15,lw=0)
axb.axvspan(2005,2015,color='grey',alpha=0.05,lw=0)
axb.axvspan(2010,2015,color='grey',alpha=0.15,lw=0)
axb.axvspan(2015,2020,color='grey',alpha=0.05,lw=0)
axb.axvspan(2020,2022,color='grey',alpha=0.15,lw=0)

projection=ccrs.NorthPolarStereo(central_longitude=0)

ax1 = fig.add_axes([leftA, btmA, widthA, heightA])
ax2 = fig.add_axes([leftB, btmB, widthB, heightB], projection=projection)
ax3 = fig.add_axes([leftC, btmC, widthB, heightB], projection=projection)
ax4 = fig.add_axes([leftD, btmC, widthB, heightB], projection=projection)
ax5 = fig.add_axes([leftE, btmC, widthB, heightB], projection=projection)
ax6 = fig.add_axes([leftF, btmC, widthB, heightB], projection=projection)

ax_list = [ax1,ax2,ax3,ax4,ax5,ax6]


ii = 0

ax = ax_list[ii] #------------------------1
subp_monthly(bmb_tt*1e-28,7)
setting(ii, 'a. BBA (10$^{28}$ mol s$^{-1}$)',y1,y2) 
ii = ii +1

#--------------------------------------------5Maps
ax = ax_list[ii] #------------------------2
subp_map_with_value('BBA',lon,lat,bmb_sd_yx)
ii = ii +1

ax = ax_list[ii] #------3------------------3
subp_map_with_value('AEROD_v',lon,lat,aod_yx)
ii = ii +1

ax = ax_list[ii] #------------------------4
subp_map_with_value('CCN3',lon,lat,ccn_yx)
ii = ii +1

ax = ax_list[ii] #------------------------5
subp_map_with_value('FSNS',lon,lat,fsns_yx)
ii = ii +1

ax = ax_list[ii] #------------------------6
subp_map_with_value('TS',lon,lat,ts_yx)

plt.show()
fig.savefig('paperfig1_cmip_compo.'+MMM+'.png', transparent=True, dpi=200)
