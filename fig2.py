#--------------------------------------------------------------------------------------
#This python code is used for generating Figure 2 in 
#"Interannual fires as a source for subarctic summer decadal climate variability
# mediated by permafrost thawing" by Ji-Eun Kim et al.,
# npj Climate and Atmospheric Science, 2023
# Author email jieunkim(at)pusan(dot)ac(dot)kr
#--------------------------------------------------------------------------------------


import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from matplotlib.ticker import MultipleLocator
plt.rcParams.update({'font.size': 12})



# For significance calculations.. random sampling and save & call netcdf files-------------
def naming(var,nyear,MMM,zz=False):
    MMMw = MMM 
    if MMM == '1M' and nyear == 1 and 'SOI' in var: MMMw = '1M.tz'
    if zz: MMMw = MMM+'.zz'

    if os.path.isfile('random001.'+var+'.'+str(nyear)+'yearave.'+MMMw+'.nc'):
        fns = sorted(glob.glob('random*.'+var+'.'+str(nyear)+'yearave.'+MMMw+'.nc'))
        pp = fns[-1].find('random')
        num = fns[-1][pp+6:pp+9]
        nextnum = str(int(num)+1).zfill(3)
    
    if not os.path.isfile('random001.'+var+'.'+str(nyear)+'yearave.'+MMMw+'.nc'):
        nextnum = '001'
    return nextnum

def random_ens(number): #sampling random 50 vs 50 ensembles
    ensSelect =  ['1001.001', '1021.002', '1041.003', '1061.004', '1081.005',\
                  '1101.006', '1121.007', '1141.008', '1161.009', '1181.010',\
                  '1231.001', '1231.002', '1231.003', '1231.004', '1231.005',\
                  '1231.006', '1231.007', '1231.008', '1231.009', '1231.010',\
                  '1251.001', '1251.002', '1251.003', '1251.004', '1251.005',\
                  '1251.006', '1251.007', '1251.008', '1251.009', '1251.010',\
                  '1281.001', '1281.002', '1281.003', '1281.004', '1281.005',\
                  '1281.006', '1281.007', '1281.008', '1281.009', '1281.010',\
                  '1301.001', '1301.002', '1301.003', '1301.004', '1301.005',\
                  '1301.006', '1301.007', '1301.008', '1301.009', '1301.010',\
                  'smbb1011.001', 'smbb1031.002', 'smbb1051.003', 'smbb1071.004', 'smbb1091.005',\
                  'smbb1111.006', 'smbb1131.007', 'smbb1151.008', 'smbb1171.009', 'smbb1191.010',\
                  'smbb1231.011', 'smbb1231.012', 'smbb1231.013', 'smbb1231.014', 'smbb1231.015',\
                  'smbb1231.016', 'smbb1231.017', 'smbb1231.018', 'smbb1231.019', 'smbb1231.020',\
                  'smbb1251.011', 'smbb1251.012', 'smbb1251.013', 'smbb1251.014', 'smbb1251.015',\
                  'smbb1251.016', 'smbb1251.017', 'smbb1251.018', 'smbb1251.019', 'smbb1251.020',\
                  'smbb1281.011', 'smbb1281.012', 'smbb1281.013', 'smbb1281.014', 'smbb1281.015',\
                  'smbb1281.016', 'smbb1281.017', 'smbb1281.018', 'smbb1281.019', 'smbb1281.020',\
                  'smbb1301.011', 'smbb1301.012', 'smbb1301.013', 'smbb1301.014', 'smbb1301.015',\
                  'smbb1301.016', 'smbb1301.017', 'smbb1301.018', 'smbb1301.019', 'smbb1301.020']

    nens = len(ensSelect)
    selected_ens = sorted(random.sample(ensSelect, number))
    for kill in selected_ens: ensSelect.remove(kill)
    unselected_ens = ensSelect
    return selected_ens, unselected_ens

def boot_single(var,nyear,MMM):  #difference calculation between random 50 vs random 50
    year1 = np.random.randint(1880, high=1990-15, size=1)[0]
    year2 = year1 + nyear - 1 

    selected_ens, unselected_ens = random_ens(50)

    if var != 'Bowen_Ratio':
        if var != 'SOILLIQ+ICE':
            #get_seletive calls monthly data of CESM2-LE.
            #Please find the data from NCAR data page.
            ds1,unit,long_name = get_selective(var,'month_1',year1,year2,selected_ens)#monthly
            ds1[var] = ds1[var].fillna(0)
            tyx1 = ds1[var].mean('ensemble')
    
            ds2,unit,long_name = get_selective(var,'month_1',year1,year2,unselected_ens)#monthly
            ds2[var] = ds2[var].fillna(0)
            tyx2 = ds2[var].mean('ensemble')
    
        if var == 'SOILLIQ+ICE':
            ds1,unit,long_name = get_selective('SOILLIQ','month_1',year1,year2,selected_ens)#monthly
            ds1['SOILLIQ'] = ds1['SOILLIQ'].fillna(0)
            tyx1liq = ds1['SOILLIQ'].mean('ensemble')
            ds1,unit,long_name = get_selective('SOILICE','month_1',year1,year2,selected_ens)#monthly
            ds1['SOILICE'] = ds1['SOILICE'].fillna(0)
            tyx1ice = ds1['SOILICE'].mean('ensemble')
            tyx1 = tyx1liq + tyx1ice

            ds2,unit,long_name = get_selective('SOILLIQ','month_1',year1,year2,unselected_ens)#monthly
            ds2['SOILLIQ'] = ds2['SOILLIQ'].fillna(0)
            tyx2liq = ds2['SOILLIQ'].mean('ensemble')
            ds2,unit,long_name = get_selective('SOILICE','month_1',year1,year2,unselected_ens)#monthly
            ds2['SOILICE'] = ds2['SOILICE'].fillna(0)
            tyx2ice = ds2['SOILICE'].mean('ensemble')
            tyx2 = tyx2liq + tyx2ice

        tyx = tyx1-tyx2
        tyx = monthly2yearly(tyx,MMM)
        yx = tyx.mean('time')

    if var == 'Bowen_Ratio':
        ds,unit,long_name = get_selective('SHFLX','month_1',year1,year2,selected_ens)#monthly
        tyx = ds['SHFLX'].mean('ensemble')
        tyx = monthly2yearly(tyx,MMM)
        smbb_sh = tyx.mean('time')

        ds,unit,long_name = get_selective('LHFLX','month_1',year1,year2,selected_ens)#monthly
        tyx = ds['LHFLX'].mean('ensemble')
        tyx = monthly2yearly(tyx,MMM)
        smbb_lh = tyx.mean('time')

        ds,unit,long_name = get_selective('SHFLX','month_1',year1,year2,unselected_ens)#monthly
        tyx = ds['SHFLX'].mean('ensemble')
        tyx = monthly2yearly(tyx,MMM)
        cmip_sh = tyx.mean('time')

        ds,unit,long_name = get_selective('LHFLX','month_1',year1,year2,unselected_ens)#monthly
        tyx = ds['LHFLX'].mean('ensemble')
        tyx = monthly2yearly(tyx,MMM)
        cmip_lh = tyx.mean('time')

        smbb_bowen = smbb_sh/smbb_lh
        cmip_bowen = cmip_sh/cmip_lh
        yx = cmip_bowen - smbb_bowen
    return yx

def boot_cyx(nboot,var,nyear,MMM):  #repeat calculations from random sampling
    for ii in np.arange(nboot):
        yx = boot_single(var,nyear,MMM)

        if ii == 0:
            cyx = yx
        if ii > 0:
            cyx = xr.concat([cyx, yx], dim='sampling')

    fname = 'random'+naming(var,nyear,MMM)+'.'+var+'.'+str(nyear)+'yearave.'+MMM+'.nc'
    cyx['sampling'] = np.arange(ii+1)
    cyx.to_netcdf(fname)
    return None

def boot_ctz(nboot,var,nyear,MMM): #repeat for soil monthly
    naming(var,nyear,MMM)

    for ii in np.arange(nboot):
        tzyx = boot_single(var,nyear,MMM)
        tz = tzyx#.mean(['lat','lon'])
        tz['time'] = xr.cftime_range(start='0001-01-01',end='0001-12-31',freq='MS')

        if ii == 0:
            ctz = tz
        if ii > 0:
            ctz = xr.concat([ctz, tz], dim='sampling')

    fname = 'random'+naming(var,nyear,MMM)+'.'+var+'.'+str(nyear)+'yearave.'+MMM+'.tz.nc'
    ctz.to_netcdf(fname)
    return None


def get_std_boot(var,nyear,MMM,region): #for significanc3
    fns = glob.glob('random001.'+var+'.'+str(nyear)+'yearave.'+MMM+'.tz.nc')
    ds = xr.open_mfdataset(fns, concat_dim='sampling')

    ds, rangeStr = subsample_region(ds, region, domain='lnd')
    ds = ds.mean(['lat','lon'])
    ds_sd = ds.std('sampling')
    tz = ds_sd[var].values
    return tz

def tz_1var(var,y1,y2,diff=True): #get values as function of time and depth
    region = '50N-70N'
    
    if var != 'SOILLIQ+ICE':
       tzyx,tz = rd_diff_tyx_tt(var,region,y1,y2,diff=diff) #monthly
    if var == 'SOILLIQ+ICE':
       tzyx1,tz1 = rd_diff_tyx_tt('SOILLIQ',region,y1,y2,diff=diff)
       tzyx2,tz2 = rd_diff_tyx_tt('SOILICE',region,y1,y2,diff=diff)
       tz = tz1 + tz2
    
    levgrnd = [0.01, 0.04, 0.09, 0.16, 0.26, 0.4, 0.58, 0.8, 1.06, 1.36, 1.7, 
               2.08, 2.5, 2.99, 3.58, 4.27, 5.06, 5.95, 6.94, 8.03, 9.795, 13.32777, 
               19.48313, 28.87072, 41.99844]
    
    
    time = np.arange(len(tz.values[:,0]))/12 + y1
    cmap = plt.cm.RdBu_r
    
    zsoil = levgrnd[:20]
    lev = (np.arange(11)-5) * 0.3
    
    #Significance--------------------------------------------------
    std = get_std_boot(var,1,'1M',region) #12month data
    tztmp = tz.values 
    hatch = np.zeros_like(tztmp) - 1
    std_tz = np.transpose(np.tile(std.T, y2-y1+1))
    hatch[np.where(abs(tztmp) < 1.96*std_tz)] = 1 #95%
    
    return time,zsoil,tztmp.T,hatch.T, tz 

def subp_tz(fig,ax,ii,var,vartext,time,zsoil,zt,std_zt,diff=True): #plotting 1 sub figure
    cmap = plt.cm.RdBu_r
    levgrnd = [0.01, 0.04, 0.09, 0.16, 0.26, 0.4, 0.58, 0.8, 1.06, 1.36, 1.7, 
               2.08, 2.5, 2.99, 3.58, 4.27, 5.06, 5.95, 6.94, 8.03, 9.795, 13.32777, 
               19.48313, 28.87072, 41.99844]
    zsoil = levgrnd[:20]
    if     diff: lev = (np.arange(11)-5) * 0.25
    if not diff: lev = (np.arange(11)) * 0.5

    cf = ax.contourf(time,zsoil,zt,levels=lev,cmap=cmap,extend='both')
    cbar =fig.colorbar(cf, orientation='vertical',fraction=0.020,aspect=20,pad=0.04, ax=ax)
    cf2 = ax.contourf(time,zsoil,std_zt,levels=[-2,0,2],colors='none',hatches=[' ',3*'.'])

    ax.set_ylim(zsoil[-1],0)
    ax.set_ylabel('Soil Depth (m)')
    ty = 7.6
    if var == 'TSOI': ty = 41
    ax.text(2017,ty,vartext, horizontalalignment='right', color='k',fontsize=14)

    if ii != 999: ax.axes.get_xaxis().set_ticklabels([])

    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_major_locator(MultipleLocator(5))
    return None


def subp_monthly(ax,tt,bigmonth, annual=False):   #plotting time series
    cc = ['darkblue',  'b','cornflowerblue','cyan','lime','green','hotpink',  'r','orange','saddlebrown','dimgray','darkviolet']
    time = tt.time
    if not annual:
        yyyy = time.dt.year
        mm = time.dt.month
        yyyy = yyyy + (mm-1)/12.

        ax.fill_between(yyyy,y1=tt.where(tt>0),y2=np.zeros(len(yyyy)),alpha=0.2,color='r')
        ax.fill_between(yyyy,y1=tt.where(tt<0),y2=np.zeros(len(yyyy)),alpha=0.2,color='b')

        for mi in range(1,13):
            ttpik = tt.where(tt.time.dt.month == mi) 
            ms = 1.5
            if mi == 7: ms = 4
            ax.plot(yyyy,ttpik,'o',markersize=ms,color=cc[mi-1],label=str(mi))
    return None

def setting(ax,ii, ylabel,y1,y2):
    ax.set_frame_on(False)
    xmin, xmax = ax.get_xaxis().get_view_interval()
    xmin, xmax = y1,y2
    ymin, ymax = ax.get_yaxis().get_view_interval()
    rr = 0.92
    ax.get_yaxis().tick_left()
    ax.add_artist(lines.Line2D((xmin, xmin), (ymin*rr, ymax*rr), color='black', linewidth=1))
    ax.add_artist(lines.Line2D((xmin, xmax), (0,0), color='black', linewidth=1))
    if ii != 5: ax.axes.get_xaxis().set_visible(False)
    ax.set_ylabel(ylabel)
    ax.set_xlim(y1,y2)
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    return None


def fig2_main():
    y1, y2 = 1995, 2022

    aa = 0.055
    rr = 0.941

    left,btm1,btm2,btm3,btm4,btm5,btm6,width,height1,height2 = \
    0.15,0.8+aa,0.6+aa,0.4+aa,0.20+aa,0.2-0.115+aa,0.2-0.215+aa,0.7,0.185,0.10

    fig, axb = plt.subplots(1,1, figsize=(6,8))
    plt.subplots_adjust(left=left, bottom=btm6, right=left+width*rr, top=0.92+aa)
    
    axb.xaxis.set_visible(False)
    axb.yaxis.set_visible(False)

    axb.set_xlim(y1,y2)
    axb.axvspan(y1,2000,color='grey',alpha=0.05,lw=0)
    axb.axvspan(2000,2005,color='grey',alpha=0.15,lw=0)
    axb.axvspan(2005,2015,color='grey',alpha=0.05,lw=0)
    axb.axvspan(2010,2015,color='grey',alpha=0.15,lw=0)
    axb.axvspan(2015,2020,color='grey',alpha=0.05,lw=0)
    axb.axvspan(2020,2022,color='grey',alpha=0.15,lw=0)

    ax1 = fig.add_axes([left, btm1, width*rr, height2])
    ax2 = fig.add_axes([left, btm2, width, height1])
    ax3 = fig.add_axes([left, btm3, width, height1])
    ax4 = fig.add_axes([left, btm4, width, height1])
    ax5 = fig.add_axes([left, btm5, width*rr, height2])
    ax6 = fig.add_axes([left, btm6, width*rr, height2])
    
    ax_list = [ax1,ax2,ax3,ax4,ax5,ax6]
    
    ii = 0

    vartext_list = [' ','b. Soil Ice','c. Soil Liquid','d. Ice+Liquid','','']
    
    ax = ax_list[ii] #------------------------1
    smbb_tyx,cmip_tyx, smbb_tt,cmip_tt = rd_bmb_sum('BB','30N-90N',y1,y2)
    bmb_tt = (cmip_tt - smbb_tt)
    subp_monthly(ax,bmb_tt*1e-28,7)
    setting(ax,ii, 'a. BBA',y1,y2)
    ii = ii +1

    ax = ax_list[ii] #------------------------2
    var = 'SOILICE'
    time,zsoil,zt,std,tz = tz_1var(var,y1,y2-1)
    subp_tz(fig,ax,ii,var,vartext_list[ii],time,zsoil,zt,std)
    ii = ii +1
    
    ax = ax_list[ii] #------------------------3
    var = 'SOILLIQ'
    time,zsoil,zt,std,tz = tz_1var(var,y1,y2-1)
    subp_tz(fig,ax,ii,var,vartext_list[ii],time,zsoil,zt,std)
    ii = ii +1
    
    ax = ax_list[ii] #------------------------4
    var = 'SOILLIQ+ICE'
    time,zsoil,zt,std,tz = tz_1var(var,y1,y2-1)
    subp_tz(fig,ax,ii,var,vartext_list[ii],time,zsoil,zt,std)
    ii = ii +1
    
    ax = ax_list[ii] #------------------------5
    tyx,tt = rd_diff_tyx_tt('TS','50N-70N_land',y1,y2)
    subp_monthly(ax,tt, 7)
    setting(ax,ii, 'e. T',y1,y2)
    ii = ii +1

    ax = ax_list[ii] #------------------------6
    tyx,soilwater_tt = rd_diff_tyx_tt('SOILWATER_10CM','50N-70N_land',y1,y2)
    subp_monthly(ax,soilwater_tt, 7)
    setting(ax,ii, 'f. Upper_Water',y1,y2)

    fig.savefig('paperfig2_multi_soil.png', transparent=True, dpi=200)
    plt.show()
    return None
fig2_main()
