#--------------------------------------------------------------------------------------
#This python code is used for generating Figure 4 in 
#"Interannual fires as a source for subarctic summer decadal climate variability
# mediated by permafrost thawing" by Ji-Eun Kim et al.,
# npj Climate and Atmospheric Science, 2023
# Author email jieunkim(at)pusan(dot)ac(dot)kr
#--------------------------------------------------------------------------------------


import numpy as np
import xarray as xr  
import glob
import matplotlib.pyplot as plt 
import matplotlib.ticker as tic
plt.rcParams.update({'font.size': 12})
import matplotlib.lines as lines

from fig2 import get_std_boot

def subp_map(var,y1,y2,MMM):
    # rd_diff_tyx_tt loads difference data BBA_CMIP6-BBA_Smooth
    # Please see NCAR CESM2 website for reference.
    tyx,tt_var = rd_diff_tyx_tt(var,'Global',y1,y2)
    lon,lat = yx.lon,yx.lat
    yx = tyx.mean('time').values
    std = get_std_boot(var,y2-y1+1,MMM)
 
    #### plotting####################################################
    levels,cmap = pltVAR(var,'diff',yx)
    cmap = 'RdBu_r'
    if var == 'TS' and MMM == '07': levels = (np.arange(19)-9)*0.13

    transform=ccrs.PlateCarree()

    theta = np.linspace(0, 2*np.pi, 100) 
    center, radius = [0.5, 0.5], 0.5 
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)
    ax.set_extent([-180,180,35,90], crs=transform)

    extend = 'neither' 
    if var == 'TS' : extend = 'max'
    if 'WATER' in var: extend = 'min'
    if 'Q' in var: extend = 'both'
    if 'FSNS' in var: extend = 'both'
    if 'RELHUM' in var: extend = 'both'
    cf = ax.contourf(lon,lat,yx,levels=levels,cmap=cmap,transform=transform,extend=extend)

    #Significance--------------------------------------------------
    hatch = np.zeros_like(yx) - 1 
    hatch[np.where(abs(yx) < std)] = 1
    londum,hatchplot = fill_edge_map(lon,np.asarray(hatch))
    cf1 = ax.contourf(lon,lat,hatchplot,levels=[-2,0,2],colors='none',hatches=[' ',4*'.'],transform=transform)
    ax.coastlines()
    fig.canvas.draw()

    gl = ax.gridlines(transform, draw_labels=False,linewidth=1,color='gray',alpha=0.2)
    gl.xlocator = mticker.FixedLocator([-999])
    gl.ylocator = mticker.FixedLocator([50,70])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlines = False
    gl.xlabels_top = False
    gl.xlabels_bottom = False
    gl.ylabels = False

    projx, projy = ax.projection.transform_point(342,50, transform)
    ax.text(projx,projy, ' {0} '.format(str(int(50))),ha='center',va='center',color='gray',size=9)
    projx, projy = ax.projection.transform_point(358,70, transform)
    ax.text(projx,projy, ' {0} '.format(str(int(70))),ha='center',va='center',color='gray',size=9)

    var_list = ['a. Temperature','b. Upper Soil Water','c. Surface SW','d. Relative Humidity','e. Low Cloud Fraction']
    unit_list = ['   \u2103','kg/m$^{2}$','W/m$^{2}$','    %','']

    ticks = levels[1::4]
    if var == 'Q': ticks = levels[1::6]
    if var == 'RELHUM': ticks = levels[1::6]

    ax.set_title(var_list[ii])
    cbar = fig.colorbar(cf, orientation='horizontal',fraction=0.035,aspect=25,pad=0.06, ax=ax,\
                        ticks=ticks)
    ax.text(31, 23, unit_list[ii], transform=transform)
    return None

#================================FIGURE MAIN=====================================
aa = -0.2
btm,diameter = 0.17, 0.60
width = 0.2
l1,l2,l3,l4,l5 = aa, aa+width, aa+width*2, aa+width*3, aa+width*4

fig, ax = plt.subplots(1,1, figsize=(14,4.2))
plt.subplots_adjust(left = 0.29, bottom = 0.297, right = 0.69, top = 0.697)
plt.axis('off')

projection=projection=ccrs.NorthPolarStereo(central_longitude=0)

ax1 = fig.add_axes([l1, btm, diameter, diameter], projection=projection)
ax2 = fig.add_axes([l2, btm, diameter, diameter], projection=projection)
ax3 = fig.add_axes([l3, btm, diameter, diameter], projection=projection)
ax4 = fig.add_axes([l4, btm, diameter, diameter], projection=projection)
ax5 = fig.add_axes([l5, btm, diameter, diameter], projection=projection)


#-----------------------------------------Maps
y1, y2 = 2000,2014 
ax_list = [ax1,ax2,ax3,ax4,ax5]
ii = 0

ax = ax_list[ii] #------------------------1
subp_map('TS',y1,y2,'07')
ii = ii +1

ax = ax_list[ii] #------------------------2
subp_map('SOILWATER_10CM',y1,y2,'07')
ii = ii +1

ax = ax_list[ii] #------------------------3
subp_map('FSNS',y1,y2,'07')
ii = ii +1

ax = ax_list[ii] #------------------------4
subp_map('RELHUM',y1,y2,'07')
ii = ii +1

ax = ax_list[ii] #------------------------5
subp_map('CLDLOW',y1,y2,'07')

plt.show()
fig.savefig('paperfig4_decadal_july.png', transparent=True, dpi=200)
