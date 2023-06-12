#--------------------------------------------------------------------------------------
#This python code is used for generating Figure 5 in 
#"Interannual fires as a source for subarctic summer decadal climate variability
# mediated by permafrost thawing" by Ji-Eun Kim et al.,
# npj Climate and Atmospheric Science, 2023
# Author email jieunkim(at)pusan(dot)ac(dot)kr
#--------------------------------------------------------------------------------------


import glob
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 17})


var = 'TS'
ens = '????.0??'
yyyy = '20?0'
runname1 = 'JulD'+yyyy+'.'+ens+'.cmip_soil'
runname2 = 'JulD'+yyyy+'.'+ens+'.smbb_soil'
mm = '07'

def rd_monthly(var,ens,runname,yyyy,mm): #Reading monthly data
    modeldt = 'cam.h0.'
    #-------------------------------------------------------------------------------------------------
    direns = runname+'/archive/b.e21.BHISTsmbb.f09_g17.LE2-'+ens+'/'
    files = sorted(glob.glob(direns+domain+'/hist/b.e21.BHISTsmbb.f09_g17.LE2-'+ens+'.'+modeldt+yyyy+'-'+mm+'.nc'))
    ensemble = range(len(files))
    
    exceptcv = [var,'time','lat','lon','lev','levsoi','levgrnd','hyam','hybm','P0']
    def process_coords(ds, except_coord_vars=exceptcv):
        coord_vars = []
        for v in np.array(ds.coords):
            if not v in except_coord_vars:
                coord_vars += [v]
        for v in np.array(ds.data_vars):
            if not v in except_coord_vars:
                coord_vars += [v]
        return ds.drop(coord_vars)
    
    ds = xr.open_mfdataset(files, combine='nested',
                           concat_dim = [[*ensemble]],
                           preprocess=process_coords,
                           decode_cf=True, decode_times=True)
 
    ds = ds.rename({'concat_dim' : 'ensemble'})
    return ds[var]
 
ds1 = rd_monthly(var,ens,runname1,yyyy,mm)
yx1 = ds1.mean(['ensemble','time'])
ds2 = rd_monthly(var,ens,runname2,yyyy,mm)
yx2 = ds2.mean(['ensemble','time'])
dff = yx1 - yx2

#For python for map plotting, see fig1.py or fig4.py.
