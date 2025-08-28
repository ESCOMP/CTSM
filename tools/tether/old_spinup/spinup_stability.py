import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob
import sys

thiscase='I1850Clm50Bgc.CPLHIST.pAD'
d='/glade/derecho/scratch/djk2120/archive/'+thiscase+'/lnd/hist/'
files=sorted(glob.glob(d+'*.h0.*'))
if len(files)<1:
    failed=True
else:
    tmp=xr.open_dataset(files[0],decode_timedelta=True)
    la=tmp.area*tmp.landfrac
    lasum=la.sum().values

    musts=['TOTECOSYSC']
    
    drift_thresholds={'TOTECOSYSC':0.02,
                      'TOTSOMC':0.02,
                      'TOTVEGC':0.02,
                      'TLAI':0.02,
                      'GPP':0.02,
                      'TWS':1,
                      'H2OSNO':1}
    
    cfs={'TOTECOSYSC':1e-9,
         'TOTSOMC':1e-9,
         'TOTVEGC':1e-9,
         'TLAI':1/lasum,
         'GPP':1e-9*24*60*60*365,
         'TWS':1/lasum,
         'H2OSNO':1/lasum}
    
    units={'TOTECOSYSC':'PgC/yr',
           'TOTSOMC':'PgC/yr',
           'TOTVEGC':'PgC/yr',
           'TLAI':'m2/m2/yr',
           'GPP':'PgC/yr',
           'TWS':'mm/yr',
           'H2OSNO':'mm/yr'}
    
    #infer if history is monthly or annual
    nsecs_per_day=24*60*60*1e9
    dt=(tmp.time_bounds.isel(time=-1,hist_interval=1)-
        tmp.time_bounds.isel(time=-1,hist_interval=0))/nsecs_per_day
    if dt<40:
        freq='monthly'
        t=xr.DataArray(xr.date_range('1999',freq='MS',periods=12),dims='time')
        dpm=xr.DataArray(t['time.daysinmonth'].values,dims='time')
    else:
        freq='annual'
        
    #open dataset and process to an orderly annual dataset
    dvs=[v for v in units]
    def pp(ds):
        return ds[dvs]
    if freq=='monthly':
        #splitting files up by year
        #too slow to have one big open_mfdataset
        fsets=[files[i:i + 12] for i in range(0, len(files), 12)]
        dsets=[]
        for i,fset in enumerate(fsets):
            ds=xr.open_mfdataset(fset,combine='by_coords',preprocess=pp,decode_timedelta=False)
            dsets.append((dpm*ds).sum(dim='time')/365)
            ds=xr.concat(dsets,dim='time')
    else:
        ds=xr.open_mfdataset(files,combine='by_coords',preprocess=pp,decode_timedelta=False)
        ds=ds.isel(time=slice(1,len(ds.time)))
        
    #set up year variables
    nyears=20  #years per met forcing cycle
    ncycles=int(len(ds.time)/nyears)
    y2=nyears*ncycles
    y1=y2-nyears
    y0=y1-nyears
    
    #evaluate global drifts
    equils={}
    drifts={}
    for v in cfs:
        cf=cfs[v]
        x=cf*(la*ds[v]).sum(dim=['lat','lon'])
        drift=abs(x.isel(time=slice(y1,y2)).mean()-x.isel(time=slice(y0,y1)).mean())/nyears
        thresh=drift_thresholds[v]
        drifts[v]=drift.values
        equils[v]=drift.values<thresh
        
    #evaluate gridded TEC disequilibrium    
    v='TOTECOSYSC'
    x=ds[v]
    dx=abs(x.isel(time=slice(y0,y1)).mean(dim='time')-x.isel(time=slice(y1,y2)).mean(dim='time'))
    pct=100*(la*(dx>nyears)).sum(dim=['lat','lon'])/lasum
    pstr=str(np.round(pct.values,1))
    equils['TEC_gridded']=pct.values<3

    #plot results
    plt.figure(figsize=[12,14])
    for j,v in enumerate(cfs):
        plt.subplot(4,2,j+2)
        cf=cfs[v]
        x=cf*(la*ds[v]).sum(dim=['lat','lon'])
        for i in range(ncycles):
            plt.plot(range(nyears),x.isel(time=np.arange(nyears)+i*nyears),label='cycle_'+str(i).zfill(3))
            plt.legend()
            plt.ylabel(v+ ' ['+units[v]+']')
        if j>5:
            plt.xlabel('year')
        dstr=str(np.round(drifts[v],4))
        thresh=drift_thresholds[v]
        if equils[v]:
            plt.title('drift tolerable: '+dstr+'<'+str(thresh)+' '+units[v])
        else:
            plt.title('drift too high: '+dstr+'>'+str(thresh)+' '+units[v])

    plt.subplot(4,2,1)
    x=ds['TOTECOSYSC']
    diseq=abs(x-x.shift(time=20))>nyears
    pct=100*(la*diseq).sum(dim=['lat','lon'])/lasum
    ix=np.arange(len(ds.time))>=nyears
    pct.where(ix).plot()
    plt.ylabel('pct landarea in TEC diseq')
    plt.xlabel('')
    plt.title(pstr+'% diseq final segment avg')
    
    plt.subplots_adjust(hspace=0.3)
    plt.savefig(thiscase+'_spinup_stability.png',dpi=300,bbox_inches='tight')
    
    with open('tmp.txt','w') as f:
        f.write(str(drifts['TOTVEGC'])+' '+units['TOTVEGC'])

    failed=False
    for v in musts:
        if ~equils[v]:
            failed=True

sys.exit(int(failed))
