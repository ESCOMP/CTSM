import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
import yaml

def parse_cf(cf,lasum):
    if type(cf)==str:
        cf=parse_cfstr(cf,lasum)
    elif type(cf)==list:
        cftot=1
        for f in cf:
            cftot*=parse_cf(f,lasum)
        cf=cftot
    return cf

def parse_cfstr(cf,lasum):
    if cf=='1/lasum':
        return 1/lasum
    else:
        return float(cf)

def parse_cfs(cfs,lasum):
    for v in cfs:
        cfs[v]=parse_cf(cfs[v],lasum)
    return cfs                  

def check_freq(tmp):
    #infer if history is monthly or annual
    nsecs_per_day=24*60*60*1e9
    dt=(tmp.time_bounds.isel(time=-1,hist_interval=1)-
        tmp.time_bounds.isel(time=-1,hist_interval=0))/nsecs_per_day
    if dt<40:
        freq='monthly'
    else:
        freq='annual'
    return freq

def get_ds(files,freq,dvs):        
    #open dataset and process to an orderly annual dataset
    def pp(ds):
        return ds[dvs]
    if freq=='monthly':
        t=xr.DataArray(xr.date_range('1999',freq='MS',periods=12),dims='time')
        dpm=xr.DataArray(t['time.daysinmonth'].values,dims='time')
        #splitting files up by year
        #too slow to have one big open_mfdataset
        fsets=[files[i:i + 12] for i in range(0, len(files), 12)]
        dsets=[]
        for i,fset in enumerate(fsets):
            ds=xr.open_mfdataset(fset,combine='by_coords',preprocess=pp,
                                 decode_timedelta=False)
            dsets.append((dpm*ds).sum(dim='time')/365)
            ds=xr.concat(dsets,dim='time')
    else:
        ds=xr.open_mfdataset(files,combine='by_coords',preprocess=pp,
                             decode_timedelta=False)
        ds=ds.isel(time=slice(1,len(ds.time)))
    return ds

def main():

    config = yaml.safe_load(open("config.yml"))
    thiscase = config['case']
    d = config['hist_dir']
    files = sorted(glob.glob(d+'/*.h0.*'))
    if len(files) < 1:
        failed = True
    else:
        #import config and parse conversion factors
        tmp = xr.open_dataset(files[0],decode_timedelta=True)
        la = tmp.area*tmp.landfrac
        lasum = la.sum().values
        thresholds=config['thresholds']
        units=config['units']
        cfs=parse_cfs(config['cfs'],lasum)    
    
        freq=check_freq(tmp)
        dvs=config['data_vars']
        ds=get_ds(files,freq,dvs)

        #set up year variables
        nyears=config['cycle_years']  #years per met forcing cycle
        ncycles=int(len(ds.time)/nyears)
        y2=nyears*ncycles
        y1=y2-nyears
        y0=y1-nyears

        #evaluate global drifts
        equils={}
        drifts={}
        for v in cfs:
            cf=cfs[v]
            if v=='TEC_gridded':
                x=cf*ds.TOTECOSYSC
            else:
                x=cf*(la*ds[v]).sum(dim=['lat','lon'])

            drift=abs(x.isel(time=slice(y1,y2)).mean(dim='time')-
                      x.isel(time=slice(y0,y1)).mean(dim='time'))/nyears
                
            if v=='TEC_gridded':
                pct=100*(la*(drift>thresholds[v])).sum(dim=['lat','lon'])/lasum
                pstr=str(np.round(pct.values,1))
                drifts[v]=pct.values
                equils[v]=pct.values<3
            else:
                drifts[v]=drift.values
                equils[v]=drift.values<thresholds[v]


        
        print(drifts)
        



    
if __name__ == '__main__':
    main()
