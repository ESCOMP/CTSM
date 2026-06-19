import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
import yaml

def parse_cf(cf,lasum):
    ''' logic for parsing each conversion factor entry in the yaml file '''
    if type(cf)==str:
        cf=parse_cfstr(cf,lasum)
    elif type(cf)==list:
        cftot=1
        for f in cf:
            cftot*=parse_cf(f,lasum)
        cf=cftot
    return cf

def parse_cfstr(cf,lasum):
    '''
    replaces the str 1/lasum with the actual number 1/lasum, where lasum is the sum of the landarea, typically (ds.area*ds.landfrac).sum(dim=['lat','lon'])
    '''
    if cf=='1/lasum':
        return 1/lasum
    else:
        return float(cf)

def parse_cfs(cfs,lasum):
    '''
    translate the conversion factors inherited from the yaml to a dictionary of conversion factors, with two logic rules, lists are returned as the group product, and strings are converted to float unless they are special reserved phrases (e.g. '1/lasum')
    '''
    for v in cfs:
        cfs[v]=parse_cf(cfs[v],lasum)
    return cfs                  

def check_freq(tmp):
    '''infer if the dataset time frequency is monthly or annual'''
    # there has to be a better way to do this
    # is there alreay a function for this?
    # tricky to handle old vs. new clm history files
    nsecs_per_day=24*60*60*1e9
    if 'nbnd' in tmp.time_bounds.dims:
        dt=(tmp.time_bounds.isel(time=-1,nbnd=1)-
            tmp.time_bounds.isel(time=-1,nbnd=0))/nsecs_per_day
    else:
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
        # this is really quite slow without dask
        # I have ideas for this, but will have to come back later
        raise NotImplementedError
    else:
        ds=xr.open_mfdataset(files,combine='by_coords',preprocess=pp,
                             decode_timedelta=False, parallel = False)
        ds=ds.isel(time=slice(1,len(ds.time)))
    return ds

def plot_drifts(xs,thiscase,ncycles,nyears,units,thresholds,drifts,equils,tpct,la,lasum):
    plt.figure(figsize=[16,12])
    for j,v in enumerate(xs):
        plt.subplot(3,3,j+1) #hard-coded, should probably revise
        if 'gridded' not in v:
            # this is one type of plot, for global variables
            for i in range(ncycles):
                plt.plot(range(nyears),xs[v].isel(time=np.arange(nyears)+i*nyears),label='cycle_'+str(i).zfill(3))
            plt.ylabel(v+' ['+units[v]+']')
            if j==5:
                plt.legend()
            dstr=str(np.round(drifts[v],3))
            tstr='drift='+dstr+units[v]+'/yr'
            gl='><'
            if v in thresholds:
                tstr+=gl[int(equils[v])]
                tstr+=str(thresholds[v])
            else:
                tstr+=' [not evaluated]'
        else:
            # this is an alternative plot for gridded landarea disequilibrium
            x=xs[v]
            if v in thresholds:
                thresh=thresholds[v]
            else:
                thresh=1
            diseq=abs(x-x.shift(time=nyears))/nyears>thresh
            pct=100*(la*diseq).sum(dim=['lat','lon'])/lasum
            ix=np.arange(len(x.time))>=nyears
            pct.where(ix).plot()
            ystr=(r'abs($\Delta$'+v.split('_')[0]+')>'+
                  str(thresh)+units[v]+'/yr'+'\n[% landarea]')
            plt.ylabel(ystr)
            plt.ylim([0,100])
            plt.xlabel('')
            dstr=str(np.round(drifts[v],1))
            tstr=dstr+'%'
            if v in equils:
                tstr+=gl[int(equils[v])]
                tstr+=str(tpct)
            else:
                tstr+=' [not evaluated]'
        if v in equils:
            if not equils[v]:
                tstr='FAILED: '+tstr

        plt.title(tstr)
    plt.subplots_adjust(hspace=0.2,wspace=0.3)
    plt.savefig(thiscase+'.png',dpi=300,bbox_inches='tight')



def main():
    cfile = sys.argv[1]
    config = yaml.safe_load(open(cfile))
    thiscase = config['case']
    d = config['hist_dir']
    files = sorted(glob.glob(d+'/*.h0*'))
    print(d)
    print(len(files))
    if len(files) < 1:
        print('no files found')
        print('matchstr: '+d+'/*.h0*')
        sys.exit(1)
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

        #abbrev dict
        stocks={'TEC':'TOTECOSYSC',
                'TSC':'TOTSOMC',
                'TVC':'TOTVEGC'}
        
        #evaluate global drifts
        drifts={}
        xs={}
        for v in cfs:
            cf=cfs[v]
            if 'gridded' in v:
                vlong=stocks[v.split('_')[0]]
                x=cf*ds[vlong]
            else:
                x=cf*(la*ds[v]).sum(dim=['lat','lon'])
            xs[v]=x
            drift=abs(x.isel(time=slice(y1,y2)).mean(dim='time')-
                      x.isel(time=slice(y0,y1)).mean(dim='time'))/nyears                
            if 'gridded' in v:
                if v in thresholds:
                    thresh=thresholds[v]
                else:
                    thresh=1
                pct=100*(la*(drift>thresh)).sum(dim=['lat','lon'])/lasum
                drifts[v]=pct.values
            else:
                drifts[v]=drift.values

        equils={}
        for v in thresholds:
            if 'gridded' in v:
                equils[v]=drifts[v]<config['pct_landarea']
            else:
                equils[v]=drifts[v]<thresholds[v]

        failed=False
        fails=[]
        for v in equils:
            if 'gridded' in v:
                print(v,drifts[v],config['pct_landarea'],equils[v])
            else:
                print(v,drifts[v],thresholds[v],equils[v])
            if ~equils[v]:
                failed=True
                fails.append(v)
        if failed:
            # legacy log behavior :)
            print(("FATAL: Your simulation is not in equilibrium, "+
                   "8 hours have been deducted from your PTO bank, try again"))
        else:
            print("Congratulations! Your simulation is in equilibrium")

        if 'pct_landarea' in config:
            pct=config['pct_landarea'];
        else:
            pct=np.nan
        plot_drifts(xs,thiscase,ncycles,nyears,units,thresholds,drifts,equils,
                    pct,la,lasum)
        sys.exit(11*int(failed))



    
if __name__ == '__main__':
    main()
