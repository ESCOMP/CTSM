load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

begin

   plot_lat=(/61.0243,60.8231,60.8328,60.9335,60.8203,60.8760,61.0866,60.5445,61.0355,60.8803,60.6652,60.6901/)
   plot_lon=(/8.12343,7.27596,7.17561,6.41504,8.70466,7.17666,6.63028,6.51468,9.07876,7.16982,6.33738,5.96487/)
   plot_name=(/"ALP1","ALP2","ALP3","ALP4","SUB1","SUB2","SUB3","SUB4","BOR1","BOR2","BOR3","BOR4"/)

   mon = ispan(0,11,1)
   mon@units  = "month"
   month =cd_string(mon,"%N")   
   print(month)

   do p=0,11,1
      clim_output   = "/cluster/shared/noresm/inputdata/atm/datm7/fates_platform/"+plot_name(p)+"/"
      clim_input    = "/cluster/shared/noresm/inputdata/atm/datm7/atm_forcing.datm7.GSWP3.0.5d.v1.c170516/"

      system("mkdir -p "+ clim_output)
 
      xfloor=plot_lon(p)-0.005
      xceil =plot_lon(p)+0.005
      yfloor=plot_lat(p)-0.005
      yceil =plot_lat(p)+0.005
      print(yceil)
      print(xceil)

      do year=1901,2010,1 
         do m=0,11,1
               
            frsds=addfile(clim_input+"/Solar/clmforc.GSWP3.c2011.0.5x0.5.Solr."+year+"-"+month(m)+".nc","r")
            fprec=addfile(clim_input+"/Precip/clmforc.GSWP3.c2011.0.5x0.5.Prec."+year+"-"+month(m)+".nc","r")
            ftphw=addfile(clim_input+"/TPHWL/clmforc.GSWP3.c2011.0.5x0.5.TPQWL."+year+"-"+month(m)+".nc","r")
      
            rsds=frsds->FSDS(:,:,:)      
            prec=fprec->PRECTmms(:,:,:)
            psrf=ftphw->PSRF(:,:,:)
            tbot=ftphw->TBOT(:,:,:)
            wind=ftphw->WIND(:,:,:)
            qbot=ftphw->QBOT(:,:,:)
            flds=ftphw->FLDS(:,:,:)
      
            lon=frsds->LONGXY(0,:)
            lat=frsds->LATIXY(:,0)
            time=frsds->time(:)
            edgew   =frsds->EDGEW
            edgee   =frsds->EDGEE
            edges   =frsds->EDGES
            edgen   =frsds->EDGEN
            longxy  =frsds->LONGXY
            latixy   =frsds->LATIXY
      
            loni=ind_nearest_coord(plot_lon(p), lon, 0)
            latj=ind_nearest_coord(plot_lat(p), lat, 0)
            print(loni)
            print(latj)

;########## Creat forcing in netcdf format

            system("rm -f "+clim_output+"clm1pt_"+plot_name(p)+"_"+year+"-"+month(m)+".nc")    ; remove if exists
            fclim  = addfile(clim_output+"clm1pt_"+plot_name(p)+"_"+year+"-"+month(m)+".nc","c")
            ntim  = dimsizes(time)                 ; get dimension sizes                                                 
            nlat  = 1  
            nlon  = 1
            nscalar = 1

            setfileoption(fclim,"DefineMode",True)
 
            fAtt               = True            ; assign file attributes
            fAtt@case_title    = "GSWP3v1: norway "  
            fAtt@conventions   = "CF-1.0" 
            fAtt@title         = "CLM single point datm input data"  
            fAtt@history       = "Original data from GSWP3v1 data"
            fAtt@creation_date = systemfunc ("date")        
            fileattdef( fclim, fAtt )            ; copy file attributes  

            dimNames = (/"scalar","lon","lat","time"/)  
            dimSizes = (/ nscalar, nlon, nlat, -1 /) 
            dimUnlim = (/ False, False, False, True/)   
            filedimdef(fclim,dimNames,dimSizes,dimUnlim)

            filevardef(fclim, "EDGEW" ,typeof(edgew),getvardims(edgew)) 
            filevardef(fclim, "EDGEE"  ,typeof(edgee),getvardims(edgee))                           
            filevardef(fclim, "EDGES"  ,typeof(edges),getvardims(edges))                          
            filevardef(fclim, "EDGEN"  ,typeof(edgen),getvardims(edgen))                          
            filevardef(fclim, "LONGXY" ,typeof(longxy)  ,getvardims(longxy))    
            filevardef(fclim, "LATIXY"  ,typeof(latixy) ,getvardims(latixy))   
            filevardef(fclim, "time",typeof(time),getvardims(time))          
            filevardef(fclim, "FSDS",typeof(rsds),getvardims(rsds))
            filevardef(fclim, "PRECTmms",typeof(prec),getvardims(prec))
            filevardef(fclim, "TBOT",typeof(tbot),getvardims(tbot))
            filevardef(fclim, "WIND",typeof(wind),getvardims(wind))
            filevardef(fclim, "PSRF",typeof(psrf),getvardims(psrf))
            filevardef(fclim, "QBOT",typeof(qbot),getvardims(qbot))
            filevardef(fclim, "FLDS",typeof(flds),getvardims(flds))
 
            filevarattdef(fclim,"EDGEW",edgew)                           
            filevarattdef(fclim,"EDGEE",edgee)                   
            filevarattdef(fclim,"EDGES",edges)                    
            filevarattdef(fclim,"EDGEN",edgen)                    
            filevarattdef(fclim,"LONGXY",longxy)                     
            filevarattdef(fclim,"LATIXY",latixy)   
            filevarattdef(fclim,"time",time)                    
            filevarattdef(fclim,"FSDS",rsds)
            filevarattdef(fclim,"PRECTmms",prec)
            filevarattdef(fclim,"TBOT",tbot)
            filevarattdef(fclim,"WIND",wind)
            filevarattdef(fclim,"PSRF",psrf)
            filevarattdef(fclim,"QBOT",qbot)
            filevarattdef(fclim,"FLDS",flds)  
 
            setfileoption(fclim,"DefineMode",False)
 
            print(time)
            fclim->time   = (/time/)        ; "(/", "/)" syntax tells NCL to only output the data values to the predefined locations on the file.
            fclim->EDGEW  = (/xfloor/)
            fclim->EDGEE  = (/xceil/) 
            fclim->EDGES  = (/yfloor/) 
            fclim->EDGEN  = (/yceil/)
            fclim->LONGXY = (/plot_lon(p)/)
            fclim->LATIXY  = (/plot_lat(p)/)
            fclim->FSDS   = (/rsds(:,latj,loni)/) 
            fclim->PRECTmms   = (/prec(:,latj,loni)/)         ; !!!! time variable can be modified in the attributes of the variable is also copied. make sure to use (/ /) to only copy data.
            fclim->TBOT   = (/tbot(:,latj,loni)/)
            fclim->WIND   = (/wind(:,latj,loni)/)
            fclim->PSRF   = (/psrf(:,latj,loni)/)
            fclim->QBOT   = (/qbot(:,latj,loni)/)
            fclim->FLDS   = (/flds(:,latj,loni)/)

            delete(time)
            delete(rsds)
            delete(prec)
            delete(tbot)
            delete(qbot)
            delete(wind)
            delete(psrf)
            delete(flds)
    
         end do
      end do
   end do
end
  
