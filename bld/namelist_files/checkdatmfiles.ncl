;
; Check that the area and fraction is the same between the datm domain files
; and the equivalent CLM land frac file.
;
;  Erik Kluzek
;  Apr/10/2008
;

  print( "Check that datm domain files and CLM land frac files are consistent" );
  resolutions = (/ "4x5", "10x15", "1.9x2.5", "0.9x1.25", "0.47x0.63", "0.23x0.31", "128x256","64x128","48x96","32x64","8x16" /);
 
  space  = "      ";
  badres = 0
  badresolutions = new( (/ 1000 /), string )
  ;type   = "fatmlndfrc";    ; Can be either fatmlndfrc or fatmgrid
  type   = "fatmgrid";    ; Can be either fatmlndfrc or fatmgrid

procedure checkit( desc:string, maxdiff:numeric, res:string, omask:string, eps:numeric )
;
;
;
begin
  if ( maxdiff .gt. eps )then
     print( space+space+space+desc+" are off by more than tolerance for "+res+"_"+omask+" resolution" );
     print( space+space+space+"maximum difference = "+maxdiff );
     badresolutions(badres) = res+"_"+omask;
     badres = badres + 1
  else
     print( space+space+space+"File OK for "+desc+"!" );
  end if
end


begin

  csmdata  = getenv("CSMDATA");
  clmroot  = getenv("CLM_ROOT");
  querynml = "bld/queryDefaultNamelist.pl -silent -justvalue ";
  if ( .not. ismissing(csmdata) )then
     querynml = querynml+" -csmdata "+csmdata;
  end if
  if ( ismissing(clmroot) )then
     querynml = "../../"+querynml;
  else
     querynml = clmroot+"/models/lnd/clm*/"+querynml;
  end if

  print( "query string="+querynml )


  do i = 0, dimsizes(resolutions)-1
     res = resolutions(i);
     print( "Resolution: "+res );
     masks = (/ "gx3v7", "gx1v6", "tx1v1", "USGS" /);
     do j = 0, dimsizes(masks)-1
        omask = masks(j);
        print( "Mask: "+omask);

        querynmlres = querynml+" -res "+res+" -options mask="+omask;
        ;
        ; Get grid filename and open it
        ;
        fracfile  = systemfunc( querynmlres+" -var fatmlndfrc" );
        if ( systemfunc("test -f "+fracfile+"; echo $?" ) .ne. 0 )then
           delete( fracfile );
           continue;
        end if
        gridfile  = systemfunc( querynmlres+" -var "+type );
        if ( systemfunc("test -f "+gridfile+"; echo $?" ) .ne. 0 )then
           delete( gridfile );
           continue;
        end if
        print( space+"Use "+type+" file:       "+gridfile );
        ncg     = addfile( gridfile,  "r" );
   
        ;
        ; Get datm filename and open it
        ;
        domfile  = systemfunc( querynmlres+" -var domainfile -namelist shr_strdata_nml" );
        if ( ismissing(domfile) )then
          print( "Missing domainfile" );
          continue;
        end if
        print( querynml+" -res "+res+" -var domainfile" );
        print( space+space+"Use dom file:         "+domfile );
        if ( systemfunc("test -f "+domfile+"; echo $?" ) .ne. 0 )then
           print( "Input domfile does not exist or not found: "+domfile );
           continue;
        end if
        ncd      = addfile( domfile,    "r" );
   
        maxdiff = max( abs(ncd->yc - ncg->LATIXY) );
        checkit( "Lats ", maxdiff, res, omask, 1.e-12 );
        maxdiff = max( abs(ncd->xc - ncg->LONGXY) );
        checkit( "Longs", maxdiff, res, omask, 1.e-12 );

        if ( type .eq. "fatmlndfrc" )then
           maxdiff = max( abs(ncd->frac - ncg->LANDFRAC) );
           checkit( "Fracs", maxdiff, res, omask, 1.0e-14 );
        else
           re = 6371.22;
           r2 = re*re;
           maxdiff = max( abs(ncd->area - ncg->AREA/r2) );
           checkit( "Area", maxdiff, res, omask, 9.0e-07 );
        end if
        delete( maxdiff );
        delete( domfile );
        delete( ncd );
   
        delete( ncg );
        delete( gridfile );

     end do

     delete( res );

  end do
  if ( badres .gt. 0 )then
     print( "badresolutions = " );
     print( badresolutions(0:badres-1) );
  end if

  print( "===============================" );
  print( "Successfully went through files" );

end

