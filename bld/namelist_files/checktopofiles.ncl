;
;
;  Erik Kluzek
;  Apr/10/2008
;
;  Check that the latitudes and longitudes between the grid and topo
;  files are within a reasonable tolerance.
;
begin

  print( "Check that lats/longs are the same between grid and topo files" );
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

  resolutions = (/ "4x5", "10x15", "1.9x2.5", "0.9x1.25", "0.47x0.63", "0.23x0.31", "128x256","64x128","48x96","32x64","8x16" /);

  space = "      ";
  eps   = 1.e-12;
  do i = 0, dimsizes(resolutions)-1
     res = resolutions(i);
     print( "Resolution: "+res );

     ;
     ; Get grid filename and open it
     ;
     gridfile  = systemfunc( querynml+" -res "+res+" -var fatmgrid" );
     print( space+"Use grid file:       "+gridfile );
     if ( systemfunc("test -f "+gridfile+"; echo $?" ) .ne. 0 )then
        print( "Input gridfile does not exist or not found: "+gridfile );
        continue
     end if
     ncg     = addfile( gridfile,  "r" );

     ;
     ; Get topo filename and open it
     ;
     topofile  = systemfunc( querynml+" -res "+res+" -var flndtopo" );
     if ( ismissing(topofile) )then
        print( "Missing topofile" );
        continue
     end if
     ;;;print( querynml+" -res "+res+" -var fatmtopo" );
     print( space+space+"Use topo file:         "+topofile );
     if ( systemfunc("test -f "+topofile+"; echo $?" ) .ne. 0 )then
        print( "Input topofile does not exist or not found: "+topofile );
        continue
     end if
     nct      = addfile( topofile,    "r" );

     maxdiff = max( abs(nct->LATIXY - ncg->LATIXY) );
     if ( maxdiff .gt. eps )then
        print( space+space+space+"Lats  are off by more than tolerance for "+res+" resolution" );
        print( space+space+space+"maximum difference = "+maxdiff );
     else
        print( space+space+space+"File OK for Lats!" );
     end if
     maxdiff = max( abs(nct->LONGXY - ncg->LONGXY) );
     if ( maxdiff .gt. eps )then
        print( space+space+space+"Longs are off by more than tolerance for "+res+" resolution" );
        print( space+space+space+"maximum difference = "+maxdiff );
     else
        print( space+space+space+"File OK for Longs!" );
     end if
     delete( maxdiff );
     delete( topofile );
     delete( nct );

     delete( ncg );
     delete( res );
     delete( gridfile );

  end do

  print( "===============================" );
  print( "Successfully went through files" );

end
