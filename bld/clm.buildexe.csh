#! /bin/csh -fx

set phys = $2
if ( "$phys" == "lnd" ) then
  set bldpath = $OBJROOT
else
  set bldpath = $1
endif
cd $bldpath/$phys/obj

foreach file ( Filepath CESM_cppdefs )
   if (-f $CASEBUILD/clmconf/$file) then
      cp $CASEBUILD/clmconf/$file ./tmp_$file
   else
      echo "clm.buildexe.csh ERROR - missing $CASEBUILD/clmconf/$file"
      exit -1
   endif
   if (-f $file) then
     cmp -s ./tmp_$file $file || mv -f ./tmp_$file $file 
   else
     mv -f ./tmp_$file $file 
   endif
end

set clmdefs = "`cat CESM_cppdefs`"
$GMAKE complib -j $GMAKE_J MODEL=clm COMPLIB=$bldpath/lib/lib${phys}.a USER_CPPDEFS="$clmdefs" -f $CASETOOLS/Makefile || exit 2

wait



