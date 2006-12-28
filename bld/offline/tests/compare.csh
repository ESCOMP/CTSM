#!/bin/csh -fx

# ORNL IBM SP: cheetah
#set workdir = /tmp/gpfs750a/$LOGNAME
# ORNL Cray X1: phoenix
#set workdir = /scratch/scr600c/hof
# NCAR
set workdir = /ptmp/$LOGNAME

set os  = `uname -s`
if ($os == UNICOS/mp) set os = UNICOSmp
set res = $1

set OS = $os
if ($os == Linux)  set OS = PC
if ($os == AIX)    set OS = aix
if ($OS == IRIX64) set OS = sgi
if ($OS == UNICOSmp) set OS = unicosmp

echo $os $res $OS

set outfile = ${workdir}/$res.dif

# do ncdumps

cd ${workdir}/debT${OS}${res}; ncdump -ff debT${OS}${res}.clm2.i.1998-01-01-00000.nc > debT${OS}${res}.clm2.i.1998-01-01-00000.ncdump
cd ${workdir}/debF${OS}${res}; ncdump -ff debF${OS}${res}.clm2.i.1998-01-01-00000.nc > debF${OS}${res}.clm2.i.1998-01-01-00000.ncdump
cd ${workdir}/test${OS}${res}; ncdump -ff test${OS}${res}.clm2.i.1998-01-01-00000.nc > test${OS}${res}.clm2.i.1998-01-01-00000.ncdump
cd ${workdir}/cont${OS}${res}; ncdump -ff cont${OS}${res}.clm2.i.1998-01-01-00000.nc > cont${OS}${res}.clm2.i.1998-01-01-00000.ncdump

cd ${workdir}/debT${OS}${res}/HISTOUT; ncdump -ff clm2.h0.1998-01-01-05400.nc > clm2.h0.1998-01-01-05400.ncdump
cd ${workdir}/debF${OS}${res}/HISTOUT; ncdump -ff clm2.h0.1998-01-01-05400.nc > clm2.h0.1998-01-01-05400.ncdump


cd ${workdir}/test${OS}${res}/INITIAL; ncdump -ff clm2.h0.1998-01-01-48600.nc > clm2.h0.1998-01-01-48600.ncdump
cd ${workdir}/test${OS}${res}/RESTART; ncdump -ff clm2.h0.1998-01-01-48600.nc > clm2.h0.1998-01-01-48600.ncdump
cd ${workdir}/cont${OS}${res}/HISTOUT; ncdump -ff clm2.h0.1998-01-01-48600.nc > clm2.h0.1998-01-01-48600.ncdump
												      
cd ${workdir}/test${OS}${res}/INITIAL; ncdump -ff clm2.h1.1998-01-01-75600.nc > clm2.h1.1998-01-01-75600.ncdump
cd ${workdir}/test${OS}${res}/RESTART; ncdump -ff clm2.h1.1998-01-01-75600.nc > clm2.h1.1998-01-01-75600.ncdump
cd ${workdir}/cont${OS}${res}/HISTOUT; ncdump -ff clm2.h1.1998-01-01-75600.nc > clm2.h1.1998-01-01-75600.ncdump

if ($res == DGVM) then
cd ${workdir}/debT${OS}${res}; ncdump -ff debT${OS}${res}.clm2.hv.1998-01-01-00000.nc > debT${OS}${res}.clm2.hv.1998-01-01-00000.ncdump
cd ${workdir}/debF${OS}${res}; ncdump -ff debF${OS}${res}.clm2.hv.1998-01-01-00000.nc > debF${OS}${res}.clm2.hv.1998-01-01-00000.ncdump
cd ${workdir}/test${OS}${res}; ncdump -ff test${OS}${res}.clm2.hv.1998-01-01-00000.nc > test${OS}${res}.clm2.hv.1998-01-01-00000.ncdump
cd ${workdir}/cont${OS}${res}; ncdump -ff cont${OS}${res}.clm2.hv.1998-01-01-00000.nc > cont${OS}${res}.clm2.hv.1998-01-01-00000.ncdump
endif

#===========================================================================================================

mkdir -p ${workdir}/Dif$res; cd ${workdir}/Dif$res

set dir1 = debT${OS}${res}
set dir2 = debF${OS}${res}

diff ${workdir}/$dir1/$dir1.clm2.i.1998-01-01-00000.ncdump  ${workdir}/$dir2/$dir2.clm2.i.1998-01-01-00000.ncdump  >! diff.debug.inic
if ($res == DGVM) then
diff ${workdir}/$dir1/$dir1.clm2.hv.1998-01-01-00000.ncdump ${workdir}/$dir2/$dir2.clm2.hv.1998-01-01-00000.ncdump >! diff.debug.hv
endif
diff ${workdir}/$dir1/HISTOUT/clm2.h0.1998-01-01-05400.ncdump   ${workdir}/$dir2/HISTOUT/clm2.h0.1998-01-01-05400.ncdump   >! diff.debug.h0

#===========================================================================================================

set dir1 = test${OS}${res}
set dir2 = cont${OS}${res}

diff ${workdir}/$dir1/$dir1.clm2.i.1998-01-01-00000.ncdump  ${workdir}/$dir2/$dir2.clm2.i.1998-01-01-00000.ncdump  >! diff.cont.inic
if ($res == DGVM) then
diff ${workdir}/$dir1/$dir1.clm2.hv.1998-01-01-00000.ncdump ${workdir}/$dir2/$dir2.clm2.hv.1998-01-01-00000.ncdump >! diff.cont.hv
endif

diff ${workdir}/$dir1/INITIAL/clm2.h0.1998-01-01-48600.ncdump   ${workdir}/$dir2/HISTOUT/clm2.h0.1998-01-01-48600.ncdump   >! diff.cont.h0
diff ${workdir}/$dir1/INITIAL/clm2.h1.1998-01-01-75600.ncdump   ${workdir}/$dir2/HISTOUT/clm2.h1.1998-01-01-75600.ncdump   >! diff.cont.h1

diff ${workdir}/$dir1/INITIAL/clm2.h0.1998-01-01-48600.ncdump   ${workdir}/$dir1/RESTART/clm2.h0.1998-01-01-48600.ncdump   >! diff.rest.h0
diff ${workdir}/$dir1/INITIAL/clm2.h1.1998-01-01-75600.ncdump   ${workdir}/$dir1/RESTART/clm2.h1.1998-01-01-75600.ncdump   >! diff.rest.h1

rm dif.$res; touch dif.$res
foreach i (diff*)
   echo $i
   cat $i >> dif.$res
end
