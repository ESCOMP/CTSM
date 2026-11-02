"""
Create a 12-month climatology of the aerosol deposition fields from daily CMEPS coupler history
(CPLHIST) files, for use by DATM.

The input files are the daily averaged atmosphere coupler history files:

    ARCHDIR/CASE/cpl/hist/CASE.cpl.hx.atm.24h.avrg.YYYY-MM-DD-00000.nc

for the years BEGYEAR to ENDYEAR. The steps are:

  1. Average the daily files of each month of each year (ncra)
  2. Average each calendar month over all the years (ncra), and put the 12 months in one
     file (ncrcat)
  3. Set the time axis to the middle of each month of YEAR (noleap calendar)

The output file is:

    OUTDIR/CASE/proc/CASE.cpl.hx.atm.monclim.BEGYEAR-ENDYEAR_cYYYYMMDD.nc

This needs the NCO commands ncra, ncrcat and ncatted (on Derecho: module load nco).
With --batch, the input is checked and then this tool is submitted to the PBS batch queue.

The general code is in cplhist_monthly_clim.py, this sets what's specific to aerosols for 1850.

A previous version of this for CESM2 is in NCL and here on Derecho by Keith Oleson:
    /glade/u/home/oleson/misc_programs/CLM5Dev/CMIP6/CreateAeroDepFile_Task5B.ncl
It did 5 year averages over the historical period and then did a linear interpolation
between the 5-year averages. It also used CPLHIST output from several ensemble members and
averged them together with ncea. This is something that will have to be done in the future
in this tool as well.

Another example for Nitrogen deposition is here from Simone Tilmes and Mike Mills:
   https://svn.code.sf.net/p/codescripts/code/trunk/ncl/cam_forcing/CreateDepositionFile.ncl

Neither were quite right for what's needed here. And both are in NCL which is well after
end of life now.

Also there is code in LDF to create monthly climatologies, but it's not setup to bring in ensemble
averages and not specialized to create datafiles for DATM aerosol forcing.
"""

import os

from ctsm.cplhist.cplhist_monthly_clim import CplhistMonthlyClim


class Aero1850MonthlyClim(CplhistMonthlyClim):
    """
    Monthly climatology of the aerosol deposition fields for 1850
    """

    DEFAULT_VARS = [
        "atmImp_Faxa_bcph1",
        "atmImp_Faxa_bcph2",
        "atmImp_Faxa_bcph3",
        "atmImp_Faxa_ocph1",
        "atmImp_Faxa_ocph2",
        "atmImp_Faxa_ocph3",
        "atmImp_Faxa_dstwet1",
        "atmImp_Faxa_dstdry1",
        "atmImp_Faxa_dstwet2",
        "atmImp_Faxa_dstdry2",
        "atmImp_Faxa_dstwet3",
        "atmImp_Faxa_dstdry3",
        "atmImp_Faxa_dstwet4",
        "atmImp_Faxa_dstdry4",
    ]
    DEFAULT_YEAR = 1850
    PBS_JOBNAME = "aero_1850_monthly_clim"
    WRAPPER = os.path.join("tools", "cplhist", "aero_1850_monthly_clim")
    TITLE = "Monthly Climatology of Daily Averaged Atmosphere Coupler History Files"
    COMMENT = "Monthly climatology created from daily averaged CPLHIST files for DATM"


def main():
    """
    Main program for aero_1850_monthly_clim
    """
    Aero1850MonthlyClim.main()
