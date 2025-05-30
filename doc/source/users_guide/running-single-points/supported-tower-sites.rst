.. include:: ../substitutions.rst

.. _supported-tower-sites:

********************************************
Supported tower sites for single-point runs
********************************************

CTSM has functionality within the ``run_tower`` tool for running single-point cases at particular supported tower sites using forcing data from those sites.

====================================================
General Information on Running Supported Tower Sites
====================================================

The ``run_tower`` capability allows users to run Community Land Model (CLM) simulations at NEON and PLUMBER tower sites in a streamlined manner by setting up the appropriate model configurations, datasets, and initial conditions. This script can run for one or more (NEON or PLUMBER) tower sites. It will do the following
    1) Create a generic base case for cloning.
    2) Make the case for the specific neon or plumber site(s).
    3) Make changes to the case, for
        a. AD spinup
        b. post-AD spinup
        c. transient
        d. SASU or Matrix spinup
    4) Build and submit the case.

The available options can be shown by running ``run_tower --help``

These options include the following:
  -h, --help            show this help message and exit
  --neon-sites          4-letter neon site code.
  --plumber-sites       six character PLUMBER2 site code (eg, AR-SLu)
  --base-case BASE_CASE_ROOT
                        Root Directory of base case build [default: None]
  --output-root OUTPUT_ROOT
                        Root output directory of cases [default:
                        CIME_OUTPUT_ROOT as defined in cime]
  --overwrite           overwrite existing case directories [default: False]
  --setup-only          Only setup the requested cases, do not build or run
                        [default: False]
  --no-input-data-check, --no-check-input-data
                        Don't check for input data. Implies --setup-only.
                        [default: False]
  --rerun               If the case exists but does not appear to be complete,
                        restart it. [default: False]
  --no-batch            Run locally, do not use batch queueing system (if
                        defined for Machine) [default: False]
  --run-type            ad, postad, transient
                        Type of run to do [default: None]
  --prism               Uses the PRISM reanaylsis precipitation data for the
                        site instead of the NEON data (only available over
                        Continental US)
  --experiment EXPERIMENT
                        Appends the case name with string for model experiment
  --run-from-postad     For transient runs only. By default start from
                        finidat, if this flag is used the postad run must be
                        available.
  --neon-version        v1,v2,v3
                        Neon data version to use for this simulation.
                        [default: use the latest data available]
  --xmlchange XMLCHANGE
                        Any xmlchanges (e.g.,
                        CLM_CO2_TYPE=constant,CCSM_CO2_PPMV=500) [default:
                        None]

Logging options:
  -d, --debug           Print debug information (very verbose) to file /glade/
                        work/tking/ctsm_project/sam_ctsm/CTSM/tools/site_and_r
                        egional/run_tower.log
  -v, --verbose         Add additional context (time and file) to log messages
  -s, --silent          Print only warnings and error messages


A `tutorial <https://ncar.github.io/CTSM-Tutorial/README.html>`_ on running ``run_tower`` is also available.

.. warning::
Note that the run_tower base case must be of same run type as a requested clone, as described by this `issue ticket <https://github.com/ESCOMP/CTSM/issues/1926>`_.

=========================================
NEON Tower Single Point Simulations
=========================================

With this tool, CLM uses gap-filled meteorology from NEON tower sites, the dominant plant species is mapped to the appropriate model plant functional type (PFT), and soil characteristics used in the simulations are updated to match observations from NEON’s soil megapits. Gap-filled NEON tower flux data are also available for model evaluation. Additionally, all the commands to run the model are combined into a script that you can easily call from a single line of code.

Currently supported NEON sites include the following: ABBY, BARR, BART, BLAN, BONA, CLBJ, CPER, DCFS, DEJU, DELA, DSNY, GRSM, GUAN, HARV, HEAL, JERC, JORN, KONA, KONZ, LAJA, LENO, MLBS, MOAB, NIWO, NOGP, OAES, ONAQ, ORNL, OSBS, PUUM, RMNP, SCBI, SERC, SJER, SOAP, SRER, STEI, STER, TALL, TEAK, TOOL, TREE, UKFS, UNDE, WOOD, WREF, YELL, all
Information on the specific sites can be found on the `NEON webpage <https://www.neonscience.org/field-sites>`_.

.. note:: An important note regarding the NEON tower site simulations is that the default run type is `transient`.


=========================================
PLUMBER Tower Single Point Simulations
=========================================

.. note:: A few important notes regarding the PLUMBER tower site simulations are that the default run type is `ad`. Additionally, PLUMBER cases all start in different years.

Currently supported PLUMBER Sites include the following: AR-SLu, AT-Neu, AU-ASM, AU-Cow, AU-Cpr, AU-Ctr, AU-Cum, AU-DaP, AU-DaS, AU-Dry, AU-Emr, AU-GWW, AU-Gin, AU-How, AU-Lit, AU-Otw, AU-Rig, AU-Rob, AU-Sam, AU-Stp, AU-TTE, AU-Tum, AU-Whr, AU-Wrr, AU-Ync, BE-Bra, BE-Lon, BE-Vie, BR-Sa3, BW-Ma1, CA-NS1, CA-NS2, CA-NS4, CA-NS5, CA-NS6, CA-NS7, CA-Qcu, CA-Qfo, CA-SF1, CA-SF2, CA-SF3, CH-Cha, CH-Dav, CH-Fru, CH-Oe1, CN-Cha, CN-Cng, CN-Dan, CN-Din, CN-Du2, CN-HaM, CN-Qia, CZ-wet, DE-Bay, DE-Geb, DE-Gri, DE-Hai, DE-Kli, DE-Meh, DE-Obe, DE-Seh, DE-SfN, DE-Tha, DE-Wet, DK-Fou, DK-Lva, DK-Ris, DK-Sor, DK-ZaH, ES-ES1, ES-ES2, ES-LMa, ES-LgS, ES-VDA, FI-Hyy, FI-Kaa, FI-Lom, FI-Sod, FR-Fon, FR-Gri, FR-Hes, FR-LBr, FR-Lq1, FR-Lq2, FR-Pue, GF-Guy, HU-Bug, ID-Pag, IE-Ca1, IE-Dri, IT-Amp, IT-BCi, IT-CA1, IT-CA2, IT-CA3, IT-Col, IT-Cpz, IT-Isp, IT-LMa, IT-Lav, IT-MBo, IT-Mal, IT-Noe, IT-Non, IT-PT1, IT-Ren, IT-Ro1, IT-Ro2, IT-SR2, IT-SRo, JP-SMF, NL-Ca1, NL-Hor, NL-Loo, PL-wet, PT-Esp, PT-Mi1, PT-Mi2, RU-Che, RU-Fyo, RU-Zot, SD-Dem, SE-Deg, UK-Gri, UK-Ham, UK-PL3, US-AR1, US-AR2, US-ARM, US-Aud, US-Bar, US-Bkg, US-Blo, US-Bo1, US-Cop, US-FPe, US-GLE, US-Goo, US-Ha1, US-Ho1, US-KS2, US-Los, US-MMS, US-MOz, US-Me2, US-Me4, US-Me6, US-Myb, US-NR1, US-Ne1, US-Ne2, US-Ne3, US-PFa, US-Prr, US-SP1, US-SP2, US-SP3, US-SRG, US-SRM, US-Syv, US-Ton, US-Tw4, US-Twt, US-UMB, US-Var, US-WCr, US-Whs, US-Wkg, ZA-Kru, ZM-Mon, all
Information on the se=pecific sites can be found `here <https://researchdata.edu.au/plumber2-forcing-evaluation-surface-models/1656048>`_.