#!/bin/bash
set -e

ssp=$1
if [[ ${ssp} == "" ]]; then
    echo "You must provide an SSP number (e.g., 245)" >&2
    exit 1
fi
ssp_punctuation="$(echo ${ssp} | sed -E "s/(.)(.)(.)/\1-\2.\3/")"

casedir="$HOME/cases_ctsm/test_cdeps_ssp-temperature"
compset="SSP${ssp}_DATM%GSWP3v1_CLM51%BGC-CROP_SICE_SOCN_MOSART_CISM2%NOEVOLVE_SWAV"
res="f10_f10_mg37"

if [[ -d "${casedir}" ]]; then
    rm -r "${casedir}"
fi

cime/scripts/create_newcase --run-unsupported --case "${casedir}" --compset ${compset} --res ${res} \
    --handle-preexisting-dirs r

cd "${casedir}"
./case.setup

#echo "anomaly_forcing = 'Anomaly.Forcing.Temperature'" >> user_nl_datm
echo "anomaly_forcing = 'Anomaly.Forcing.cmip5.rcp45'" >> user_nl_datm
echo "flanduse_timeseries = '/glade/campaign/cesm/cesmdata/cseg/inputdata/lnd/clm2/surfdata_esmf/ctsm5.2.0/landuse.timeseries_10x15_SSP${ssp_punctuation}_78_CMIP6_1850-2100_c230517.nc'" >> user_nl_clm

echo "Generating namelists..."
./preview_namelists  1>/dev/null

./xmlquery COMPSET
grep "<file>" CaseDocs/datm.streams.xml | grep anomaly_forcing | grep '.tas.'

exit 0
