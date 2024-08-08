#!/bin/bash
set -e

module load nco

indir=/glade/campaign/cesm/cesmdata/cseg/inputdata/atm/datm7/anomaly_forcing
outdir=/glade/work/samrabin/ctsm_cdeps_ssp-temperature

rcp=$1
if [[ ${rcp} == "" ]]; then
    echo "You must provide rcp (e.g., 45)" >&2
    exit 1
fi

if [[ ! -d "${indir}" ]]; then
    echo "Input directory not found: ${indir}" >&2
    exit 1
fi

restoffilename="ccsm4.rcp${rcp}.2006-2300.nc"

pushd "${indir}" 1>/dev/null
varlist="$(ls af.*.${restoffilename} | grep -v allvars | cut -f2 -d.)"
if [[ "${varlist}" == "" ]]; then
    echo "No variables found" >&2
    exit 1
fi
popd 1>/dev/null

newfile="${outdir}/af.allvars.${restoffilename}"

i=0
for var in ${varlist}; do
    echo $var
    i=$((i + 1))
    thisfile="${indir}/af.${var}.${restoffilename}"
    if [[ ${i} -eq 1 ]]; then
        nccopy -k cdf5 "${thisfile}" "${newfile}"
        continue
    fi
    ncks -A -v ${var} "${thisfile}" "${newfile}"
done

echo Done

exit 0