#!/bin/bash
set -e

cropcals_python_dir="/glade/u/home/samrabin/CTSM_cropcals_hist/crop_calendars"

# Where did the run happen? If not provided, the assumption is
# that this script is being called as POSTRUN_SCRIPT. Therefore,
# the outputs will still be in the run directory. We also assume
# this is being run as part of a test suite, so we'll save the
# generate_gdds.py outputs to the case directory (so subsequent
# tests can find them easily).
run_outputs="$1"
if [[ "${run_outputs}" == "" ]]; then
    run_outputs="$(./xmlquery --value RUNDIR)"
    save_to="."
else
    save_to="${run_outputs}"
fi

# Set ${y1} to third year in run
y1=$(ls -1 *clm2.h1* | grep -oE "h1.[0-9]{4}" | sed "s/h1.//" | head -n 1)
y1=$((y1 + 2))

# Set ${yN} to last full year – 1
yN=$(ls -1 *clm2.h1* | grep -oE "h1.[0-9]{4}" | sed "s/h1.//" | tail -n 1)
yN=$((yN - 1))

# Get ${sdates_file} from namelists
sdates_file="$(realpath "$(grep "sdate" CaseDocs/* | grep -oE "'.*'" | sed "s/'//g" | sed "s/.fill1//")")"
# Get ${hdates_file} from ${sdates_file}
hdates_file="${sdates_file/sdates/hdates}"

#source ~/.bash_profile

set +e
module is-loaded conda
conda_not_loaded=$?
set -e
if [[ ${conda_not_loaded} -ne 0 ]]; then
    module unload python
    module load conda
fi
conda activate npl
date
echo "Starting..."
python "${cropcals_python_dir}"/generate_gdds.py -r "${run_outputs}" -sd "${sdates_file}" -hd "${hdates_file}" -1 ${y1} -n ${yN} --dont-save-figs -o "${save_to}"

exit 0
