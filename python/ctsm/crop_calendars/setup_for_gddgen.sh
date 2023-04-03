#!/bin/bash
set -e

# This script provides a general method to set up a GDD-generating run from
# an existing case directory. Call it from the case directory, preferably
# on a compute node, as the "If needed, generate a surface dataset file"
# step can be intensive.

cropcals_python_dir="/glade/u/home/samrabin/CTSM_cropcals_hist/crop_calendars"

# Set up
cp user_nl_clm.orig user_nl_clm
./check_case 1>/dev/null

# Ensure that run will go for at least 4 years; otherwise error.
stop_option=$(./xmlquery --value STOP_OPTION)
stop_n=$(./xmlquery --value STOP_N)
min_Nyears=4
errMsg="GDD-generating runs must be at least ${min_Nyears} years long."
if [[ "${stop_option}" == "nyear"* ]]; then
    if [[ "${stop_n}" -lt ${min_Nyears} ]]; then
        echo ${errMsg} >&2
        exit 1
    fi
elif [[ "${stop_option}" == "nmonth"* ]]; then
    if [[ "${stop_n}" -lt $((min_Nyears * 12)) ]]; then
        echo ${errMsg} >&2
        exit 1
    fi
elif [[ "${stop_option}" == "nday"* ]]; then
    if [[ "${stop_n}" -lt $((min_Nyears * 365)) ]]; then
        echo ${errMsg} >&2
        exit 1
    fi
else
    echo "STOP_OPTION '${stop_option}' not recognized by preprocess_inputs.sh" >&2
    exit 1
fi

# If needed, generate a surface dataset file with no crops missing years
flanduse_timeseries="$(grep -h -i "flanduse_timeseries" CaseDocs/* | sed "s/ flanduse_timeseries = //" | sed "s/'//g")"
if [[ "${flanduse_timeseries}" != "" ]]; then
    if [[ ! -e "${flanduse_timeseries}" ]]; then
        echo "flanduse_timeseries not found: ${flanduse_timeseries}" >&2
        exit 1
    fi

    # Make new fsurdat file (in case directory)
    module unload python
    module load conda
    conda activate npl
    fsurdat="$(grep -h -i "fsurdat" CaseDocs/* | sed "s/ fsurdat = //" | sed "s/'//g")"
    paramfile="$(grep -h -i "paramfile" CaseDocs/* | sed "s/ paramfile = //" | sed "s/'//g")"
    set +e
    python "${cropcals_python_dir}/make_surface_for_gddgen.py" --flanduse_timeseries "${flanduse_timeseries}" --fsurdat "${fsurdat}" --paramfile "${paramfile}" > make_surface_for_gddgen.log
    if [[ $? -ne 0 ]]; then
        cat make_surface_for_gddgen.log
        exit 1
    fi
    set -e
    new_fsurdat="$PWD/$(cat make_surface_for_gddgen.log)"
    rm make_surface_for_gddgen.log

    # Change fsurdat to point to that file and disable transient crops
    echo "fsurdat = '${new_fsurdat}'" >> user_nl_clm
    echo "do_transient_crops = .false." >> user_nl_clm
    echo "flanduse_timeseries = ''" >> user_nl_clm
fi

# user_nl_clm: Replace MESHFILE_PLACEHOLDER with the appropriate file for this resolution
file_mesh="$(./xmlquery --value LND_DOMAIN_MESH)"
sed -i "s@MESHFILE_PLACEHOLDER@${file_mesh}@" user_nl_clm

# user_nl_clm: Replace SDATEFILE_PLACEHOLDER with the appropriate file for this resolution.
# EVENTUALLY: Generate sdate (and hdate) file at case resolution from GGCMI files; save to case directory and use that.
lnd_grid=$(./xmlquery --value LND_GRID)
blessed_crop_dates_dir="/glade/work/samrabin/crop_dates_blessed"
if [[ "${lnd_grid}" == "10x15" ]]; then
    file_sdates="${blessed_crop_dates_dir}/sdates_ggcmi_crop_calendar_phase3_v1.01_nninterp-f10_f10_mg37.2000-2000.20230330_165301.fill1.nc"
elif [[ "${lnd_grid}" == "1.9x2.5" ]]; then
    file_sdates="${blessed_crop_dates_dir}/sdates_ggcmi_crop_calendar_phase3_v1.01_nninterp-f19_g17.2000-2000.20230102_175625.fill1.nc"
else
    echo "Land resolution (LND_GRID) '${lnd_grid}' not recognized by preprocess_inputs.sh" >&2
    exit 1
fi
sed -i "s@SDATEFILE_PLACEHOLDER@${file_sdates}@" user_nl_clm

# Rebuild namelists
./check_case 1>/dev/null

exit 0
