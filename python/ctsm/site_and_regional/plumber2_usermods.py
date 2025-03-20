#! /usr/bin/env python3

"""

Reads in .csv files with PLUMBER2 site information
Creates individual usermod_dirs for each PLUMBER2 site with shell_commands

"""

#  Import libraries
from __future__ import print_function

import os
import tqdm

import pandas as pd


# Big ugly function to create usermod_dirs for each site
def write_usermods(
    *,
    lat,
    lon,
    site,
    start_year,
    end_year,
    start_date,
    start_year_actual,
    start_tod,
    atm_ncpl,
    stop_n,
):
    """
    Write information to be added to user mods
    """

    site_dir = os.path.join("../../cime_config/usermods_dirs/clm/PLUMBER2/", site)

    if not os.path.isdir(site_dir):
        os.makedirs(site_dir, exist_ok=True)

    # create files in each directory
    include = os.path.join(site_dir, "include_user_mods")
    i_file = open(include, "w")  # or 'a' to add text instead of truncate
    i_file.write("../defaults")
    i_file.close()

    # pylint: disable=anomalous-backslash-in-string
    lai_stream = (
        "\$DIN_LOC_ROOT/lnd/clm2/lai_streams/PLUMBER2/"
        + site
        + "/LAI_stream_"
        + site
        + "_"
        + str(start_year)
        + "-"
        + str(end_year)
        + ".nc"
    )
    shell = os.path.join(site_dir, "shell_commands")
    s_file = open(shell, "w")  # or 'a' to add text instead of truncate
    # pylint: disable=line-too-long
    s_file.write(
        # TODO turn on following line after cdeps changes are added
        #'./xmlchange PLUMBER2SITE='+site + '\n' \
        "./xmlchange PTS_LON=" + str(lon) + "\n"
        "./xmlchange PTS_LAT=" + str(lat) + "\n"
        "./xmlchange DATM_YR_END=" + str(end_year) + "\n"
        "./xmlchange DATM_YR_START_FILENAME=" + str(start_year) + "\n"
        "./xmlchange START_TOD=" + str(start_tod) + "\n"
        "./xmlchange ATM_NCPL=" + str(atm_ncpl) + "\n"
        "\n"
        'echo "presaero.SSP3-7.0:year_first=' + str(start_year) + '" >> user_nl_datm_streams \n'
        'echo "presaero.SSP3-7.0:year_last=' + str(end_year) + '" >> user_nl_datm_streams \n'
        'echo "presaero.SSP3-7.0:year_align=' + str(start_year) + '" >> user_nl_datm_streams \n'
        "\n"
        'echo "presndep.SSP3-7.0:year_first=' + str(start_year) + '" >> user_nl_datm_streams \n'
        'echo "presndep.SSP3-7.0:year_last=' + str(end_year) + '" >> user_nl_datm_streams \n'
        'echo "presndep.SSP3-7.0:year_align=' + str(start_year) + '" >> user_nl_datm_streams \n'
        "\n"
        'echo "co2tseries.SSP3-7.0:year_first=' + str(start_year) + '" >> user_nl_datm_streams \n'
        'echo "co2tseries.SSP3-7.0:year_last=' + str(end_year) + '" >> user_nl_datm_streams \n'
        'echo "co2tseries.SSP3-7.0:year_align=' + str(start_year) + '" >> user_nl_datm_streams \n'
        "\n"
        "compset=`./xmlquery COMPSET --value` \n"
        "CLM_USRDAT_NAME=`./xmlquery CLM_USRDAT_NAME --value` \n"
        "TEST=`./xmlquery TEST --value` \n"
        "\n"
        "# For a transient case run the whole length and do not cycle \n"
        "if  [[ $compset =~ ^HIST ]]; then \n"
        "  # Number of years that can be run for the full transient case \n"
        '  if [[ $TEST != "TRUE" ]]; then  \n'
        "    ./xmlchange STOP_N=" + str(stop_n) + "\n"
        "  fi \n"
        "  # set start date for transient case with historical compset \n"
        "  ./xmlchange RUN_STARTDATE=" + str(start_date) + "\n"
        "  ./xmlchange DATM_YR_ALIGN=" + str(start_year_actual) + "\n"
        "  ./xmlchange DATM_YR_START=" + str(start_year_actual) + "\n"
        "else \n"
        "  # for spinup case with I2000 compset \n"
        "  ./xmlchange RUN_STARTDATE=0001-01-01" + "\n"
        "  ./xmlchange DATM_YR_ALIGN=" + str(1) + "\n"
        "  ./xmlchange DATM_YR_START=" + str(start_year) + "\n"
        "fi \n"
        "\n"
        "# Turn on LAI streams for a SP case \n"
        "if [[ $compset =~ .*CLM[0-9]+%[^_]*SP.* ]]; then \n"
        "  echo \"stream_fldfilename_lai='" + lai_stream + "'\" >> user_nl_clm \n"
        '  echo "stream_year_last_lai=' + str(end_year) + '" >> user_nl_clm \n'
        "  if  [[ $compset =~ ^HIST ]]; then \n"
        "    # for transient case with a historical compset \n"
        '    echo "model_year_align_lai=' + str(start_year_actual) + '" >> user_nl_clm \n'
        '    echo "stream_year_first_lai=' + str(start_year_actual) + '" >> user_nl_clm \n'
        "  else \n"
        "    # for a spinup case with a i2000 compset \n"
        '    echo "model_year_align_lai=1" >> user_nl_clm \n'
        '    echo "stream_year_first_lai=' + str(start_year) + '" >> user_nl_clm \n'
        "  fi \n"
        "fi \n"
        "\n"
    )
    # pylint: enable=line-too-long, anomalous-backslash-in-string

    s_file.close()

    # add baseflow_scalar = 0 to user_nl_clm for wetland sites
    wetland = [
        "CZ-wet",
        "DE-SfN",
        "FI-Kaa",
        "FI-Lom",
        "RU-Che",
        "SE-Deg",
        "US-Los",
        "US-Myb",
        "US-Tw4",
        "PL-wet",
    ]
    if any(x == site for x in wetland):
        s_file = open(shell, "a")  # or 'a' to add text instead of truncate
        s_file.write(
            "\n"
            "# set baseflow scalar to zero for wetland site \n"
            'echo "baseflow_scalar = 0" >> user_nl_clm'
        )
        s_file.close()


# End write_usermods function


def main():
    """
    Iterate through plumber2 sites and create usermod_dirs
    """

    # For now we can just run the 'main' program as a loop
    plumber2_sites = pd.read_csv("PLUMBER2_sites.csv", skiprows=4)

    for _, row in tqdm.tqdm(plumber2_sites.iterrows()):
        lat = row["Lat"]
        lon = row["Lon"]
        site = row["Site"]
        start_year = row["start_year"]
        end_year = row["end_year"]
        start_date = row["RUN_STARTDATE"]
        start_year_actual = start_date[:4]
        start_tod = row["START_TOD"]
        atm_ncpl = row["ATM_NCPL"]
        stop_n = 1 + end_year - start_year

        write_usermods(
            lat=lat,
            lon=lon,
            site=site,
            start_year=start_year,
            end_year=end_year,
            start_date=start_date,
            start_year_actual=start_year_actual,
            start_tod=start_tod,
            atm_ncpl=atm_ncpl,
            stop_n=stop_n,
        )


if __name__ == "__main__":
    main()
