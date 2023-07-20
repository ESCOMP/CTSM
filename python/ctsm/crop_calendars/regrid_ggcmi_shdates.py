from subprocess import run
import os
import shutil
import glob

def run_and_check(cmd):
    result = run(
        cmd,
        shell = True,
        capture_output = True,
        text = True,
        )
    if result.returncode != 0:
        raise RuntimeError(f"Trouble running `{result.args}` in shell:\n{result.stderr}")

def main():

    # Ensure we can call necessary shell scripts
    for cmd in ["ncks", "ncrename", "ncpdq", "cdo"]:
        run_and_check(f"{cmd} --help")
        
    templatefile0 = "/Users/Shared/CESM_work/CropEvalData_ssr/danica_timeseries-cmip6_i.e21.IHIST.f09_g17/month_1/ssr_trimmed_annual.nc"

    templatefile = "template.nc"
    
    # For some reason, doing ncks -v directly doesn't work. Have to copy it over first.
    shutil.copyfile(templatefile0, os.path.basename(templatefile0))
    templatefile0 = os.path.basename(templatefile0)

    if os.path.exists(templatefile):
        os.remove(templatefile)

    for v in ["planting_day", "maturity_day", "growing_season_length"]:
        run_and_check(f"ncks -A -v area '{templatefile0}' '{templatefile}'")
        run_and_check(f"ncrename -v area,{v} '{templatefile}'")
    os.remove(templatefile0)
    run_and_check(f"ncpdq -O -h -a -lat '{templatefile}' '{templatefile}'")


    input_files = glob.glob("/Users/Shared/GGCMI/AgMIP.input/phase3/ISIMIP3/crop_calendar/*nc4")
    input_files.sort()
    for f in input_files:

        f2 = os.path.basename(f)
        f3 = f2.replace(".nc4", "_nninterp-f09_g17.nc4")
        print(f3)

        if os.path.exists(f3):
            os.remove(f3)
        run_and_check(f"cdo -L -remapnn,'{templatefile}' -setmisstonn '{f}' '{f3}'")

    os.remove(templatefile)

if __name__ == "__main__":
    main()