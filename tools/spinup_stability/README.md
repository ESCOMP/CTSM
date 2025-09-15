# Spinup stability
## a python script to assess CLM spinup stability

## Requirements:

1. ctsm_pylib
    - see zqz for details

2. a completed CLM simulation 
    - needs area and landfrac variables
    - should have completed at least two forcing cycles

3. spinup stability configuration file
    - more info below        

4. missing bits
    - only implemented for annualized output presently
    - only implemented for regular grids
        - no sparsegrid
        - no spectral element

## Files

1. spinup_stability.py
    - a python script to calculate and evaluate model drift
    - will create a set of diagnostic plots
    - exits with 0 if the simulation meets all criteria
    - exits with 11 if the simulation is drifting

2. default.yml
    - a configuration file for spinup_stability.py
    - includes our standard stability criteria
    - points to the case that will be evaulated
    - thresholds with _gridded, are interpreted as pct landarea gridded disequilibrium
    - threshoulds without _gridded, are interpreted as global drifts
    - cf = conversion factor
        - land area (la) multiplication is also already included
        - i.e. should be interpreted as GPP=cf*(la*ds.GPP).sum(dim=['lat','lon'])
    - lists of conversion factors will be multiplied together (see for example GPP)
    
3. batch.job
    - an example job submission script in lieu of running from the command line

## Usage

1. load a suitable python environment.
2. python spinup_stability.py my_config.yml

May be preferred as a batch job, e.g.
qsub batch.job

Output:
 - The script will write some text to stdout detailing spinup status.
 - The script will exit 0 for spun up, 11 for drifting, 1 for other problems.
 - The script will also create a set of plots $CASE.png
            