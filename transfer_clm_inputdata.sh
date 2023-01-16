#!/bin/bash -l
#SBATCH --job-name="transfer"
#SBATCH --time=02:00:00
#SBATCH --partition=xfer
#SBATCH --hint=nomultithread
#SBATCH --nodes=1
#SBATCH --job-name="store"
#SBATCH --output=transfer_clm_inputdata.out
#SBATCH --error=transfer_clm_inputdata.err

# Transfer CLM inputdata from the project folder to SCRATCH
if [ $PROJ == sm61 ]; then
 ORIGIN=/project/$PROJ/shared/CCLM2_inputdata/
elif [ $PROJ == sm62 ]; then
 ORIGIN=/project/$PROJ/CCLM2_inputdata/
else
    echo "Unknown project"
fi

TARGET=$SCRATCH/CCLM2_inputdata/

if [ ! -d $TARGET ]
then
  mkdir -p $TARGET
fi

rsync -av --progress ${ORIGIN} ${TARGET}

# Ensure that the production job can begin execution only after the stage job (transfer) has successfully executed
#srun rsync -av $1 $2
#if [ -n '$3' ]; then
#    sbatch --dependency=afterok:$SLURM_JOB_ID $3
#fi
