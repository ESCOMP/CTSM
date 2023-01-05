#!/bin/bash -l
#SBATCH --job-name="transfer"
#SBATCH --account="sm62"
#SBATCH --time=06:00:00
#SBATCH --partition=xfer
#SBATCH --hint=nomultithread
#SBATCH --nodes=1
#SBATCH --job-name="store"
#SBATCH --output=transfer_clm_inputdata.out
#SBATCH --error=transfer_clm_inputdata.err

# Transfer CLM inputdata back and forth from the project folder


ORIGIN=/project/sm62/CCLM2_inputdata/
TARGET=/scratch/snx3000/$USER/CCLM2_inputdata/


if [ ! -d $TARGET ]
then
  mkdir -p $TARGET
fi

rsync -av --progress ${ORIGIN} ${TARGET}

#srun rsync -av $1 $2
#if [ -n '$3' ]; then
#    sbatch --dependency=afterok:$SLURM_JOB_ID $3
#fi
