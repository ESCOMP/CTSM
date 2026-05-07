#!/bin/bash
# Debug helper: submit a PBS job that runs ./driver directly (no
# verify.sh grep filter) so any GPU runtime error is fully visible.
#
# Assumes ./driver was already built on the login node (e.g. via:
#   make clean && make EXTRA_FFLAGS="-acc=gpu -gpu=cc80"
# ).
#
# Usage:
#   ./debug_gpu.sh                  # runs ./driver --fast
#   ./debug_gpu.sh --all            # runs ./driver --all
#   ./debug_gpu.sh --fast           # explicit
#
# Output is written to ./sbgc_dbg.o<jobid> (gitignored) and cat'd here.

set -euo pipefail
cd "$(dirname "$0")"

driver_args="${*:---fast}"

job_id=$(qsub -W block=true \
              -A ucsg0003 -q develop \
              -l select=1:ncpus=1:ngpus=1 \
              -l walltime=00:05:00 \
              -N sbgc_dbg \
              -j oe \
              <<EOF
#!/bin/bash
cd "\$PBS_O_WORKDIR"
. ../env.sh
echo "=== nvidia-smi ==="
nvidia-smi || echo "(nvidia-smi failed)"
echo "=== ./driver $driver_args ==="
./driver $driver_args
echo "=== exit status: \$? ==="
EOF
)

job_num=${job_id%%.*}
out_file="sbgc_dbg.o${job_num}"

echo "=== job: $job_id (output: $out_file) ==="
if [ -f "$out_file" ]; then
    cat "$out_file"
else
    echo "(output file $out_file not found)"
    exit 1
fi
