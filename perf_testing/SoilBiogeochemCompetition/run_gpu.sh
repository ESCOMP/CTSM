#!/bin/bash
# Submit a non-interactive PBS job to a GPU node, run verify.sh there
# (with optional flags passed through), wait for completion, and print
# the job's stdout/stderr.
#
# Defaults: account ucsg0003, queue 'tutorial', 1 GPU, 5 min walltime.
# Override walltime via env var (e.g. WALLTIME=00:30:00 ./run_gpu.sh).
#
# Examples:
#   ./run_gpu.sh                                        # default verify on GPU node
#   ./run_gpu.sh EXTRA_FFLAGS="-acc=gpu -gpu=cc80"      # GPU build flags
#   ./run_gpu.sh INNER_TIMING=1 EXTRA_FFLAGS="-acc=gpu -gpu=cc80"
#   WALLTIME=00:30:00 ./run_gpu.sh EXTRA_FFLAGS="-acc=gpu -gpu=cc80"
#
# All script args are passed through verbatim to verify.sh inside the job.
# Job stdout+stderr are joined and written to ./sbgc_gpu.o<jobid>
# (gitignored). The script cats the output file before exiting.

set -euo pipefail
cd "$(dirname "$0")"

walltime="${WALLTIME:-00:05:00}"

# Re-quote each arg with printf %q so values containing spaces (e.g.
# EXTRA_FFLAGS="-acc=gpu -gpu=cc80") survive the heredoc round-trip.
quoted_args=""
for arg in "$@"; do
    quoted_args+=" $(printf '%q' "$arg")"
done

# qsub -W block=true (PBS Pro) submits and waits for the job to finish,
# returning the jobid on stdout and the job's exit status as its own.
job_id=$(qsub -W block=true \
              -A ucsg0003 -q tutorial \
              -l select=1:ncpus=1:ngpus=1 \
              -l walltime="$walltime" \
              -N sbgc_gpu \
              -j oe \
              <<EOF
#!/bin/bash
cd "\$PBS_O_WORKDIR"
. ../env.sh
./verify.sh$quoted_args
EOF
)

job_num=${job_id%%.*}
out_file="sbgc_gpu.o${job_num}"

echo "=== job: $job_id (output: $out_file) ==="
if [ -f "$out_file" ]; then
    cat "$out_file"
else
    echo "(output file $out_file not found — check qstat / job state)"
    exit 1
fi
