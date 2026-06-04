#!/bin/bash

#SBATCH --time 18:00:00
#SBATCH --ntasks=1      # Processes per task
#SBATCH --cpus-per-task=2 # Allocate num CPUs (threads) to each task
#SBATCH --mem-per-cpu=4G
#SBATCH -J "run_snakefile"   # job name
#SBATCH -o logs/snake_%A.log

REPO_DIR="$(cd "$(dirname "$0")" && pwd)"

# Load OUTPUT_FOLDER from .env if present
if [ -f "$REPO_DIR/.env" ]; then
    set -a
    source "$REPO_DIR/.env"
    set +a
fi

# Cleanup old logs and temporary files to prevent space issues
echo "Cleaning up old logs and temporary files..."
find "$HOME/.snakemake" -type f -mtime +7 -delete 2>/dev/null || true

# Setup output directory with group permissions
if [ -n "$OUTPUT_FOLDER" ]; then
    echo "Setting up output directory: $OUTPUT_FOLDER"
    mkdir -p "$OUTPUT_FOLDER"
    find -L "$OUTPUT_FOLDER" -type d -exec chmod g+s {} +
fi

# Load system profile to ensure sbatch and module commands are available
if [ -f /etc/profile ]; then
    source /etc/profile
fi

# The slurm account and slurm partition are essential for grouping
# Choosing the solver helps snakemake find it
echo "Starting Snakemake"
pixi run --no-progress snakemake -s "$REPO_DIR/Snakefile" \
    --configfile "$REPO_DIR/config.yaml" \
    --scheduler-ilp-solver COIN_CMD \
    --executor slurm \
    --default-resources slurm_account=srp33 \
    --jobs 600 \
    --group-components batch_real_group=2 batch_simulation_group=3 class_imbalance_group=2 \
    --max-jobs-per-second 20 \
    --max-status-checks-per-second 10 \
    --latency-wait 40 \
    --rerun-incomplete
    # --keep-going
    # --quiet
