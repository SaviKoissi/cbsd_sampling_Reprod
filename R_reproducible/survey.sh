#!/bin/bash
#SBATCH --job-name=cbsd_surv
#SBATCH --nodelist=node06
#SBATCH --output=R_reproducible/outputs/surv_%j.log
#SBATCH --error=R_reproducible/outputs/surv_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=05:00:00

# 1. Path Configuration
# Move into the folder where your script and simulation.rds are
cd /home/savi/project/SurveyCBSD/SurveyCBSD/R_reproducible

echo "JOB START: $(date)"
echo "Using Node: $SLURM_NODELIST"

# 2. Execution using Absolute Path
# This bypasses the need for 'module load'
/usr/local/bioinfo/R-4.5.2/bin/Rscript complet_section2.R

echo "JOB END: $(date)"