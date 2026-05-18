#!/bin/bash
#SBATCH --job-name=cbsd_surv
#SBATCH --output=R_reproducible/outputs/surv_%j.log
#SBATCH --error=R_reproducible/outputs/surv_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=05:00:00

# 1. Initialize Environment
source /etc/profile

# 2. Load Base Modules
module purge
module load bioinfo-wave
module load R/4.5.2

# 3. ENVIRONMENT ISOLATION (Fixes GDAL / Driver Mismatch)
# Directs R to your isolated library space and overrides conflicting systemic links
export R_LIBS_USER="/home/savi/project/R_local_libs"
mkdir -p $R_LIBS_USER

OPENBLAS_DIR="/usr/local/bioinfo/miniconda3-26.1.1.1-1/pkgs/libopenblas-0.3.33-pthreads_h94d23a6_0/lib"
TBB_DIR="/usr/local/bioinfo/R-4.5.2/lib64/R/library/RcppParallel/lib"
R_LIB_DIR="/usr/local/bioinfo/R-4.5.2/lib64/R/lib"
export LD_LIBRARY_PATH=$OPENBLAS_DIR:$TBB_DIR:$R_LIB_DIR:/usr/local/lib64:/usr/lib64:$LD_LIBRARY_PATH

# 4. Navigation
cd /home/savi/project/SurveyCBSD/SurveyCBSD/R_reproducible

# 5. Logging and Execution
echo "===================================================="
echo "JOB START:    $(date)"
echo "JOB ID:       $SLURM_JOB_ID"
echo "RUNNING ON:   $SLURM_NODELIST"
echo "WORKDIR:      $(pwd)"
echo "===================================================="

echo "----------------------------------------------------"
echo "STEP 1: REGENERATING CLEAN RDS DATA MATRICES..."
echo "----------------------------------------------------"
# This script reads the base map via low-level binary streams, stripping out the NaN corruption
/usr/local/bioinfo/R-4.5.2/bin/Rscript generate_clean_data.R

echo "----------------------------------------------------"
echo "STEP 2: RUNNING BIO-ECONOMIC SIMULATION..."
echo "----------------------------------------------------"
# Ensure your complet_section2.R script reads the .rds files instead of .tif!
/usr/local/bioinfo/R-4.5.2/bin/Rscript complet_section2.R

echo "----------------------------------------------------"
echo "JOB END:      $(date)"
echo "===================================================="
