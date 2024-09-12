#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=carbon-kp
#SBATCH --partition=carbon-kp
#SBATCH -o ./logs/csl_regrid_%A.out

source activate noaa_csl
cd /uufs/chpc.utah.edu/common/home/u0890904/NOAA_CSL/noaa_csl

python -u regrid_data.py

# Extract the path from the 4th line of the log file
details_path=$(sed -n '4s/^Saving regridded .nc files to //p' ./slurm/logs/csl_regrid_${SLURM_JOB_ID}.out)/details

# Copy the log file to the path returned by the Python script
echo "Copying log file to $details_path"
cp ./slurm/logs/csl_regrid_${SLURM_JOB_ID}.out "$details_path"
