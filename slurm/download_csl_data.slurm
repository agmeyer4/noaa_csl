#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=lin-kp
#SBATCH --partition=lin-kp
#SBATCH -o ./logs/csl_download_%A.out

source activate atmos_column1
cd /uufs/chpc.utah.edu/common/home/u0890904/NOAA_CSL/noaa_csl
python -u data_download.py