#!/bin/bash
#SBATCH --job-name=MLDD-sac
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cbowers@stanford.edu
#SBATCH --partition=cee
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=11
#SBATCH --mem=240G
#SBATCH --output=_scripts/fitrf_sacramento_1981_haz.log

## load modules
ml R/4.1.2
ml physics
ml proj gdal geos udunits

## execute script(s)
cd /scratch/users/cbowers/6-MLDD
Rscript _scripts/fitrf_sacramento_1981_haz.R 1981_haz