#!/bin/bash
#SBATCH --job-name=MLDD-sac
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cbowers@stanford.edu
#SBATCH --partition=cee
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=240G
#SBATCH --output=_scripts/MLDD_results_sacramento.log

## load modules
ml R/4.1.2
ml physics
ml proj gdal geos udunits

## execute script(s)
cd /scratch/users/cbowers/6-MLDD
Rscript _scripts/MLDD_results_sacramento.R