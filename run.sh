#!/bin/bash

# For using on local machine just start script directly
# Rscript ./calc_enrichment.R

## call the script
singularity run --app Rscript /projectnb2/wax-es/routines/singularity/r442.sif ./calc_enrichment.R
