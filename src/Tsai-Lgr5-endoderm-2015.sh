
#! /bin/bash
# First execute tuxedo pipeline to generate FPKM data (this will take a long time)
./tuxedo_complete.sh
# Then execute the R script
Rscript expression_analysis.R
