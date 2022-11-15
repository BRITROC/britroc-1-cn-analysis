#!/bin/bash

set -e

#Rscript metadata_formatting/generate_meta_data.R
#Rscript -e "rmarkdown::render('generate_britroc_cohort_patient_data.Rmd')"
Rscript -e "rmarkdown::render('copy_number_signatures/britroc_30kb_ds_preanalysis_qc.Rmd')"

