#!/bin/bash

set -e

echo -e "Generate metadata"
#Rscript metadata_formatting/generate_meta_data.R
echo -e "Generate patient data"
#Rscript -e "rmarkdown::render('generate_britroc_cohort_patient_data.Rmd')"
echo -e "Pre-signature QC"
#Rscript -e "rmarkdown::render('copy_number_signatures/britroc_30kb_ds_preanalysis_qc.Rmd')"
echo -e "generate signatures"
#Rscript -e "rmarkdown::render('copy_number_signatures/britroc_30kb_cnsig_gen.Rmd')"
echo -e "OV signature analysis"
#Rscript -e "rmarkdown::render('copy_number_signatures/britroc_30kb_copy_number_signature_analysis.Rmd')"
echo -e "Focal copy number analysis"
#Rscript -e "rmarkdown::render('copy_number_analysis/focal_analysis/britroc_30kb_abs_cna_analysis_focal.Rmd')"
echo -e "Broad copy number analysis"
#Rscript -e "rmarkdown::render('copy_number_analysis/broad_analysis/britroc_30kb_abs_cna_analysis_broad.Rmd')"
echo -e "Ploidy change QC"
Rscript -e "rmarkdown::render('copy_number_analysis/broad_analysis/ploidy_change_patients_scoring.Rmd')"

