#!/bin/bash

set -e

echo -e "Generate metadata"
Rscript metadata_formatting/generate_meta_data.R
echo -e "Generate patient data"
Rscript -e "rmarkdown::render('generate_britroc_cohort_patient_data.Rmd')"
echo -e "Pre-signature QC"
Rscript -e "rmarkdown::render('copy_number_signatures/britroc_30kb_ds_preanalysis_qc.Rmd')"
echo -e "generate signatures"
Rscript -e "rmarkdown::render('copy_number_signatures/britroc_30kb_cnsig_gen.Rmd')"
echo -e "OV signature analysis"
Rscript -e "rmarkdown::render('copy_number_signatures/britroc_30kb_copy_number_signature_analysis.Rmd')"
echo -e "Focal copy number analysis"
Rscript -e "rmarkdown::render('copy_number_analysis/focal_analysis/britroc_30kb_abs_cna_analysis_focal.Rmd')"
echo -e "Broad copy number analysis"
Rscript -e "rmarkdown::render('copy_number_analysis/broad_analysis/britroc_30kb_abs_cna_analysis_broad.Rmd')"
echo -e "Ploidy change QC"
Rscript -e "rmarkdown::render('copy_number_analysis/broad_analysis/ploidy_change_patients_scoring.Rmd')"
echo -e "Heterogeneity analysis"
Rscript -e "rmarkdown::render('copy_number_analysis/heterogeneity/britroc_30kb_heterogeneity_analysis.Rmd')"
echo -e "sites of relapse analysis"
Rscript -e "rmarkdown::render('copy_number_analysis/sites_of_relapse/britroc_30kb_sites_of_relapse.Rmd')"
echo -e "Focal copy number analysis - patient specific"
Rscript -e "rmarkdown::render('copy_number_analysis/focal_analysis/britroc_30kb_abs_cna_analysis_focal_pat_specific.Rmd')"
echo -e "Focal copy number analysis - pt treatment"
Rscript -e "rmarkdown::render('copy_number_analysis/focal_analysis/britroc_30kb_abs_cna_analysis_primary_resistance.Rmd')"
echo -e "Focal copy number analysis - BRCA-HRD analysis"
Rscript -e "rmarkdown::render('copy_number_analysis/focal_analysis/britroc_30kb_cna_analysis_brca_and_hrd_stratification.Rmd')"
echo -e "Pan-cancer signature analysis"
Rscript -e "rmarkdown::render('copy_number_signatures/britroc_30kb_pan_cancer_copy_number_signature_analysis.Rmd')"
echo -e "Signature-immune analysis"
Rscript -e "rmarkdown::render('copy_number_signatures/britroc_30kb_cn_signatures_immune_analysis.Rmd')"
echo -e "Pipeline end"
