#!/bin/bash

script="install_env"
ENV="britroc"
# VARS
#CONDA_VERSION=4.8.2
#CONDA_VERSION_N=$(sed 's/\.//g' <<< "${CONDA_VERSION}")
#INSTALLED_CONDA_VERSION=$(conda -V | sed 's/conda //' | sed 's/\.//g')

# Check conda available
if ! [ -x "$(command -v conda)" ]; then
	echo -e "[${script}] Error: conda has not been installed or is not available on PATH"
	exit 1
fi

# Check conda version (rudamentary)
#if [ "${INSTALLED_CONDA_VERSION}" -lt "${CONDA_VERSION_N}" ]; then
#	echo -e "[${script}] Error - conda/miniconda is older than the required version"
#	echo -e "[${script}] Required: conda ${CONDA_VERSION} / Installed: $(conda -V)"
#	exit	
#fi

if ! [ -x "$(command -v wget)" ]; then
        echo -e "[${script}] Error: wget is not installed or is not available on PATH"
        exit 1
fi

if ! [ -x "$(command -v unzip)" ]; then
        echo -e "[${script}] Error: unzip is not installed or is not available on PATH"
        exit 1
fi

# Download additional files
echo -e "[${script}] Downloading additional data"
wget -qc https://zenodo.org/record/7573784/files/britroc_1_additional_data.zip
unzip -qq britroc_1_additional_data.zip

# Organising additional files
echo -e "[${script}] Moving additional data"
mv zenodo/component_parameters.rds copy_number_signatures/data/
mv zenodo/britroc_30kb_ds_absCopyNumber.rds absolute_POST_down_sampling/
mv zenodo/britroc_smoothed_copyNumbersSegmented.rds absolute_PRE_down_sampling/
mv zenodo/abcel_absoluteCN_30kbps_refitted.rds absolute_PRE_down_sampling/
mv zenodo/prenorm_downsampled_fitsbritroc_30kb_ds_absCopyNumber.rds absolute_PRE_down_sampling/
mv zenodo/clonality_results/ absolute_PRE_down_sampling/
# File clean up
echo -e "[${script}] Cleaning up additional data"
rm -rf zenodo/
rm britroc_1_additional_data.zip

# Check provided conda directory
if [ "$#" -lt 1 ]; then
	echo -e "[${script}] Error - conda/miniconda directory missing"
	echo -e "[${script}] Usage example './install_env.sh /home/user/miniconda3/'"
	exit 1
fi

# Set conda directory
CONDA_DIR=$1

echo -e "[${script}] Creating conda env"
# Create conda env
conda env create -q -f config/conda_env.yml
# Intitialise conda env
DIR=${CONDA_DIR}etc/profile.d/conda.sh
if [ -f "${DIR}" ]; then
	echo -e "[${script}] Initialising conda env"
	source ${DIR}
else
	echo -e "[${script}] Error: Unable to find conda intialisation script"
	echo -e "[${script}] Error: Make sure to provide the miniconda directory - e.g. '/home/user/miniconda3/'"
	echo -e "[${script}] Usage example './install_env.sh /home/user/miniconda3/'"
	exit
fi

# Activate env
echo -e "[${script}] Activating conda env"
conda activate ${ENV}
# Install modified QDNAseq pkg using devtools
echo -e "[${script}] Installing modified QDNAseq package"
Rscript -e 'devtools::install_github(repo = "markowetzlab/QDNAseqmod",quiet=TRUE,upgrade=FALSE)'
echo -e "[${script}] Installing CINSignatureQuantification package"
Rscript -e 'devtools::install_github(repo = "markowetzlab/cinsignaturequantification",ref="dev",quiet=TRUE,upgrade=FALSE)'
echo -e "[${script}] Installing logisticPCA package"
Rscript -e 'devtools::install_github(repo = "andland/logisticPCA",quiet=TRUE,upgrade=FALSE)'
# test all packages run and are available
#echo -e "[${script}] Testing package installation"
#R_LIB_PATH=$(Rscript resources/libpath.R)
#Rscript config/package_load.R
echo -e "[${script}] Conda env ready and all packages installed!"
echo -e "[${script}] Activate with 'conda activate ${ENV}'"

