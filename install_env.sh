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
# test all packages run and are available
#echo -e "[${script}] Testing package installation"
#R_LIB_PATH=$(Rscript resources/libpath.R)
#Rscript config/package_load.R
echo -e "[${script}] Conda env ready and all packages installed!"
echo -e "Activate with `conda activate ${ENV}"
