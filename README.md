# BriTROC-1 copy number analysis pipeline
This repositiory contains the code required to reproduce the copy number analysis presented in the publication [The copy number and mutational landscape of recurrent ovarian high-grade serous carcinoma](https://doi.org/10.1038/s41467-023-39867-7) (Smith & Bradley _et al._ 2023 - _Nature Communications_).

## Copy number analysis
### Setup
Clone this git repository and `cd` into the directory
```
git clone https://github.com/BRITROC/britroc-1-cn-analysis.git
cd britroc-1-cn-analysis
```
#### Install the conda environment
Provided conda is installed and available on the `$PATH` run the following line. This will downloaded the pre-processed data from [zenodo](https://zenodo.org/record/7573784#.ZC60rnbMJhE), install the required conda environment, and install packages which are not available via conda channels.
```
./install_env.sh {CONDA_INSTALL}
```
Where `{CONDA_INSTALL}` is the path to your conda directory, such as `/home/username/resources/miniconda/`.

#### Activate the newly installed environment
```
conda activate britroc
```

### Run analysis
Given all the previous steps were completed without issue, the entire copy number analysis can be reproduced by running the following bash script
```
./run_analysis_pipeline.sh
```
This will generate html markdown documents for all analyses performed as part of this publication excluding the compilation of main figures which is dependent on additional repositories.

### Author
[@phil9s](https://github.com/Phil9S)
