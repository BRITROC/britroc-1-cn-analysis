# BriTROC
This repositiory contains the code required to reproduce the copy number analysis presented in the publication [The copy number and mutational landscape of recurrent ovarian high-grade serous carcinoma: the BriTROC-1 study](https://www.medrxiv.org/content/10.1101/2022.10.21.22280992v1) (Smith & Bradley _et al._ 2023 - under review).

## Copy number analysis
### Setup
Clone this git repository and `cd` into the directory
```
git clone https://github.com/Phil9S/britroc-cn-analysis.git
cd britroc-cn-analysis
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
#### Clinical DB login (To be reimplemented)

Currently patient and sample level data from the Britroc-1 study is accessed via an internally hosted [PostgreSQL database](https://github.com/TBradley27/britroc1_db) and as such access to the clinical data is restricted. Two requirements must be met to successfully access the clinical database utilised in this analysis pipeline.
- An connection to the CRUK CI internal network (either by local area or VPN).
- An available config file in the top-level directory containing valid credentials to log in to the postgreSQL server

The config file should be named `config.yml`, be placed in the top-level directory of this repository, and contain the following lines;
```
default:
  clinDB:
    dbname: 'britroc1'
    host: 'hostconn.cam.org'
    port: 5432
    user: 'username'
    password: 'pass phrase'
```
Where `host`, `username`, and `pass phrase` match the database host name, your log in, and pass phrase, respectively.

### Run analysis
Given all the previous steps were completed without issue, the entire copy number analysis can be reproduced by running the following bash script
```
./run_analysis_pipeline.sh
```
This will generate html markdown documents for all analyses performed as part of this publication excluding the compilation of main figures which is dependent on additional repositories.

### Author
[@phil9s](https://github.com/Phil9S)
