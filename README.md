# metabinner_cami2b
GitHub repository for binning CAMI2b challenge dataset using MetaBinner

## <a name="started"></a>Getting Started

### <a name="docker"></a>Conda

We recommend using conda to run MetaBinner. Download [here](https://www.continuum.io/downloads)

### <a name="docker"></a>Obtain SolidBin and create an environment
After installing Anaconda (or miniconda), fisrt obtain MetaBinner:

```sh
git clone https://github.com/sufforest/SolidBin
```
Then simply create a metabinner_cami2b environment 

```sh
cd metabinner_cami2b
conda env create -f metabinner_cami2b_env.yaml
conda activate metabinner_cami2b_env
```

### <a name="docker"></a>Install checkM (python3 version) like this

(please make sure you have installed openssl)

```sh
cd CheckM-1.0.18
python setup.py install
```
Install checkM database:

CheckM relies on a number of precalculated data files which can be downloaded from https://data.ace.uq.edu.au/public/CheckM_databases/. (More details are available at https://github.com/Ecogenomics/CheckM/wiki/Installation#how-to-install-checkm):

```sh
mkdir <checkm_data_dir>
cd <checkm_data_dir>
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar xzf checkm_data_2015_01_16.tar.gz 
checkm data setRoot .
```
