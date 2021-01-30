# metabinner_cami2b
GitHub repository for binning CAMI2b challenge dataset using MetaBinner

## <a name="started"></a>Getting Started

### <a name="docker"></a>Conda

We recommend using conda to run MetaBinner.

### <a name="docker"></a>Obtain codes and create an environment
After installing Anaconda (or miniconda), fisrt obtain MetaBinner:

```sh
git clone https://github.com/ziyewang/metabinner_cami2b
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

An example:
```sh
#Filter short contigs and generate kmer profiles:
python scripts/filter_tooshort_for_contig_file.py test_data/Sim40_20_ori/input/final_contigs.fa 999
python scripts/gen_kmer.py test_data/Sim40_20_ori/input/final_contigs.fa 999 4

#path to the input files for metabinner and the output dir:
contig_file=test_data/Sim40_20_ori/input/final_contigs_f1k_999.fa
kmer_files=test_data/Sim40_20_ori/input/kmer_4_f999.csv
coverage_profiles=test_data/Sim40_20_ori/input/Coverage_f1k.tsv
output_dir=test_data/Sim40_20_ori/output

mkdir ${output_dir}/metabinner_res

python metabinner.py \
--contig_file ${contig_file} \
--coverage_profiles ${coverage_profiles} \
--composition_profiles ${kmer_files} \
--output ${output_dir}/metabinner_res/result.tsv \
--log ${output_dir}/metabinner_res/result.log \
--threads 40
```

