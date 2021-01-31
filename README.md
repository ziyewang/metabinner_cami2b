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
python scripts/filter_tooshort_for_contig_file.py test_data/final_contigs.fa 999
python scripts/gen_kmer.py test_data/final_contigs.fa 999 4

#path to the input files for metabinner and the output dir:
contig_file=test_data/final_contigs_999.fa
kmer_files=test_data/kmer_4_f999.csv
coverage_profiles=test_data/Coverage_f1k.tsv
output_dir=test_data/output

mkdir ${output_dir}/metabinner_res

bash ${metabinner_path}/code_for_cami2b/metabinner_cami2b_pipeline_v1.2.sh ${contig_file} ${output_dir} ${coverage_profiles} ${kmer_profile} ${metabinner_path}


#The file "final_result_combo_my_pipeline2.tsv" in the "${output_dir}/metabinner_res" is the final output.
```


## <a name="preprocessing"></a>Contacts and bug reports
Please send bug reports or questions (such as the appropriate modes for your datasets) to
Ziye Wang: zwang17@fudan.edu.cn and Dr. Shanfeng Zhu: zhusf@fudan.edu.cn

## <a name="preprocessing"></a>References

[1] Lu, Yang Young, et al. "COCACOLA: binning metagenomic contigs using sequence COmposition, read CoverAge, CO-alignment and paired-end read LinkAge." Bioinformatics 33.6 (2017): 791-798.

[2] https://github.com/dparks1134/UniteM.

[3] Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. 2015. "CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes." Genome Research, 25: 1043–1055.

[4] Graham ED, Heidelberg JF, Tully BJ. (2017) "BinSanity: unsupervised clustering of environmental microbial assemblies using coverage and affinity propagation." PeerJ 5:e3035
