import os
from collections import defaultdict
import pickle
import sys

"""
#megahit short read
seed_file_prefix ='/home/wangzy/data/cami2019/cami2b/output/use_solidbin_scriptes_get_kmer_cov/rhimgCAMI2_short_read_pooled_megahit_map/output/new_ensemble_strategy/rhimgCAMI2_short_read_pooled_megahit_999.fa.bacar.seed'
cov_file ='/home/wangzy/data/cami2019/cami2b/output/use_solidbin_scriptes_get_kmer_cov/rhimgCAMI2_short_read_pooled_megahit_map/coverage_f999.tsv'
dict_out ='/home/wangzy/data/cami2019/cami2b/output/use_solidbin_scriptes_get_kmer_cov/rhimgCAMI2_short_read_pooled_megahit_map/output/new_ensemble_strategy/rhimgCAMI2_short_read_pooled_megahit_999.fa.bacar.seed.seed_contig_dict.pkl'
"""
seed_file_prefix =sys.argv[1]
marker_set = sys.argv[2]

dict_out = seed_file_prefix +'.seed_contig_dict.pkl'

if marker_set=='bacar':
    bacar_markers_num = 40
elif marker_set=='marker':
    bacar_markers_num = 107
else:
    print("The marker set is not included in this method.")


seed_contig_dict = defaultdict()
for i in range(bacar_markers_num): # 40 bacar markers
    file_name = seed_file_prefix + '.' + str(i)
    marker_id = os.path.basename(file_name)
    with open (file_name) as f:
        contigs_list = [x.strip() for x in f]
        seed_contig_dict[marker_id] = contigs_list



with open (dict_out, 'wb') as f: #
    pickle.dump(seed_contig_dict, f)


#with open (dict_out, 'rb') as f: #
#   t3 = pickle.load(f)
#  print(t3)