#######pooled
contig_file=$1
#/home/wangzy/data/cami2019/cami2b/output/use_solidbin_scriptes_get_kmer_cov/rhimgCAMI2_short_read_pooled_gsa_map/rhimgCAMI2_short_read_pooled_gsa_999.fa
output_dir=$2
#/home/wangzy/data/cami2019/cami2b/output/use_solidbin_scriptes_get_kmer_cov/rhimgCAMI2_short_read_pooled_gsa_map/output
coverage_profiles=$3
#${data_files_path}/coverage_f999.tsv
kmer_profile=$4
metabinner_path=$5
#/home/wangzy/tools/Metabinner2.0

contig_file_name=`basename ${contig_file}`
path=${output_dir}/metabinner_res/ensemble
############


Method_name=combo_sort_matches_maxaddsum_q_quality_filter_flag_False_updatesel_q_50.0_cont_weight_3_reprot_minq_10.0_add_matches_2_bins
components_path=${output_dir}/metabinner_res/component_result

mkdir ${output_dir}/metabinner_res

cd ${metabinner_path}
python metabinner.py \
--contig_file ${contig_file} \
--coverage_profiles ${coverage_profiles} \
--composition_profiles ${kmer_profile} \
--output ${output_dir}/metabinner_res/result.tsv \
--log ${output_dir}/metabinner_res/result.log \
--threads 40


##################
#get contig namelist information for each marker gene
mkdir ${output_dir}/metabinner_res/ensemble

for i in {0..39}
do
   echo $i
   perl ${metabinner_path}/auxiliary/get_whole_marker/getmarker_idx.pl ${contig_file}.bacar_marker.hmmout \
${contig_file} 1000 \
${output_dir}/metabinner_res/ensemble/${contig_file_name}.bacar.seed.${i} ${i}
done

for i in {0..106}
do
   echo $i
   perl ${metabinner_path}/auxiliary/get_whole_marker/getmarker_idx.pl ${contig_file}.marker.hmmout \
${contig_file} 1000 \
${output_dir}/metabinner_res/ensemble/${contig_file_name}.marker.seed.${i} ${i}
done


###
#get seed_contig_dict

python ${metabinner_path}/code_for_cami2b/gen_seed_contig_dict_bacar.py \
${output_dir}/metabinner_res/ensemble/${contig_file_name}.bacar.seed bacar

python ${metabinner_path}/code_for_cami2b/gen_seed_contig_dict_bacar.py \
${output_dir}/metabinner_res/ensemble/${contig_file_name}.marker.seed marker

#########
###the demo is suitable when the sample number is larger than 5
components_path=${output_dir}/metabinner_res/component_result

cp -r ${output_dir}/metabinner_res/intermediate_result \
${components_path}

cd ${components_path}
rm -rf kmeans_length_weight_X_*_result.tsv

ls *tsv > ${components_path}.res_namelist.tsv

cat ${components_path}.res_namelist.tsv | while read LINE
do
echo $LINE;
python ${metabinner_path}/code_for_cami2b/gen_bins_from_tsv.py \
-f ${contig_file} \
-r ${components_path}/${LINE} \
-o ${components_path}/${LINE}_bins
done

bash ${metabinner_path}/code_for_cami2b/gen_bin_dir_files.sh ${components_path}


#can be changed!!
mkdir ${output_dir}/metabinner_res/ensemble/em_res_X_t_logtrans
mkdir ${output_dir}/metabinner_res/ensemble/em_res_X_t_notrans
mkdir ${output_dir}/metabinner_res/ensemble/em_res_X_cov_logtrans
mkdir ${output_dir}/metabinner_res/ensemble/em_res_X_cov_notrans


python ${metabinner_path}/code_for_cami2b/ensemble_for_comps.py combo ${output_dir}/metabinner_res/ensemble/${contig_file_name}.bacar.seed.seed_contig_dict.pkl \
${output_dir}/metabinner_res/ensemble/${contig_file_name}.marker.seed.seed_contig_dict.pkl \
${components_path}.X_t_logtrans_bin_dirs.tsv \
${output_dir}/metabinner_res/ensemble/em_res_X_t_logtrans &

python ${metabinner_path}/code_for_cami2b/ensemble_for_comps.py combo ${output_dir}/metabinner_res/ensemble/${contig_file_name}.bacar.seed.seed_contig_dict.pkl \
${output_dir}/metabinner_res/ensemble/${contig_file_name}.marker.seed.seed_contig_dict.pkl \
${components_path}.X_t_notrans_bin_dirs.tsv \
${output_dir}/metabinner_res/ensemble/em_res_X_t_notrans &


python ${metabinner_path}/code_for_cami2b/ensemble_for_comps.py combo ${output_dir}/metabinner_res/ensemble/${contig_file_name}.bacar.seed.seed_contig_dict.pkl \
${output_dir}/metabinner_res/ensemble/${contig_file_name}.marker.seed.seed_contig_dict.pkl \
${components_path}.X_cov_logtrans_bin_dirs.tsv \
${output_dir}/metabinner_res/ensemble/em_res_X_cov_logtrans &

python ${metabinner_path}/code_for_cami2b/ensemble_for_comps.py combo ${output_dir}/metabinner_res/ensemble/${contig_file_name}.bacar.seed.seed_contig_dict.pkl \
${output_dir}/metabinner_res/ensemble/${contig_file_name}.marker.seed.seed_contig_dict.pkl \
${components_path}.X_cov_notrans_bin_dirs.tsv \
${output_dir}/metabinner_res/ensemble/em_res_X_cov_notrans &

wait

# em_mode = sys.argv[1] #'combo'
#bacar_dict_out = sys.argv[2]
#bacteria_dict_out = sys.argv[3]
#bin_dirs_file = sys.argv[4]
#output_dir= sys.argv[5]



#gen_callback

python ${metabinner_path}/code_for_cami2b/maxbin_to_bin_zy.py \
--paths ${path}/em_res_X_t_notrans/${Method_name}/ \
-o ${path}/em_res_X_t_notrans/${Method_name}_res.tsv

python ${metabinner_path}/code_for_cami2b/maxbin_to_bin_zy.py \
--paths ${path}/em_res_X_t_logtrans/${Method_name}/ \
-o ${path}/em_res_X_t_logtrans/${Method_name}_res.tsv

python ${metabinner_path}/code_for_cami2b/maxbin_to_bin_zy.py \
--paths ${path}/em_res_X_cov_logtrans/${Method_name}/ \
-o ${path}/em_res_X_cov_logtrans/${Method_name}_res.tsv

python ${metabinner_path}/code_for_cami2b/maxbin_to_bin_zy.py \
--paths ${path}/em_res_X_cov_notrans/${Method_name}/ \
-o ${path}/em_res_X_cov_notrans/${Method_name}_res.tsv



sed -i "s/,/\t/g" ${path}/em_res_X_t_notrans/${Method_name}_res.tsv
sed -i "s/,/\t/g" ${path}/em_res_X_t_logtrans/${Method_name}_res.tsv
sed -i "s/,/\t/g" ${path}/em_res_X_cov_logtrans/${Method_name}_res.tsv
sed -i "s/,/\t/g" ${path}/em_res_X_cov_notrans/${Method_name}_res.tsv


python ${metabinner_path}/code_for_cami2b/callback_contigs_after_ensemble.py \
--contig_file ${contig_file} \
--coverage_profiles ${coverage_profiles} \
--composition_profiles ${kmer_profile} \
--output ${output_dir}/metabinner_res/callback_contigs_after_ensemble.tsv \
--log ${output_dir}/metabinner_res/callback_contigs_after_ensemble.log \
--threads 40 --das_tool_result ${path}/em_res_X_t_notrans/${Method_name}_res.tsv \
--intermediate_result ${output_dir}/metabinner_res/intermediate_result/partial_seed_kmeans_bacar_marker_seed_length_weight_3quarter_X_t_notrans_result.tsv &
python ${metabinner_path}/code_for_cami2b/callback_contigs_after_ensemble.py \
--contig_file ${contig_file} \
--coverage_profiles ${coverage_profiles} \
--composition_profiles ${kmer_profile} \
--output ${output_dir}/metabinner_res/callback_contigs_after_ensemble.tsv \
--log ${output_dir}/metabinner_res/callback_contigs_after_ensemble.log \
--threads 40 --das_tool_result ${path}/em_res_X_t_logtrans/${Method_name}_res.tsv \
--intermediate_result ${output_dir}/metabinner_res/intermediate_result/partial_seed_kmeans_bacar_marker_seed_length_weight_3quarter_X_t_logtrans_result.tsv &
python ${metabinner_path}/code_for_cami2b/callback_contigs_after_ensemble.py \
--contig_file ${contig_file} \
--coverage_profiles ${coverage_profiles} \
--composition_profiles ${kmer_profile} \
--output ${output_dir}/metabinner_res/callback_contigs_after_ensemble.tsv \
--log ${output_dir}/metabinner_res/callback_contigs_after_ensemble.log \
--threads 40 --das_tool_result ${path}/em_res_X_cov_notrans/${Method_name}_res.tsv \
--intermediate_result ${output_dir}/metabinner_res/intermediate_result/partial_seed_kmeans_bacar_marker_seed_length_weight_3quarter_X_cov_notrans_result.tsv &
python ${metabinner_path}/code_for_cami2b/callback_contigs_after_ensemble.py \
--contig_file ${contig_file} \
--coverage_profiles ${coverage_profiles} \
--composition_profiles ${kmer_profile} \
--output ${output_dir}/metabinner_res/callback_contigs_after_ensemble.tsv \
--log ${output_dir}/metabinner_res/callback_contigs_after_ensemble.log \
--threads 40 --das_tool_result ${path}/em_res_X_cov_logtrans/${Method_name}_res.tsv \
--intermediate_result ${output_dir}/metabinner_res/intermediate_result/partial_seed_kmeans_bacar_marker_seed_length_weight_3quarter_X_cov_logtrans_result.tsv &

wait


######
#genbins

res_path=${path}/em_res_X_cov_notrans/${Method_name}_res.tsv.add_remained_after_dastool.tsv

python ${metabinner_path}/code_for_cami2b/gen_bins_from_tsv.py \
-f ${contig_file} \
-r ${res_path} \
-o ${res_path}_bins

res_path=${path}/em_res_X_cov_logtrans/${Method_name}_res.tsv.add_remained_after_dastool.tsv

python ${metabinner_path}/code_for_cami2b/gen_bins_from_tsv.py \
-f ${contig_file} \
-r ${res_path} \
-o ${res_path}_bins

res_path=${path}/em_res_X_t_logtrans/${Method_name}_res.tsv.add_remained_after_dastool.tsv

python ${metabinner_path}/code_for_cami2b/gen_bins_from_tsv.py \
-f ${contig_file} \
-r ${res_path} \
-o ${res_path}_bins

res_path=${path}/em_res_X_t_notrans/${Method_name}_res.tsv.add_remained_after_dastool.tsv

python ${metabinner_path}/code_for_cami2b/gen_bins_from_tsv.py \
-f ${contig_file} \
-r ${res_path} \
-o ${res_path}_bins


binsA=${path}/em_res_X_cov_notrans/${Method_name}_res.tsv.add_remained_after_dastool.tsv_bins
binsB=${path}/em_res_X_cov_logtrans/${Method_name}_res.tsv.add_remained_after_dastool.tsv_bins
binsC=${path}/em_res_X_t_logtrans/${Method_name}_res.tsv.add_remained_after_dastool.tsv_bins
binsD=${path}/em_res_X_t_notrans/${Method_name}_res.tsv.add_remained_after_dastool.tsv_bins

mkdir ${path}/${Method_name}
cd ${path}/${Method_name}

mkdir addcallback
cd addcallback

cp -r ${binsA} X_cov_notrans_${Method_name}_addcallback_bins
cp -r ${binsB} X_cov_logtrans_${Method_name}_addcallback_bins
cp -r ${binsC} X_t_logtrans_${Method_name}_addcallback_bins
cp -r ${binsD} X_t_notrans_${Method_name}_addcallback_bins

echo "X_t_logtrans,${binsC}" >> bins_dir.tsv
echo "X_t_notrans,${binsD}" >> bins_dir.tsv
echo "X_cov_logtrans,${binsB}" >> bins_dir.tsv
echo "X_cov_notrans,${binsA}" >> bins_dir.tsv
sed -i "s/,/\t/g" bins_dir.tsv

#############
#callback -> binning refiner
#binning refiner

binsA=X_cov_logtrans_${Method_name}_addcallback_bins
binsB=X_cov_notrans_${Method_name}_addcallback_bins
binsC=X_t_logtrans_${Method_name}_addcallback_bins
binsD=X_t_notrans_${Method_name}_addcallback_bins


python ${metabinner_path}/code_for_cami2b/binning_refiner_fa.py -1 ${binsA} -2 ${binsB} -3 ${binsC} -o Refined_ABC > Refined_ABC.out &
python ${metabinner_path}/code_for_cami2b/binning_refiner_fa.py -1 ${binsA} -2 ${binsB} -3 ${binsD} -o Refined_ABD > Refined_ABD.out &
python ${metabinner_path}/code_for_cami2b/binning_refiner_fa.py -1 ${binsA} -2 ${binsC} -3 ${binsD} -o Refined_ACD > Refined_ACD.out &
python ${metabinner_path}/code_for_cami2b/binning_refiner_fa.py -1 ${binsB} -2 ${binsC} -3 ${binsD} -o Refined_BCD > Refined_BCD.out &
python ${metabinner_path}/code_for_cami2b/binning_refiner_fa.py -1 ${binsA} -2 ${binsB} -o Refined_AB > Refined_AB.out &
python ${metabinner_path}/code_for_cami2b/binning_refiner_fa.py -1 ${binsA} -2 ${binsC} -o Refined_AC > Refined_AC.out &
python ${metabinner_path}/code_for_cami2b/binning_refiner_fa.py -1 ${binsA} -2 ${binsD} -o Refined_AD > Refined_AD.out &
python ${metabinner_path}/code_for_cami2b/binning_refiner_fa.py -1 ${binsB} -2 ${binsC} -o Refined_BC > Refined_BC.out &
python ${metabinner_path}/code_for_cami2b/binning_refiner_fa.py -1 ${binsB} -2 ${binsD} -o Refined_BD > Refined_BD.out &
python ${metabinner_path}/code_for_cami2b/binning_refiner_fa.py -1 ${binsC} -2 ${binsD} -o Refined_CD > Refined_CD.out &

wait

mv Refined_ABC/Refined Refined_ABC/Refined_ABC
mv Refined_ABD/Refined Refined_ABD/Refined_ABD
mv Refined_ACD/Refined Refined_ACD/Refined_ACD
mv Refined_BCD/Refined Refined_BCD/Refined_BCD

cd ${path}/${Method_name}

cd addcallback

cp -r bins_dir.tsv addrefined3comps_bins_dir.tsv
echo "Refined_ABC,${path}/${Method_name}/addcallback/Refined_ABC/Refined_ABC" >> addrefined3comps_bins_dir.tsv
echo "Refined_ABD,${path}/${Method_name}/addcallback/Refined_ABD/Refined_ABD" >> addrefined3comps_bins_dir.tsv
echo "Refined_ACD,${path}/${Method_name}/addcallback/Refined_ACD/Refined_ACD" >> addrefined3comps_bins_dir.tsv
echo "Refined_BCD,${path}/${Method_name}/addcallback/Refined_BCD/Refined_BCD" >> addrefined3comps_bins_dir.tsv

sed -i "s/,/\t/g" addrefined3comps_bins_dir.tsv


mv Refined_AB/Refined Refined_AB/Refined_AB
mv Refined_AC/Refined Refined_AC/Refined_AC
mv Refined_AD/Refined Refined_AD/Refined_AD
mv Refined_BC/Refined Refined_BC/Refined_BC
mv Refined_BD/Refined Refined_BD/Refined_BD
mv Refined_CD/Refined Refined_CD/Refined_CD

cp -r addrefined3comps_bins_dir.tsv addrefined2and3comps_bins_dir.tsv
echo "Refined_AB,${path}/${Method_name}/addcallback/Refined_AB/Refined_AB" >> addrefined2and3comps_bins_dir.tsv
echo "Refined_AC,${path}/${Method_name}/addcallback/Refined_AC/Refined_AC" >> addrefined2and3comps_bins_dir.tsv
echo "Refined_AD,${path}/${Method_name}/addcallback/Refined_AD/Refined_AD" >> addrefined2and3comps_bins_dir.tsv
echo "Refined_BC,${path}/${Method_name}/addcallback/Refined_BC/Refined_BC" >> addrefined2and3comps_bins_dir.tsv
echo "Refined_BD,${path}/${Method_name}/addcallback/Refined_BD/Refined_BD" >> addrefined2and3comps_bins_dir.tsv
echo "Refined_CD,${path}/${Method_name}/addcallback/Refined_CD/Refined_CD" >> addrefined2and3comps_bins_dir.tsv

sed -i "s/,/\t/g" addrefined2and3comps_bins_dir.tsv

###########
#ensemble_second_step
mkdir addrefined3comps
mkdir addrefined2and3comps


python ${metabinner_path}/code_for_cami2b/ensemble_second_step_mypipeline2.py \
${output_dir}/metabinner_res/ensemble/${contig_file_name}.bacar.seed.seed_contig_dict.pkl \
${output_dir}/metabinner_res/ensemble/${contig_file_name}.marker.seed.seed_contig_dict.pkl \
${path}/${Method_name}/addcallback/addrefined2and3comps_bins_dir.tsv \
${path}/${Method_name}/addcallback/addrefined2and3comps

second_method_name=my_pipeline2_sort_matches_maxaddsum_q_quality_filter_flag_True_dynamicsel_q_50.0_cont_weight_3_reprot_minq_10.0_sel_min_comp_0.0_sel_max_cont_100.0_add_matches_3_bins

python ${metabinner_path}/code_for_cami2b/maxbin_to_bin_zy.py \
--paths ${path}/${Method_name}/addcallback/addrefined2and3comps/${second_method_name} \
-o ${path}/${Method_name}/addcallback/addrefined2and3comps/${second_method_name}_res.tsv

cp -r ${path}/${Method_name}/addcallback/addrefined2and3comps/${second_method_name}_res.tsv ${output_dir}/metabinner_res/final_result_combo_my_pipeline2.tsv
sed -i "s/,/\t/g" ${output_dir}/metabinner_res/final_result_combo_my_pipeline2.tsv


echo "finished!"


