# -*- coding: utf-8 -*-
# @Author  : ZWang
# @FileName: Metabinner.py
# scikit-learn == 0.22.1
# python 3.7

import numpy as np
import pandas as pd
import functools
import sys
import time
import gzip
import csv
import mimetypes
import os
import re

import logging
import argparse

from Bio import SeqIO
import scipy.sparse as sp
from sklearn.cluster.k_means_ import euclidean_distances, stable_cumsum, KMeans, check_random_state, row_norms

import shutil
import subprocess


logger = logging.getLogger('Metabinner final version not finished')

logger.setLevel(logging.INFO)

# logging
formatter = logging.Formatter('%(asctime)s - %(message)s')

console_hdr = logging.StreamHandler()
console_hdr.setFormatter(formatter)

logger.addHandler(console_hdr)

def arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('--contig_file', type=str, help=("The contigs file."))
    parser.add_argument('--coverage_profiles', type=str, help=(
        "The coverage profiles, containing a table where each row correspond to a contig, and each column correspond to a sample. All values are separated with tabs."))
    parser.add_argument('--composition_profiles', type=str, help=(
        "The composition profiles, containing a table where each row correspond to a contig, and each column correspond to the kmer composition of particular kmer. All values are separated with comma."))
    parser.add_argument('--output', type=str, help="The output file, storing the binning result.")
    parser.add_argument('--log', type=str, help="Specify where to store log file")
    parser.add_argument('--clusters', default=0, type=int,
                        help="Specify the number of clusters. If not specified, the cluster number is estimated by single-copy genes. If the number is small, the cluster number is estimated by single-copy genes ")
    parser.add_argument('--estimated_k', default=0, type=int,
                        help="Specify the number of clusters. If not specified, the cluster number is estimated by single-copy genes. If the number is small, the cluster number is estimated by single-copy genes ")
    parser.add_argument('--binscore', default=0.3, type=float,
                        help="Specify the score threshold for das_tool.")
    parser.add_argument('--threads', default=20, type=int,
                        help="the number of threads. default is 20.")
    parser.add_argument('--dataset_scale', type=str,default="large", help=(
        "The scale of the dataset (for bin number identification),large or small, default is large"))
    parser.add_argument('--search_engine_for_dastool', type=str, default="usearch", help=(
        "search engine for dastool, default is usearch. You can choose another engine with 'diamond' or 'blast'."))

    args = parser.parse_args()
    if not (args.contig_file and args.coverage_profiles and args.composition_profiles and args.output):
        parser.error(
            "Data is missing, add file(s) using --contig_file <contig_file> and/or --coverage_profiles <abund_profiles> and/or --composition_profiles <comp_profiles> and/or --output <out_file>")
        sys.exit(0)
    return args

def gen_X(com_file, cov_file):
    covHeader = pd.read_csv(cov_file, sep='\t', nrows=1)
    covMat = pd.read_csv(cov_file, sep='\t', usecols=range(1, covHeader.shape[1])).values
    namelist = pd.read_csv(cov_file, sep='\t', usecols=range(1)).values[:, 0]
    mapObj = dict(zip(namelist, range(len(namelist))))

    compositHeader = pd.read_csv(com_file, sep=',', nrows=1)
    shuffled_compositMat = pd.read_csv(com_file, sep=',', usecols=range(1, compositHeader.shape[1])).values
    shuffled_namelist = pd.read_csv(com_file, sep=',', usecols=range(1)).values[:, 0]

    covIdxArr = np.empty(len(mapObj), dtype=np.int)
    for contigIdx in range(len(shuffled_namelist)):
        if shuffled_namelist[contigIdx] in mapObj:
            covIdxArr[mapObj[shuffled_namelist[contigIdx]]] = contigIdx
    compositMat = shuffled_compositMat[covIdxArr]

    covMat = covMat + 1e-2
    covMat = covMat / covMat.sum(axis=0)[None, :]
    if covMat.shape[1] > 1:
        covMat = covMat / covMat.sum(axis=1)[:, None]
    compositMat = compositMat + 1
    compositMat = compositMat / compositMat.sum(axis=1)[:, None]
    X_t = np.hstack((covMat, compositMat)) # del * 1e1
    return X_t, namelist, mapObj, covMat, compositMat

def gen_second_covprofiles(another_cov_file, mapObj):
    second_cov_Header = pd.read_csv(another_cov_file, sep='\t', nrows=1)
    shuffled_second_cov = pd.read_csv(another_cov_file, sep='\t', usecols=range(1, second_cov_Header.shape[1])).values
    shuffled_namelist = pd.read_csv(another_cov_file, sep='\t', usecols=range(1)).values[:, 0]

    covIdxArr = np.empty(len(mapObj), dtype=np.int)
    for contigIdx in range(len(shuffled_namelist)):
        if shuffled_namelist[contigIdx] in mapObj:
            covIdxArr[mapObj[shuffled_namelist[contigIdx]]] = contigIdx
    second_cov_mat = shuffled_second_cov[covIdxArr]

    second_cov_mat = second_cov_mat + 1e-2
    second_cov_mat = second_cov_mat / second_cov_mat.sum(axis=0)[None, :]
    if second_cov_mat.shape[1] > 1:
        second_cov_mat = second_cov_mat / second_cov_mat.sum(axis=1)[:, None]

    return second_cov_mat

def gen_seed(contig_file, threads, marker_name="marker",quarter="3quarter"):
    fragScanURL = os.path.join(os.getcwd(), 'auxiliary', 'FragGeneScan1.19', 'run_FragGeneScan.pl')
    os.system("chmod 777 " + fragScanURL)
    os.system("chmod 777 "+ os.path.join(os.getcwd(), 'auxiliary', 'FragGeneScan1.19', 'FGS_gff.py'))
    os.system("chmod 777 " + os.path.join(os.getcwd(), 'auxiliary', 'FragGeneScan1.19', 'FragGeneScan'))

    hmmExeURL = os.path.join(os.getcwd(), 'auxiliary', 'hmmer-3.1b1', 'bin', 'hmmsearch')
    os.system("chmod 777 " + hmmExeURL)
    markerExeURL = os.path.join(os.getcwd(), 'auxiliary', 'test_getmarker_'+quarter+'.pl')
    os.system("chmod 777 " + markerExeURL)
    markerURL = os.path.join(os.getcwd(), 'auxiliary', marker_name + '.hmm')
    seedURL = contig_file +"." + marker_name +"."+quarter +".seed"
    fragResultURL = contig_file + ".frag.faa"
    hmmResultURL = contig_file + '.'+ marker_name + ".hmmout"

    if not (os.path.exists(fragResultURL)):
        fragCmd = fragScanURL + " -genome=" + contig_file + " -out=" + contig_file + ".frag -complete=0 -train=complete -thread=" + str(threads)+" 1>" + contig_file + ".frag.out 2>" + contig_file + ".frag.err"
        logger.info("exec cmd: " + fragCmd)
        os.system(fragCmd)

    if os.path.exists(fragResultURL):
        if not (os.path.exists(hmmResultURL)):
            hmmCmd = hmmExeURL + " --domtblout " + hmmResultURL + " --cut_tc --cpu " + str(threads) +" "+ markerURL + " " + fragResultURL + " 1>" + hmmResultURL + ".out 2>" + hmmResultURL + ".err"
            logger.info("exec cmd: " + hmmCmd)
            os.system(hmmCmd)

        if os.path.exists(hmmResultURL):
            if not (os.path.exists(seedURL)):
                markerCmd = markerExeURL + " " + hmmResultURL + " " + contig_file + " 1001 " + seedURL
                logger.info("exec cmd: " + markerCmd)
                os.system(markerCmd)

            if os.path.exists(seedURL):
                candK = file_len(seedURL)
            else:
                logger.info("markerCmd failed! Not exist: " + markerCmd)
                candK = 0
        else:
            logger.info("Hmmsearch failed! Not exist: " + hmmResultURL)
            sys.exit()
    else:
        logger.info("FragGeneScan failed! Not exist: " + fragResultURL)
        sys.exit()
    return candK

#estimate bin_number from candk
def estimate_bin_number(X_mat,candK,dataset_scale="large",len_weight=None):
    if dataset_scale=="small":
        candK = max(candK, 2)
        maxK = 4 * candK
        stepK = 2
    else:
        candK = max(candK, 2)
        maxK = 3 * candK
        stepK = 5
    bestK = candK
    bestSilVal = 0
    t = time.time()
    for k in range(candK, maxK, stepK):
        if k < len(X_mat):
            kmeans = KMeans(n_clusters=k, init='k-means++', random_state=7,n_init=30, n_jobs=-1)
            kmeans.fit(np.log(X_mat),sample_weight=len_weight)
            silVal = silhouette(np.log(X_mat), kmeans.cluster_centers_, kmeans.labels_,len_weight)
            logger.info("k:" + str(k) + "\tsilhouette:" + str(silVal) + "\telapsed time:" + str(time.time() - t))
            t = time.time()

            if silVal > bestSilVal:
                bestSilVal = silVal
                bestK = k
            else:
                break
        else:
            break
    candK = bestK + 2*stepK
    bestSilVal_2nd = 0
    for k in range(candK, maxK, stepK):
        if k < len(X_mat):
            kmeans = KMeans(n_clusters=k, init='k-means++', random_state=7,n_init=30, n_jobs=-1)
            kmeans.fit(np.log(X_mat),sample_weight=len_weight)
            silVal_2nd = silhouette(np.log(X_mat), kmeans.cluster_centers_, kmeans.labels_,len_weight)
            logger.info("k:" + str(k) + "\tsilhouette:" + str(silVal_2nd) + "\telapsed time:" + str(time.time() - t))
            t = time.time()
            if silVal_2nd > bestSilVal_2nd:
                bestSilVal_2nd = silVal_2nd
                bestK = k
            else:
                break
        else:
            break
    if bestSilVal_2nd > bestSilVal:
        bestSilVal = bestSilVal_2nd
    else:
        bestK = candK - 2*stepK
    logger.info("bestk:" + str(bestK) + "\tsilVal:" + str(bestSilVal))
    return bestK


def silhouette(X, W, label,len_weight):
    X_colsum = np.sum(X ** 2, axis=1)
    X_colsum = X_colsum.reshape(len(X_colsum), 1)
    W_colsum = np.sum(W ** 2, axis=1)
    W_colsum = W_colsum.reshape(len(W_colsum), 1)

    Dsquare = np.tile(X_colsum, (1, W.shape[0])) + np.tile(W_colsum.T, (X.shape[0], 1)) - 2 * X.dot(W.T)
    # avoid error caused by accuracy
    Dsquare[Dsquare < 0] = 0
    D = np.sqrt(Dsquare)
    aArr = D[np.arange(D.shape[0]), label]
    D[np.arange(D.shape[0]), label] = np.inf
    bArr = np.min(D, axis=1)
    tmp = (bArr - aArr) / np.maximum(aArr, bArr)
    return np.average(tmp,weights=len_weight)



def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def get_length(fastx_file):
    file_type = mimetypes.guess_type(fastx_file)[1]
    if file_type == 'gzip':
        f = gzip.open(fastx_file, "rt")
    elif not file_type:
        f = open(fastx_file, "rt")
    else:
        raise RuntimeError("Unknown type of file: '{}".format(fastx_file))
    length = {}
    if os.path.getsize(fastx_file) == 0:
        return length
    file_format = None
    line = f.readline()
    if line.startswith('@'):
        file_format = "fastq"
    elif line.startswith(">"):
        file_format = "fasta"
    f.seek(0)
    if not file_format:
        raise RuntimeError("Invalid sequence file: '{}".format(fastx_file))
    for seq_record in SeqIO.parse(f, file_format):
        length[seq_record.id] = len(seq_record.seq)

    f.close()
    return length


def gen_seed_idx(seedURL,contig_id_list):
    seed_list = []
    with open(seedURL) as f:
        for line in f:
            if line.rstrip('\n') in contig_id_list:
                seed_list.append(line.rstrip('\n'))
    name_map = dict(zip(contig_id_list, range(len(contig_id_list))))
    seed_idx = [name_map[seed_name] for seed_name in seed_list]
    return seed_idx


def save_result(result, filepath, namelist):
    filedir, filename = os.path.split(filepath)
    if not filename:
        filename = "result.tsv"
    if not os.path.exists(filedir):
        os.makedirs(filedir)
    f = open(filepath, 'w')
    for contigIdx in range(len(result)):
        f.write(namelist[contigIdx] + "\t" + str(result[contigIdx].item(0)) + "\n")
    f.close()

def save_ensemble_result(result, filepath, namelist):
    filedir, filename = os.path.split(filepath)
    if not filename:
        filename = "result.tsv"
    if not os.path.exists(filedir):
        os.makedirs(filedir)
    f = open(filepath, 'w')
    for contigIdx in range(len(result)):
        f.write(namelist[contigIdx] + "\t" + str(result[contigIdx]) + "\n")
    f.close()


# change from sklearn.cluster.kmeans
def partial_seed_init(X, n_clusters, random_state, seed_idx, n_local_trials=None):
    print('Using partial seed')

    random_state = check_random_state(random_state)
    x_squared_norms = row_norms(X, squared=True)

    n_samples, n_features = X.shape

    centers = np.empty((n_clusters, n_features), dtype=X.dtype)

    # Set the number of local seeding trials if none is given
    if n_local_trials is None:
        # This is what Arthur/Vassilvitskii tried, but did not report
        # specific results for other than mentioning in the conclusion
        # that it helped.
        n_local_trials = 2 + int(np.log(n_clusters))

    # Pick first center randomly

    center_id = seed_idx[0]

    if sp.issparse(X):
        centers[0] = X[center_id].toarray()
    else:
        centers[0] = X[center_id]

    # Initialize list of closest distances and calculate current potential
    closest_dist_sq = euclidean_distances(
        centers[0, np.newaxis], X, Y_norm_squared=x_squared_norms,
        squared=True)

    for c, center_id in enumerate(seed_idx[1:], 1):
        if sp.issparse(X):
            centers[c] = X[center_id].toarray()
        else:
            centers[c] = X[center_id]
        closest_dist_sq = np.minimum(closest_dist_sq,
                                     euclidean_distances(
                                         centers[c, np.newaxis], X, Y_norm_squared=x_squared_norms,
                                         squared=True))
    current_pot = closest_dist_sq.sum()

    # Pick the remaining n_clusters-1 points
    for c in range(len(seed_idx), n_clusters):
        # Choose center candidates by sampling with probability proportional
        # to the squared distance to the closest existing center
        rand_vals = random_state.random_sample(n_local_trials) * current_pot
        candidate_ids = np.searchsorted(stable_cumsum(closest_dist_sq),
                                        rand_vals)
        # XXX: numerical imprecision can result in a candidate_id out of range
        np.clip(candidate_ids, None, closest_dist_sq.size - 1,
                out=candidate_ids)

        # Compute distances to center candidates
        distance_to_candidates = euclidean_distances(
            X[candidate_ids], X, Y_norm_squared=x_squared_norms, squared=True)

        # Decide which candidate is the best
        best_candidate = None
        best_pot = None
        best_dist_sq = None
        for trial in range(n_local_trials):
            # Compute potential when including center candidate
            new_dist_sq = np.minimum(closest_dist_sq,
                                     distance_to_candidates[trial])
            new_pot = new_dist_sq.sum()

            # Store result if it is the best local trial so far
            if (best_candidate is None) or (new_pot < best_pot):
                best_candidate = candidate_ids[trial]
                best_pot = new_pot
                best_dist_sq = new_dist_sq

        # Permanently add best center candidate found in local tries
        if sp.issparse(X):
            centers[c] = X[best_candidate].toarray()
        else:
            centers[c] = X[best_candidate]
        current_pot = best_pot
        closest_dist_sq = best_dist_sq

    return centers

def seed_kmeans_combo(seed_idx, output, X_mat, bin_number, prefix, length_weight, marker_name="marker1",quarter="3quarter"):

    # run partial seed kmeans marker1_seed length weight
    logger.info("run partial seed kmeans "+marker_name+ " seed length weight with:\t" + quarter+'_'+prefix)
    output_temp = os.path.dirname(
        output) + '/intermediate_result' + '/partial_seed_kmeans_'+marker_name+'_seed_length_weight_' +quarter+'_'+ prefix + '_result.tsv'
    if not (os.path.exists(output_temp)):
        km = KMeans(n_clusters=bin_number, n_jobs=-1, random_state=7, n_init=30,
                    init=functools.partial(partial_seed_init, seed_idx=seed_idx))
        km.fit(X_mat,sample_weight=length_weight)
        idx = km.labels_
        save_result(idx, output_temp, namelist)

def my_kmeans(X_mat, namelist, bin_number, marker1_seed_num, bacar_marker_seed_num, length_weight, output, prefix='X_t_notrans',quarter="3quarter",contig_file=None):

    if marker1_seed_num > 0:
        seed_marker1_url = contig_file + ".marker"+"."+quarter+".seed"
        seed_marker1_idx = gen_seed_idx(seed_marker1_url, contig_id_list=namelist)
    if bacar_marker_seed_num > 0:
        seed_bacar_marker_url = contig_file + ".bacar_marker"+"."+quarter+".seed"
        seed_bacar_marker_idx = gen_seed_idx(seed_bacar_marker_url, contig_id_list=namelist)

    #run kmeans length weight
    logger.info("run kmeans length weight with:\t"+prefix)
    output_temp = os.path.dirname(output) + '/intermediate_result' + '/kmeans_length_weight_'+prefix+'_result.tsv'
    if not (os.path.exists(output_temp)):
        km = KMeans(n_clusters=bin_number, init='k-means++', n_jobs=-1, n_init=30, random_state=7)
        km.fit(X_mat, sample_weight=length_weight)  # add log transform
        idx = km.labels_
        save_result(idx, output_temp, namelist)

    # run kmeans with marker1 seed
    if marker1_seed_num > 0:
        seed_kmeans_combo(seed_marker1_idx,output, X_mat, bin_number, prefix,length_weight,quarter=quarter)
    # run kmeans with bacarmarker seed
    if bacar_marker_seed_num > 0:
        seed_kmeans_combo(seed_bacar_marker_idx,output, X_mat, bin_number, prefix,length_weight,marker_name="bacar_marker",quarter=quarter)


def run_dastool(intermediate_result_dir, output, dir_name="/das_tool_output",search_engine=None,partial_flag=None):
    result_files = os.listdir(intermediate_result_dir)
    result_files_del_ori=[]
    for i in range(len(result_files)):
        if 'kmeans_length' not in result_files[i]:
            result_files_del_ori.append(result_files[i])

    das_tool_dir = os.path.dirname(output) + dir_name
    if not (os.path.exists(das_tool_dir)):
        os.mkdir(das_tool_dir)
        if partial_flag:
            result_files_partial=[]
            for i in range(len(result_files_del_ori)):
                if partial_flag in result_files_del_ori[i]:
                    result_files_partial.append(result_files_del_ori[i])

            input_str = ''
            for i in range(len(result_files_partial)):
                if (i < (len(result_files_partial) - 1)):
                    input_str = input_str + intermediate_result_dir + '/' + result_files_partial[i] + ','
                else:
                    input_str = input_str + intermediate_result_dir + '/' + result_files_partial[i]

            label_str = ''
            for i in range(len(result_files_partial)):
                if (i < (len(result_files_partial) - 1)):
                    label_str = label_str + result_files_partial[i] + ','
                else:
                    label_str = label_str + result_files_partial[i]
        else:
            input_str = ''
            for i in range(len(result_files_del_ori)):
                if (i < (len(result_files_del_ori) - 1)):
                    input_str = input_str + intermediate_result_dir + '/' + result_files_del_ori[i] + ','
                else:
                    input_str = input_str + intermediate_result_dir + '/' + result_files_del_ori[i]

            label_str = ''
            for i in range(len(result_files_del_ori)):
                if (i < (len(result_files_del_ori) - 1)):
                    label_str = label_str + result_files_del_ori[i] + ','
                else:
                    label_str = label_str + result_files_del_ori[i]
        if partial_flag:
            das_toolCmd = (
                    "DAS_Tool -i " + input_str + " -l " + label_str + " -c " + contig_file + " -o " + das_tool_dir + "/das_tool_goodbins --threads " + str(
                threads) + " --score_threshold " + str(
                binscore) + " --search_engine " + str(
                search_engine) + " --proteins " +os.path.dirname(das_tool_dir)+"/das_tool_output_all/das_tool_goodbins_proteins.faa")
        else:
            das_toolCmd = (
                    "DAS_Tool -i " + input_str + " -l " + label_str + " -c " + contig_file + " -o " + das_tool_dir + "/das_tool_goodbins --threads " + str(
                threads) + " --score_threshold " + str(
                binscore) + " --search_engine " + str(
                search_engine))

        logger.info("exec cmd: " + das_toolCmd)
        os.system(das_toolCmd)
    else:
        logger.info("das_tool output exist")

def save_result_refine(result, filepath, namelist, unclassified_contigs_id_number):
    filedir, filename = os.path.split(filepath)
    if not filename:
        filename = "result.tsv"
    if not os.path.exists(filedir):
        os.makedirs(filedir)
    f = open(filepath, 'w')
    for Idx in range(len(result)):
        f.write(namelist[unclassified_contigs_id_number[Idx]] + "\t" + str(result[Idx].item(0)) + "\n")
    f.close()

def gen_remained_fasta_file(fastafile, remained_contig_id, outputdir, prefix_str):
    # read fasta file
    logger.info("Processing file:\t{}".format(fastafile))
    sequences = {}
    if fastafile.endswith("gz"):
        with gzip.open(fastafile, 'r') as f:
            for line in f:
                line = str(line, encoding="utf-8")
                if line.startswith(">"):
                    if " " in line:
                        seq, others = line.split(' ', 1)
                        sequences[seq] = ""
                    else:
                        seq = line.rstrip("\n")
                        sequences[seq] = ""
                else:
                    sequences[seq] += line.rstrip("\n")
    else:
        with open(fastafile, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    if " " in line:
                        seq, others = line.split(' ', 1)
                        sequences[seq] = ""
                    else:
                        seq = line.rstrip("\n")
                        sequences[seq] = ""
                else:
                    sequences[seq] += line.rstrip("\n")
    dic = {}
    cluster_name = 'remained'
    for contig_name in remained_contig_id:
        try:
            dic[cluster_name].append(contig_name)
        except:
            dic[cluster_name] = []
            dic[cluster_name].append(contig_name)
    logger.info("Writing bins:\t{}".format(outputdir))
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    bin_name = 0
    for _, cluster in dic.items():
        binfile = os.path.join(outputdir, "{}_{}.fa".format(prefix_str, bin_name))
        with open(binfile, "w") as f:
            for contig_name in cluster:
                contig_name = ">" + contig_name
                try:
                    sequence = sequences[contig_name]
                except:
                    bin_name += 1
                    continue
                f.write(contig_name + "\n")
                f.write(sequence + "\n")
                bin_name += 1

# change from binsanity
def checkm_analysis(file_, suffix_str, output):
    file_object = open(file_, 'r')
    lines = file_object.readlines()
    file_deal = open(file_ + '_deal.txt', 'w')
    try:
        for line in lines:
            if line.startswith(' '):
                print(line)
                file_deal.writelines(line)
    finally:
        file_object.close()
    file_deal.close()  # 要不会影响读�?

    checkm = list(csv.reader(open(file_ + '_deal.txt', 'r')))
    new = []
    for list_ in checkm:
        x = re.sub(' +', ' ', str(re.split(r'\t+', list_[0].rstrip('\t'))))
        new.append(x)

    del new[0]

    checkm_info_list = [list_.strip("['']") for list_ in new]
    checkm_info_list = [x.split() for x in checkm_info_list]

    good_bins = []
    High_completion_high_contamination = []
    # low_completion=[]
    others = []

    for list_ in checkm_info_list:
        if ((float(list_[12]) > 70 and (float(list_[13]) < 15)) or (float(list_[12]) - 5 * float(list_[13])) > 50):
            good_bins.append(list_[0])
        elif (float(list_[12]) > 70 and (float(list_[13]) > 50)):
            High_completion_high_contamination.append(list_[0])
        else:
            others.append(list_[0])

    if os.path.isdir(os.path.dirname(output) + "/good_bins") is False:
        os.makedirs(os.path.dirname(output) + "/good_bins")
    if os.path.isdir(os.path.dirname(output) + "/High_completion_high_contamination") is False:
        os.makedirs(os.path.dirname(output) + "/High_completion_high_contamination")
    if os.path.isdir(os.path.dirname(output) + "/others") is False:
        os.makedirs(os.path.dirname(output) + "/others")

    for name in good_bins:
        shutil.move((os.path.dirname(output) + '/' + str(name) + suffix_str),
                    os.path.dirname(output) + "/good_bins")
    for name in High_completion_high_contamination:
        shutil.move((os.path.dirname(output) + '/' + str(name) + suffix_str),
                    os.path.dirname(output) + "/High_completion_high_contamination")
    for name in others:
        shutil.move((os.path.dirname(output) + '/' + str(name) + suffix_str), os.path.dirname(output) + "/others")

def gen_bins(fastafile, resultfile, outputdir, prefix_str):
    # read fasta file
    logger.info("Processing file:\t{}".format(fastafile))
    sequences = {}
    if fastafile.endswith("gz"):
        with gzip.open(fastafile, 'r') as f:
            for line in f:
                line = str(line, encoding="utf-8")
                if line.startswith(">"):
                    if " " in line:
                        seq, others = line.split(' ', 1)
                        sequences[seq] = ""
                    else:
                        seq = line.rstrip("\n")
                        sequences[seq] = ""
                else:
                    sequences[seq] += line.rstrip("\n")
    else:
        with open(fastafile, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    if " " in line:
                        seq, others = line.split(' ', 1)
                        sequences[seq] = ""
                    else:
                        seq = line.rstrip("\n")
                        sequences[seq] = ""
                else:
                    sequences[seq] += line.rstrip("\n")
    logger.info("Reading Map:\t{}".format(resultfile))
    dic = {}
    with open(resultfile, "r") as f:
        for line in f:
            contig_name, cluster_name = line.strip().split('\t')  # change from split(',')
            try:
                dic[cluster_name].append(contig_name)
            except:
                dic[cluster_name] = []
                dic[cluster_name].append(contig_name)
    logger.info("Writing bins:\t{}".format(outputdir))
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    bin_name = 0
    for _, cluster in dic.items():
        binfile = os.path.join(outputdir, "{}_{}.bin".format(prefix_str, bin_name))
        with open(binfile, "w") as f:
            for contig_name in cluster:
                contig_name = ">" + contig_name
                try:
                    sequence = sequences[contig_name]
                except:
                    bin_name += 1
                    continue
                f.write(contig_name + "\n")
                f.write(sequence + "\n")
                bin_name += 1


def recluster_hh_bins(high_com_p_high_cont_path, mapObj, X_t, length_weight, namelist):
    hh_files = os.listdir(high_com_p_high_cont_path)
    for file_ in hh_files:
        if file_.endswith('.bin'):
            hh_contigs_id = []
            hh_contig_file = high_com_p_high_cont_path + '/' + file_
            for seq_record in SeqIO.parse(hh_contig_file, "fasta"):
                hh_contigs_id.append(seq_record.id)
            hh_contigs_id_number = [mapObj[x] for x in hh_contigs_id]
            X_t_hh_unclustered = X_t[hh_contigs_id_number]
            hh_weight = []
            for i in range(len(hh_contigs_id_number)):
                hh_weight.append(length_weight[hh_contigs_id_number[i]])

            seed_hh_num = gen_seed(hh_contig_file, threads, marker_name="bacar_marker")
            bin_number = estimate_bin_number(X_t_hh_unclustered, seed_hh_num, dataset_scale="small",len_weight=hh_weight)
            #bin_number = gen_bestk(hh_contig_file, X_t_hh_unclustered, 0,hh_weight)


            # seedurl may not exits??
            seedURL = hh_contig_file + ".bacar_marker.3_quarter.seed"
            # global seed_idx
            if os.path.exists(seedURL):
                seed_list = []
                with open(seedURL) as f:
                    for line in f:
                        if line.rstrip('\n') in namelist:
                            seed_list.append(line.rstrip('\n'))
                name_map = dict(zip(hh_contigs_id, range(len(hh_contigs_id))))
                seed_idx = [name_map[seed_name] for seed_name in seed_list]
                km = KMeans(n_clusters=bin_number, n_jobs=-1, n_init=30, random_state=7,
                            init=functools.partial(partial_seed_init, seed_idx=seed_idx))
            else:
                km = KMeans(n_clusters=bin_number, n_jobs=-1, n_init=30, random_state=7)

            km.fit(X_t_hh_unclustered, sample_weight=hh_weight)
            idx = km.labels_
            save_result_refine(idx, hh_contig_file + ".reclustered.tsv",
                               namelist, hh_contigs_id_number)
            gen_bins(hh_contig_file, hh_contig_file + ".reclustered.tsv",
                     os.path.dirname(high_com_p_high_cont_path) + '/good_bins', file_ + "_reclustered")

def recluster_other_contigs(not_clustered_path, X_t, namelist, mapObj, length_weight):
    files = os.listdir(not_clustered_path)
    other_contig_file = not_clustered_path + '/init_unclustered_contigs.fa'
    ofile = open(other_contig_file, 'w')
    # 遍历读取所有文件，并写入到输出文件
    for fr in files:
        if fr != 'init_unclustered_contigs.fa':
            for txt in open(not_clustered_path + '/' + fr, 'r'):
                ofile.write(txt)
    ofile.close()

    unclassified_contigs_id = []
    for seq_record in SeqIO.parse(other_contig_file, "fasta"):
        unclassified_contigs_id.append(seq_record.id)
    unclassified_contigs_id_number = [mapObj[x] for x in unclassified_contigs_id]
    X_t_unclustered = X_t[unclassified_contigs_id_number]

    unclassified_contigs_weight = []
    for i in range(len(unclassified_contigs_id_number)):
        unclassified_contigs_weight.append(length_weight[unclassified_contigs_id_number[i]])

    seed_o_num = gen_seed(other_contig_file, threads, marker_name="bacar_marker")
    bin_number = estimate_bin_number(X_t_unclustered, seed_o_num, dataset_scale="small",len_weight=unclassified_contigs_weight)
    logger.info("bin_number for other contigs: %d", bin_number)

    seedURL = other_contig_file +  ".bacar_marker.3_quarter.seed"
    # global seed_idx
    if os.path.exists(seedURL):
        seed_list = []
        with open(seedURL) as f:
            for line in f:
                if line.rstrip('\n') in namelist:
                    seed_list.append(line.rstrip('\n'))
        name_map = dict(zip(unclassified_contigs_id, range(len(unclassified_contigs_id))))
        seed_idx = [name_map[seed_name] for seed_name in seed_list]
        km = KMeans(n_clusters=bin_number, n_jobs=-1, n_init=30, random_state=7,
                    init=functools.partial(partial_seed_init, seed_idx=seed_idx))
    else:
        km = KMeans(n_clusters=bin_number, n_jobs=-1, n_init=30, random_state=7)

    logger.info("Start bin the other bins.")
    #
    km.fit(X_t_unclustered, sample_weight=unclassified_contigs_weight)
    idx = km.labels_
    not_clustered_path_output = not_clustered_path + 'reclustered_result.tsv'
    save_result_refine(idx, not_clustered_path_output, namelist, unclassified_contigs_id_number)
    gen_bins(other_contig_file, not_clustered_path_output, os.path.dirname(not_clustered_path) + '/good_bins',
             "reclustered")

def convert(paths, output_file):
    files = os.listdir(paths)
    fasta_files = []
    for file in files:
        if file.endswith(('.fasta', '.fa', '.fna', '.bin')):
            fasta_files.append(file)
    with open(output_file, 'w') as write_handler:
        for bin_id, fasta_file in enumerate(fasta_files):
            for sequence_id in read_fasta_file(paths + '/' + fasta_file):
                write_handler.write("%s\t%s\n" % (sequence_id, bin_id))  # change from ","

def read_fasta_file(fasta_file):
    with open(fasta_file, 'r') as read_handler:
        for line in read_handler:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                yield line[1:]

def checkm_post_precess(kmeans_initial_output,kmeans_initial_output_dir,prefix=None):
    if not os.path.exists(kmeans_initial_output_dir):
        os.mkdir(kmeans_initial_output_dir)
        gen_bins(contig_file, kmeans_initial_output, kmeans_initial_output_dir,
                 "kmeans_length_weight_X_t_notrans_result")
        checkm_out_dir = kmeans_initial_output_dir + '/checkm_out'
        os.mkdir(checkm_out_dir)
        output = args.output

        checkm_file = checkm_out_dir + "/checkm_analysis_init.txt"
        checkm_file_output = open((checkm_file), "w")
        threads = 20
        subprocess.call(["checkm", "lineage_wf", "-x", "bin", "-t", str(threads), kmeans_initial_output_dir,
                         checkm_out_dir], stdout=checkm_file_output)

        suffix_str = '.bin'
        checkm_analysis(checkm_file, suffix_str, checkm_out_dir)
        #
        goodbin_path = kmeans_initial_output_dir + '/good_bins'

        # 处理high_com_p_high_cont文件
        logger.info("Recluster the contigs from high_com_p_high_cont bins")
        high_com_p_high_cont_path = kmeans_initial_output_dir + "/High_completion_high_contamination"
        recluster_hh_bins(high_com_p_high_cont_path, mapObj, X_t, length_weight, namelist)

        # recluster other contigs
        logger.info("Recluster other contigs.")
        not_clustered_path = kmeans_initial_output_dir + "/others"
        recluster_other_contigs(not_clustered_path, X_t, namelist, mapObj, length_weight)

        convert(kmeans_initial_output_dir + '/good_bins', os.path.dirname(
            args.output) + '/kmeans_'+prefix+'_seed_partial_ori_with_postprocess.tsv')
        convert(kmeans_initial_output_dir + '/good_bins', intermediate_result_dir+ '/kmeans_' + prefix + '_seed_partial_ori_with_postprocess.tsv')

def get_unique_bin_number(result_path):
    clustered_contig_id = []
    result_id = []
    with open(result_path) as f:
        for line in f:
            clustered_contig_id.append(line.rstrip('\n').split('\t')[0])
            result_id.append(line.rstrip('\n').split('\t')[1])
    unique_binnumber = len(np.unique(result_id))
    return unique_binnumber, clustered_contig_id, result_id

if __name__ == '__main__':
    args = arguments()

    if args.log:
        handler = logging.FileHandler(args.log)
        handler.setLevel(logging.INFO)
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    logger.info("Input arguments:")
    logger.info("Contig_file:\t" + args.contig_file)
    logger.info("Coverage_profiles:\t" + args.coverage_profiles)
    logger.info("Composition_profiles:\t" + args.composition_profiles)
    logger.info("Output file path:\t" + args.output)
    logger.info("Predefined Clusters:\t" + (str(args.clusters) if args.clusters > 0 else "Auto"))
    #logger.info("Another_coverage_profiles:\t" + args.another_coverage_profiles if args.another_coverage_profiles else None)
    #logger.info("Binscore threshold for Das_tool:\t" + str(args.binscore))
    logger.info("The number of threads:\t" + str(args.threads))

    com_file = args.composition_profiles
    cov_file = args.coverage_profiles

    X_t, namelist, mapObj, X_cov, X_com = gen_X(com_file, cov_file)

    contigNum = X_t.shape[0]
    contig_file = args.contig_file
    output = args.output

    logger.info("The number of contigs:\t" + str(contigNum))

    clusters = args.clusters
    binscore = args.binscore
    threads = args.threads

    logger.info("gen marker seed")
    marker1_1quarter_seed_num = gen_seed(contig_file, threads, marker_name="marker",quarter="1quarter")
    logger.info("marker1_1quarter_seed_num:\t"+str(marker1_1quarter_seed_num))
    marker1_2quarter_seed_num = gen_seed(contig_file, threads, marker_name="marker",quarter="2quarter")
    logger.info("marker1_2quarter_seed_num:\t"+str(marker1_2quarter_seed_num))
    marker1_3quarter_seed_num = gen_seed(contig_file, threads, marker_name="marker",quarter="3quarter")
    logger.info("marker1_3quarter_seed_num:\t"+str(marker1_3quarter_seed_num))

    logger.info("gen bacar marker seed")
    bacar_marker_1quarter_seed_num = gen_seed(contig_file, threads, marker_name="bacar_marker",quarter="1quarter")
    logger.info("bacar_marker_1quarter_seed_num:\t" + str(bacar_marker_1quarter_seed_num))
    bacar_marker_2quarter_seed_num = gen_seed(contig_file, threads, marker_name="bacar_marker",quarter="2quarter")
    logger.info("bacar_marker_2quarter_seed_num:\t" + str(bacar_marker_2quarter_seed_num))
    bacar_marker_3quarter_seed_num = gen_seed(contig_file, threads, marker_name="bacar_marker",quarter="3quarter")
    logger.info("bacar_marker_3quarter_seed_num:\t" + str(bacar_marker_3quarter_seed_num))

    logger.info("start calculate contig length")
    lengths = get_length(contig_file)
    length_weight = []
    for seq_id in namelist:
        length_weight.append(lengths[seq_id])

    #set k0 using the large seed number from the two single-copy marker sets
    if args.estimated_k:
        bin_number = args.estimated_k
    else:
        dataset_scale = args.dataset_scale
        candK = max(marker1_3quarter_seed_num, bacar_marker_3quarter_seed_num) + 1
        logger.info("start estimate_bin_number")
        bin_number = estimate_bin_number(X_t, candK, dataset_scale=dataset_scale,len_weight=length_weight)
        logger.info("estimated_bin_number:\t" + str(bin_number))


    if args.clusters:
        bin_number = max(clusters,bin_number) #只有当用户指定k大于估计k时才会被采用

    intermediate_result_dir =os.path.dirname(output) + '/intermediate_result'
    if not (os.path.exists(intermediate_result_dir)):
        os.mkdir(intermediate_result_dir)
        # prefix denote which feature
        # X_t
        my_kmeans(X_t, namelist, bin_number, marker1_1quarter_seed_num, bacar_marker_1quarter_seed_num, length_weight, output,
                  prefix='X_t_notrans',quarter="1quarter",contig_file=contig_file)

        my_kmeans(np.log(X_t), namelist, bin_number, marker1_1quarter_seed_num, bacar_marker_1quarter_seed_num, length_weight, output,
                  prefix='X_t_logtrans',quarter="1quarter",contig_file=contig_file)

        my_kmeans(X_t, namelist, bin_number, marker1_2quarter_seed_num, bacar_marker_2quarter_seed_num, length_weight, output,
                  prefix='X_t_notrans',quarter="2quarter",contig_file=contig_file)

        my_kmeans(np.log(X_t), namelist, bin_number, marker1_2quarter_seed_num, bacar_marker_2quarter_seed_num, length_weight, output,
                  prefix='X_t_logtrans',quarter="2quarter",contig_file=contig_file)

        my_kmeans(X_t, namelist, bin_number, marker1_3quarter_seed_num, bacar_marker_3quarter_seed_num, length_weight, output,
                  prefix='X_t_notrans',quarter="3quarter",contig_file=contig_file)

        my_kmeans(np.log(X_t), namelist, bin_number, marker1_3quarter_seed_num, bacar_marker_3quarter_seed_num, length_weight, output,
                  prefix='X_t_logtrans',quarter="3quarter",contig_file=contig_file)


        if len(X_cov[0]) > 5:
            my_kmeans(X_cov, namelist, bin_number, marker1_1quarter_seed_num, bacar_marker_1quarter_seed_num,
                      length_weight, output,
                      prefix='X_cov_notrans', quarter="1quarter", contig_file=contig_file)

            my_kmeans(np.log(X_cov), namelist, bin_number, marker1_1quarter_seed_num, bacar_marker_1quarter_seed_num,
                      length_weight, output,
                      prefix='X_cov_logtrans', quarter="1quarter", contig_file=contig_file)

            my_kmeans(X_cov, namelist, bin_number, marker1_2quarter_seed_num, bacar_marker_2quarter_seed_num,
                      length_weight,
                      output,
                      prefix='X_cov_notrans', quarter="2quarter", contig_file=contig_file)

            my_kmeans(np.log(X_cov), namelist, bin_number, marker1_2quarter_seed_num, bacar_marker_2quarter_seed_num,
                      length_weight, output,
                      prefix='X_cov_logtrans', quarter="2quarter", contig_file=contig_file)

            my_kmeans(X_cov, namelist, bin_number, marker1_3quarter_seed_num, bacar_marker_3quarter_seed_num,
                      length_weight,
                      output,
                      prefix='X_cov_notrans', quarter="3quarter", contig_file=contig_file)

            my_kmeans(np.log(X_cov), namelist, bin_number, marker1_3quarter_seed_num, bacar_marker_3quarter_seed_num,
                      length_weight, output,
                      prefix='X_cov_logtrans', quarter="3quarter", contig_file=contig_file)


    # add post process for
    kmeans_initial_output = intermediate_result_dir + '/kmeans_length_weight_X_t_notrans_result.tsv'
    kmeans_initial_output_dir = os.path.dirname(
        args.output) + '/kmeans_length_weight_X_t_notrans_result'
    checkm_post_precess(kmeans_initial_output, kmeans_initial_output_dir, prefix='notrans')

    kmeans_initial_output = intermediate_result_dir + '/kmeans_length_weight_X_t_logtrans_result.tsv'
    kmeans_initial_output_dir = os.path.dirname(
        args.output) + '/kmeans_length_weight_X_t_logtrans_result'
    checkm_post_precess(kmeans_initial_output, kmeans_initial_output_dir, prefix='logtrans')

    # add post process for X_cov output
    if len(X_cov[0])>=6:
        kmeans_initial_output = intermediate_result_dir + '/kmeans_length_weight_X_cov_notrans_result.tsv'
        kmeans_initial_output_dir = os.path.dirname(
            args.output) + '/kmeans_length_weight_X_cov_notrans_result'
        checkm_post_precess(kmeans_initial_output, kmeans_initial_output_dir, prefix='notrans_X_cov')

        kmeans_initial_output = intermediate_result_dir + '/kmeans_length_weight_X_cov_logtrans_result.tsv'
        kmeans_initial_output_dir = os.path.dirname(
            args.output) + '/kmeans_length_weight_X_cov_logtrans_result'
        checkm_post_precess(kmeans_initial_output, kmeans_initial_output_dir, prefix='logtrans_X_cov')


    run_dastool(intermediate_result_dir, output, dir_name="/das_tool_output_all",
                search_engine=args.search_engine_for_dastool,partial_flag=None)
    run_dastool(intermediate_result_dir, output, dir_name="/das_tool_output_notrans",
                search_engine=args.search_engine_for_dastool,partial_flag="notrans")
    run_dastool(intermediate_result_dir, output, dir_name="/das_tool_output_logtrans",
                search_engine=args.search_engine_for_dastool,partial_flag="logtrans")

    # dir_name = "/das_tool_output"
    dir_name = "/das_tool_output_all"
    das_tool_all_result_path = os.path.dirname(output) + '/das_tool_output_all/das_tool_goodbins_DASTool_scaffolds2bin.txt'
    das_tool_notrans_result_path = os.path.dirname(output) + '/das_tool_output_notrans/das_tool_goodbins_DASTool_scaffolds2bin.txt'
    das_tool_logtrans_result_path = os.path.dirname(output) + '/das_tool_output_logtrans/das_tool_goodbins_DASTool_scaffolds2bin.txt'

    unique_binnumber_all,clustered_contig_id_all,result_id_all=get_unique_bin_number(das_tool_all_result_path)
    unique_binnumber_notrans,clustered_contig_id_notrans, result_id_notrans = get_unique_bin_number(das_tool_notrans_result_path)
    unique_binnumber_logtrans, clustered_contig_id_logtrans, result_id_logtrans = get_unique_bin_number(das_tool_logtrans_result_path)

    unique_binnumber_part=[unique_binnumber_notrans,unique_binnumber_logtrans]
    max_idx= unique_binnumber_part.index(max(unique_binnumber_part))
    logger.info("max unique_binnumber for dastool of three feature transforms : " + str(max(unique_binnumber_part)))
    logger.info("unique_binnumber for dastool of all the transforms : " + str(unique_binnumber_all))
    if max_idx==0:
        max_tran='notrans'
    if max_idx == 1:
        max_tran = 'logtrans'
    partial_result = intermediate_result_dir + '/partial_seed_kmeans_bacar_marker_seed_length_weight_3quarter_X_t_'+ max_tran+'_result.tsv'

    logger.info(max_tran+" has max unique_binnumber for dastool of three feature transform." )

    clustered_contig_id = clustered_contig_id_all
    result_id = result_id_all
    das_tool_output= das_tool_all_result_path
    logger.info("Das_tool result for all is saved as the result file.")

    f = open(args.output+'.2.tsv', 'w')
    for Idx in range(len(result_id)):
        f.write(clustered_contig_id[Idx] + "\t" + result_id[Idx] + "\n")
    f.close()

    logger.info("Run postprcessing for saved das_tool result.")


    output = args.output+'.2.tsv'+'.remained_after_dastool.tsv'
    fo = open(output, 'w')
    with open(partial_result) as f:
        for line in f:
            contig_id = line.rstrip('\n').split('\t')[0]
            if contig_id not in clustered_contig_id:
                fo.write(contig_id + "\t" + str(line.rstrip('\n').split('\t')[1]) + "\n")
    fo.close()

    output = args.output+'.2.tsv'+'.add_remained_after_dastool.tsv'
    fo = open(output, 'w')
    for Idx in range(len(result_id)):
        fo.write(clustered_contig_id[Idx] + "\t" + result_id[Idx] + "\n")
    with open(partial_result) as f:
        for line in f:
            contig_id = line.rstrip('\n').split('\t')[0]
            if contig_id not in clustered_contig_id:
                fo.write(contig_id + "\t" + str(line.rstrip('\n').split('\t')[1]) + "\n")
    fo.close()

