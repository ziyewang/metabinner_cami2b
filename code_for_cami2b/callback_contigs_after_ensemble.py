# -*- coding: utf-8 -*-
# @Author  : ZWang
# @FileName: Metabinner.py
# scikit-learn == 0.20.4
# python 2.7

import numpy as np
import pandas as pd
import sys
import os
import time
import functools

import mimetypes
import gzip

from Bio import SeqIO
import logging
import argparse

import sklearn.cluster.k_means_ as kmeans
import scipy.sparse as sp
from sklearn.cluster.k_means_ import euclidean_distances, stable_cumsum, KMeans, check_random_state, row_norms
from sklearn.metrics import pairwise_distances

from scipy.sparse import csr_matrix
from collections import defaultdict
from scipy.sparse import csc_matrix
from scipy import linalg
from sklearn.preprocessing import normalize

from sklearn.cluster import SpectralClustering
from sklearn.manifold import spectral_embedding

import csv
import mimetypes
import os
import re
import shutil
import subprocess

from collections import defaultdict, Counter


logger = logging.getLogger('Metabinner 2.0 python 3')

logger.setLevel(logging.INFO)

# logging
formatter = logging.Formatter('%(asctime)s - %(message)s')

console_hdr = logging.StreamHandler()
console_hdr.setFormatter(formatter)

logger.addHandler(console_hdr)

def arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('--contig_file', type=str, help=("The contigs file."))
    parser.add_argument('--intermediate_result', type=str, help=("The intermediate_result (bacer 3quarter) file."))
    parser.add_argument('--das_tool_result', type=str, help=("The das_tool_result."))
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



def read_fasta_file(fasta_file):
    with open(fasta_file, 'r') as read_handler:
        for line in read_handler:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                yield line[1:]



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
    logger.info("Binscore threshold for Das_tool:\t" + str(args.binscore))
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

    intermediate_result=args.intermediate_result
    das_tool_result=args.das_tool_result

    unique_binnumber,clustered_contig_id,result_id=get_unique_bin_number(das_tool_result)


    logger.info("Run postprcessing for saved das_tool result.")

    output = das_tool_result+'.add_remained_after_dastool.tsv'
    fo = open(output, 'w')
    for Idx in range(len(result_id)):
        fo.write(clustered_contig_id[Idx] + "\t" + result_id[Idx] + "\n")
    with open(intermediate_result) as f:
        for line in f:
            contig_id = line.rstrip('\n').split('\t')[0]
            if contig_id not in clustered_contig_id:
                fo.write(contig_id + "\tcallbackBin_" + str(line.rstrip('\n').split('\t')[1]) + "\n")
    fo.close()
