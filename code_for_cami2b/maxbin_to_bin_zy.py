#!/usr/bin/env python

import argparse
import os

def read_fasta_file(fasta_file):
    with open(fasta_file, 'r') as read_handler:
        for line in read_handler:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                yield line[1:]


def convert(paths, output_file):
    files = os.listdir(paths)
    fasta_files =[]
    for file in files :
        if file.endswith(('.fasta', '.fa','.fna','.bin','_0.bin','_1.bin')):
            fasta_files.append(file)
    with open(output_file, 'w') as write_handler:
        for bin_id, fasta_file in enumerate(fasta_files):
            for sequence_id in read_fasta_file(paths+'/'+fasta_file):
                write_handler.write("%s,%s\n" % (sequence_id, bin_id))


def main():
    parser = argparse.ArgumentParser(description="Convert bins file to one result file")
    parser.add_argument("--paths", help="FASTA files path")
    parser.add_argument("-o", "--output_file", required=False, help="Output file")
    args = parser.parse_args()
    convert(args.paths, args.output_file)


if __name__ == "__main__":
    main()
