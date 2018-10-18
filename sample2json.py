#!/usr/bin/env python3


import json
import os
import csv
import re
from os.path import join
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("fastq_dir", help="Required. the FULL path to the fastq folder")
parser.add_argument("meta", help="Required. the FULL path to the tab delimited meta file")
args = parser.parse_args()

assert args.fastq_dir is not None, "please provide the path to the fastq folder"
assert args.meta is not None, "please provide the path to the meta file"


## collect all the fastq.gz full path in to a list
fastq_paths = []

for root, dirs, files in os.walk(args.fastq_dir):
    for file in files:
        if file.endswith("fq.gz") or file.endswith("fastq.gz"):
            full_path = join(root, file)
            fastq_paths.append(full_path)


FILES = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

with open(args.meta, "r") as f:
    reader = csv.reader(f, delimiter = "\t")
    # skip the header
    header = next(reader)
    for row in reader:
        sample_name = row[0].strip()
        fastq_name = row[1].strip()
        factor = row[2].strip()
        # forward or reverse, R1 or R2
        reads = row[3].strip()
        ## now just assume the file name in the metafile contained in the fastq file path
        fastq_full_path = [x for x in fastq_paths if fastq_name in x]
        if fastq_full_path:
            FILES[sample_name][factor][reads].extend(fastq_full_path)
        else:
            print("sample {sample_name} missing {factor} {reads} {fastq_name} fastq files".format(sample_name = sample_name, factor = factor, reads = reads, fastq_name = fastq_name))


print()
sample_num = len(FILES.keys())
print ("total {} unique samples will be processed".format(sample_num))
print ("------------------------------------------")
for sample_name in sorted(FILES.keys()):
    factors = sorted(FILES[sample_name].keys())
    print ("{sample_name} has {n} marks: {factors}".format(sample_name = sample_name, n = len(FILES[sample_name]), factors= " ".join(factors)))
print ("------------------------------------------")
print("check the samples.json file for fastqs belong to each sample")
print()

js = json.dumps(FILES, indent = 4, sort_keys=True)
open('samples.json', 'w').writelines(js)

