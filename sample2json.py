#!/usr/bin/env python3


import json
import os
import csv
import re
from os.path import join
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--fastq_dir", help="Required. the FULL path to the fastq folder")
parser.add_argument("--meta", help="Required. the FULL path to the tab delimited meta file")
args = parser.parse_args()

assert args.fastq_dir is not None, "please provide the path to the fastq folder"
assert args.meta is not None, "please provide the path to the meta file"


## collect all the fastq.gz full path in to a list
fastq_paths = []

for root, dirs, files in os.walk(args.fastq_dir):
    for file in files:
        if file.endswith("fastq.gz"):
            full_path = join(root, file)
            fastq_paths.append(full_path)


FILES = defaultdict(lambda: defaultdict(list))

with open(args.meta, "r") as f:
    reader = csv.reader(f, delimiter = "\t")
    # skip the header
    header = next(reader)
    for row in reader:
    	sample_name = row[0].strip()
    	fastq_name = row[1].strip()
    	sample_type = row[2].strip()
    	## now just assume the file name in the metafile contained in the fastq file path
    	fastq_full_path = [x for x in fastq_paths if fastq_name in x]
    	if fastq_full_path:
    		FILES[sample_name][sample_type].extend(fastq_full_path)
    	else:
    		print("sample {sample_name} missing {sample_type} {fastq_name} fastq files".format(sample_name = sample_name, sample_type = sample_type, fastq_name = fastq_name))


print()
sample_num = len(FILES.keys())
print ("total {} unique samples will be processed".format(sample_num))
print ("------------------------------------------")
for sample_name in sorted(FILES.keys()):
	for sample_type in FILES[sample_name]:
		fastq_file = "".join(FILES[sample_name][sample_type])
		print("sample {sample_name}'s {sample_type} fastq path is {fastq_file}".format(sample_name = sample_name, sample_type = sample_type, fastq_file = fastq_file))
print ("------------------------------------------")
for sample in FILES.keys():
	print ("{sample} has {n} marks".format(sample = sample, n = len(FILES[sample])))
print ("------------------------------------------")
print("check the samples.json file for fastqs belong to each sample")
print()

js = json.dumps(FILES, indent = 4, sort_keys=True)
open('samples.json', 'w').writelines(js)
