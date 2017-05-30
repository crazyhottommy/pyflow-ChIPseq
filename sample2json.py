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


## read in the meta file pairing IP and input
## the first column (IP) and the second column (Input) of the meta.txt file stroes the
## name information.

## a list of lists to hold the IP, Input pairs 
pair = []
with open(args.meta, "r") as f:
    reader = csv.reader(f, delimiter = "\t")
    # skip the header
    header = next(reader)
    for row in reader:
        pair.append([row[0], row[1]])


## default dictionary is quite useful!
FILES = defaultdict(lambda: defaultdict(list))

## build the dictionary with full path for each fastq.gz file
for root, dirs, files in os.walk(args.fastq_dir):
    for file in files:
        if file.endswith("fastq.gz"):
            full_path = join(root, file)
            find_sample = lambda x : x in file
            ## python3 returns a filter object(iterable) not a list anymore
            ## below is the same as [ x[0] for x in pair if x[0] in file]
            ## but when the filter requirments are more complex, a function of find_sample can
            ## be more flexiable.
            IP = list(filter(find_sample, [x[0] for x in pair]))
            Input = list(filter(find_sample, [x[1] for x in pair]))
            if IP:
                FILES[IP[0]]["IP"].append(full_path)
            if Input:
                ## multiple IPs could use the same Input
                IPs = [x[0] for x in pair if x[1] == Input[0]]
                for IP in IPs:
                    FILES[IP]["Input"].append(full_path)
            else:
                print("sample {file} is not described in the meta file".format(file = file))

print()
print ("total {} unique IP Input pairs will be processed".format(len(FILES.keys())))
print ("------------------------------------------")
for sample in FILES.keys():
	Input = [x[1] for x in pair if x[0] == sample][0]
	print ("IP {sample} pairs with Input {input}".format(sample = sample, input = Input))
print ("------------------------------------------")
print("check the samples.json file for fastqs belong to each sample")
print()

js = json.dumps(FILES, indent = 4, sort_keys=True)
open('samples.json', 'w').writelines(js)
