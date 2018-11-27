#! /bin/bash


# invoke Snakemake in cluster mode with custom wrapper scripts

#snakemake --rerun-incomplete -j 1000 -k \
#		  --jobscript ./jobscript.sh \
#		  --latency-wait 120 \
#		  --cluster-config cluster.json --cluster './sbatch_cluster.py' \
#		  "$@"

#! /bin/bash

## make a folder for bsub_log if not exist
folder="sbatch_log"


if [ ! -d "$folder" ]; then
    mkdir -p "$folder"
fi

snakemake --jobs 99999 \
    --latency-wait 240 \
    --cluster-config cluster.json \
    --cluster "sbatch --partition {cluster.p} --mem {cluster.mem} -N {cluster.N} -n {cluster.n} --time {cluster.time} --job-name {cluster.name} --output {cluster.output} --error {cluster.error}" \
    "$@"
	 

