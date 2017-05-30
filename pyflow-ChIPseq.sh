#! /bin/bash


# invoke Snakemake in cluster mode with custom wrapper scripts

snakemake --rerun-incomplete -j 1000 --jobscript ./jobscript.sh \
		  --timestamp --latency-wait 120 --immediate-submit \
		  --cluster-config cluster.json --cluster './bsub_cluster.py {dependencies}' \
		  "$@"

