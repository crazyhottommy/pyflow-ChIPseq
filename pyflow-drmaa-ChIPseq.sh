#! /bin/bash

## use drmma for better job control, do not use --immediate-submit, start snakemake in a screen session

## make a folder for bsub_log if not exist
folder="bsub_log"


if [ ! -d "$folder" ]; then
    mkdir -p "$folder"
fi


snakemake --rerun-incomplete -j 1000 --jobscript ./jobscript.sh \
		  --timestamp --latency-wait 120 --jobname {rulename}.{jobid} \
		  --cluster-config cluster.json \
		  --drmaa " -n {cluster.cpu} -W {cluster.time} -u {cluster.email} -q {cluster.queue} -M {cluster.MaxMem} -R rusage[mem={cluster.MaxMem}] -o bsub_log/ -e bsub_log/" \
		  "$@"
