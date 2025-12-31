"""
pyflow-ChIPseq: A Snakemake pipeline for ChIP-seq analysis
Modernized version (2025) with Python 3, MACS3, and updated tools
"""

import json
import csv
import os

configfile: "config.yaml"

# localrules will run locally rather than submitting to cluster
localrules: all

# Load sample configuration
FILES = json.load(open(config['SAMPLES_JSON']))

SAMPLES = sorted(FILES.keys())

## list all samples by sample_name and sample_type
MARK_SAMPLES = []
for sample in SAMPLES:
    for sample_type in FILES[sample].keys():
        MARK_SAMPLES.append(sample + "_" + sample_type)


# which sample_type is used as control for calling peaks: e.g. Input, IgG...
CONTROL = config["control"]
CONTROLS = [sample for sample in MARK_SAMPLES if CONTROL in sample]
CASES = [sample for sample in MARK_SAMPLES if CONTROL not in sample]

## multiple samples may use the same control input/IgG files
CONTROLS_UNIQUE = list(set(CONTROLS))


## list BAM files
CONTROL_BAM = expand("03aln/{sample}.sorted.bam", sample=CONTROLS_UNIQUE)
CASE_BAM = expand("03aln/{sample}.sorted.bam", sample=CASES)

## peaks and bigwigs
ALL_PEAKS = []
ALL_inputSubtract_BIGWIG = []
ALL_SUPER = []

for case in CASES:
    sample = "_".join(case.split("_")[0:-1])
    control = sample + "_" + CONTROL
    if control in CONTROLS:
        ALL_PEAKS.append("08peak_macs3/{}_vs_{}_narrow_peaks.narrowPeak".format(case, control))
        ALL_PEAKS.append("09peak_macs3/{}_vs_{}_broad_peaks.broadPeak".format(case, control))
        ALL_inputSubtract_BIGWIG.append("06bigwig_inputSubtract/{}_subtract_{}.bw".format(case, control))
        ALL_SUPER.append("11superEnhancer/{}_vs_{}-super/".format(case, control))

ALL_SAMPLES = CASES + CONTROLS_UNIQUE
ALL_BAM     = CONTROL_BAM + CASE_BAM
ALL_DOWNSAMPLE_BAM = expand("04aln_downsample/{sample}-downsample.sorted.bam", sample = ALL_SAMPLES)
ALL_INDEX = expand("03aln/{sample}.sorted.bam.bai", sample = ALL_SAMPLES)
ALL_DOWNSAMPLE_INDEX = expand("04aln_downsample/{sample}-downsample.sorted.bam.bai", sample = ALL_SAMPLES)
ALL_FLAGSTAT = expand("03aln/{sample}.sorted.bam.flagstat", sample = ALL_SAMPLES)
ALL_BIGWIG = expand("07bigwig/{sample}.bw", sample = ALL_SAMPLES)
ALL_QC = ["10multiQC/multiQC_log.html"]

# Conditional outputs based on config
if config.get("run_phantompeakqual", True):
    ALL_PHATOM = expand("05phantompeakqual/{sample}_phantom.txt", sample = ALL_SAMPLES)
else:
    ALL_PHATOM = []

TARGETS = []
TARGETS.extend(ALL_BAM)
TARGETS.extend(ALL_DOWNSAMPLE_BAM)
TARGETS.extend(ALL_INDEX)
TARGETS.extend(ALL_PHATOM)
TARGETS.extend(ALL_PEAKS)
TARGETS.extend(ALL_BIGWIG)
TARGETS.extend(ALL_inputSubtract_BIGWIG)
TARGETS.extend(ALL_FLAGSTAT)
TARGETS.extend(ALL_QC)
TARGETS.extend(ALL_SUPER)

## sometimes if if you have TF ChIP-seq data, do not include it to chromHMM, or you want
## only a subset of the histone marks be included in the chromHMM call

if config["chromHMM"]:
    HISTONE_INCLUDED = config["histone_for_chromHMM"].split(" ")
    HISTONE_CASES = [sample for sample in MARK_SAMPLES if sample.split("_")[-1] in HISTONE_INCLUDED ]
    ALL_BED = expand("12bed/{sample}.bed", sample = HISTONE_CASES + CONTROLS)
    CHROMHMM = ["13chromHMM/MYOUTPUT", "13chromHMM/binarizedData"]
    CHROMHMM_TABLE = ["12bed/cellmarkfiletable.txt"]
    TARGETS.extend(ALL_BED)
    TARGETS.extend(CHROMHMM)
    TARGETS.extend(CHROMHMM_TABLE)


localrules: all
rule all:
    input: TARGETS


## get a list of fastq.gz files for the same mark, same sample
def get_input_files(wildcards):
    sample_name = "_".join(wildcards.sample.split("_")[0:-1])
    mark = wildcards.sample.split("_")[-1]
    if not config["paired_end"]:
        return FILES[sample_name][mark]['R1']
    ## start with bam files
    if not config["from_fastq"]:
        return FILES[sample_name][mark]

def get_R1_files(wildcards):
    sample_name = "_".join(wildcards.sample.split("_")[0:-1])
    mark = wildcards.sample.split("_")[-1]
    if config["from_fastq"]:
        if config["paired_end"]:
            return  FILES[sample_name][mark]['R1']

def get_R2_files(wildcards):
    sample_name = "_".join(wildcards.sample.split("_")[0:-1])
    mark = wildcards.sample.split("_")[-1]
    if config["from_fastq"]:
        if config["paired_end"]:
            return  FILES[sample_name][mark]['R2']

## when there are multiple fastq.gz for the same sample, same mark, merge them.
if config["paired_end"]:
    rule merge_fastqs:
        input:
            r1 = get_R1_files,
            r2 = get_R2_files
        output:
            temp("01seq/{sample}_R1.fastq.gz"),
            temp("01seq/{sample}_R2.fastq.gz")
        log: "00log/{sample}_merge_fastq.log"
        params:
            jobname = "{sample}"
        threads: 1
        resources:
            mem_mb=4000,
            runtime=120
        message: "processing fastqs for {wildcards.sample}"
        run:
            # R1 handling: symlink if single file, merge if multiple
            if len(input.r1) == 1:
                shell("ln -sf $(realpath {input.r1}) {output[0]} 2> {log}")
            else:
                shell("gunzip -c {input.r1} | gzip > {output[0]} 2> {log}")

            # R2 handling: symlink if single file, merge if multiple
            if len(input.r2) == 1:
                shell("ln -sf $(realpath {input.r2}) {output[1]} 2>> {log}")
            else:
                shell("gunzip -c {input.r2} | gzip > {output[1]} 2>> {log}")
else:
    rule merge_fastqs:
        input: get_input_files
        output: temp("01seq/{sample}.fastq.gz")
        log: "00log/{sample}_unzip"
        threads: 1
        params: jobname = "{sample}"
        resources:
            mem_mb=4000,
            runtime=120
        message: "processing fastqs for {wildcards.sample}"
        run:
            # Symlink if single file, merge if multiple
            if len(input) == 1:
                shell("ln -sf $(realpath {input}) {output} 2> {log}")
            else:
                shell("gunzip -c {input} | gzip > {output} 2> {log}")




if config["paired_end"]:
    rule fastqc:
        input:  "01seq/{sample}_R1.fastq.gz", "01seq/{sample}_R2.fastq.gz"
        output: "02fqc/{sample}_R1_fastqc.zip", "02fqc/{sample}_R2_fastqc.zip"
        log:    "00log/{sample}_fastqc"
        threads: 1
        params : jobname = "{sample}"
        conda: "envs/qc.yaml"
        resources:
            mem_mb=4000,
            runtime=60
        message: "fastqc {input}: {threads} threads"
        shell:
            """
            fastqc -o 02fqc -f fastq --noextract {input[0]} {input[1]} 2> {log}
            """
else:
    rule fastqc:
        input:  "01seq/{sample}.fastq.gz"
        output: "02fqc/{sample}_fastqc.zip", "02fqc/{sample}_fastqc.html"
        log:    "00log/{sample}_fastqc"
        threads: 1
        params : jobname = "{sample}"
        conda: "envs/qc.yaml"
        resources:
            mem_mb=4000,
            runtime=60
        message: "fastqc {input}: {threads}"
        shell:
            """
            fastqc -o 02fqc -f fastq --noextract {input} 2> {log}
            """

# get the duplicates marked sorted bam, remove unmapped reads by samtools view -F 4 and dupliated reads by samblaster -r
# samblaster should run before samtools sort

## align from fastq files
## deal with paired end and single end data
## see https://bitbucket.org/snakemake/snakemake/pull-requests/148/examples-for-paired-read-data-and/diff
## some hints:
## https://bitbucket.org/snakemake/snakemake/issues/37/add-complex-conditional-file-dependency
## conditional rules?
## https://groups.google.com/forum/#!msg/snakemake/qX7RfXDTDe4/cKZBfc_PAAAJ


## when feeding @RG to bwa mem,  literal tab will cause the resulting bam file header violating the
## SAM specification at @PG line. https://github.com/lh3/bwa/issues/83. This will interfere with IGV and picard
## to process the bam files.
## and https://github.com/Duke-GCB/bespin-cwl/issues/6
## BWA after version BWA-0.7.16a (r1181) with verbose >=1 will  this problem.
## now feed bwa with read group by escaping the tab
## Before: -R "@RG\tID:1\tLB:LIBRARY\tPL:illumina\tSM:sample1\tPU:AB1234"
## Now: -R "@RG\\tID:1\\tLB:LIBRARY\\tPL:illumina\\tSM:sample1\\tPU:AB1234"


## ============================================================================
## PAIRED-END ALIGNMENT (Split into 3 steps for easier debugging)
## ============================================================================

if config["from_fastq"] and config["paired_end"]:
    # Step 1: BWA alignment only (produces unsorted BAM)
    rule bwa_align_pe:
        input:
            r1 = "01seq/{sample}_R1.fastq.gz",
            r2 = "01seq/{sample}_R2.fastq.gz"
        output:
            bam = temp("03aln/{sample}.unsorted.bam")
        threads: config.get("max_threads_align", 8)
        params:
            jobname = "{sample}",
            rg = "@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA"
        conda: "envs/alignment.yaml"
        resources:
            mem_mb=16000,
            runtime=360
        log: "00log/{sample}.bwa_align"
        message: "BWA alignment for {wildcards.sample}: {threads} threads"
        run:
            if config["long_reads"]:
                # Long reads: use bwa mem
                shell(
                    r"""
                    bwa mem -t {threads} -M -v 1 -R '{params.rg}' {config[ref_fa]} \
                        {input.r1} {input.r2} 2> {log} \
                        | samtools view -@ {threads} -b -o {output.bam}
                    """)
            else:
                # Short reads: use bwa aln + sampe
                shell(
                    r"""
                    mkdir -p 01aln
                    bwa aln -t {threads} {config[ref_fa]} {input.r1} 2> {log} > 01aln/{wildcards.sample}_R1.sai
                    bwa aln -t {threads} {config[ref_fa]} {input.r2} 2>> {log} > 01aln/{wildcards.sample}_R2.sai
                    bwa sampe -r '{params.rg}' {config[ref_fa]} \
                        01aln/{wildcards.sample}_R1.sai 01aln/{wildcards.sample}_R2.sai \
                        {input.r1} {input.r2} 2>> {log} \
                        | samtools view -@ {threads} -b -o {output.bam}
                    rm 01aln/{wildcards.sample}_R1.sai 01aln/{wildcards.sample}_R2.sai
                    """)

    # Step 2: Remove duplicates, filter, and sort (PE)
    rule process_bam_pe:
        input:
            bam = "03aln/{sample}.unsorted.bam"
        output:
            bam = "03aln/{sample}.sorted.bam"
        threads: config.get("max_threads_align", 8)
        params:
            jobname = "{sample}"
        conda: "envs/alignment.yaml"
        resources:
            mem_mb=16000,
            runtime=360
        log:
            markdup = "00log/{sample}.markdup",
            sort = "00log/{sample}.sort"
        message: "Processing PE BAM for {wildcards.sample}: remove duplicates, filter, sort"
        shell:
            r"""
            samtools view -h {input.bam} \
                | samblaster -r 2> {log.markdup} \
                | samtools view -@ {threads} -b -F 4 - \
                | samtools sort -m 2G -@ {threads} -T {output.bam}.tmp -o {output.bam} 2> {log.sort}
            """

    # Step 3: Index BAM (PE)
    rule index_bam_pe:
        input:
            bam = "03aln/{sample}.sorted.bam"
        output:
            bai = "03aln/{sample}.sorted.bam.bai"
        threads: 2
        params:
            jobname = "{sample}"
        conda: "envs/alignment.yaml"
        resources:
            mem_mb=4000,
            runtime=60
        log: "00log/{sample}.index"
        message: "Indexing BAM for {wildcards.sample}"
        shell:
            """
            samtools index -@ {threads} {input.bam} 2> {log}
            """


## ============================================================================
## SINGLE-END ALIGNMENT (Split into 3 steps for easier debugging)
## ============================================================================

if config["from_fastq"] and not config["paired_end"]:
    # Step 1: BWA alignment only (produces unsorted BAM)
    rule bwa_align_se:
        input:
            fastq = "01seq/{sample}.fastq.gz"
        output:
            bam = temp("03aln/{sample}.unsorted.bam")
        threads: config.get("max_threads_align", 8)
        params:
            jobname = "{sample}",
            rg = "@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA"
        conda: "envs/alignment.yaml"
        resources:
            mem_mb=16000,
            runtime=360
        log: "00log/{sample}.bwa_align"
        message: "BWA alignment for {wildcards.sample}: {threads} threads"
        run:
            if config["long_reads"]:
                # Long reads: use bwa mem
                shell(
                    r"""
                    bwa mem -t {threads} -M -v 1 -R '{params.rg}' {config[ref_fa]} \
                        {input.fastq} 2> {log} \
                        | samtools view -@ {threads} -b -o {output.bam}
                    """)
            else:
                # Short reads: use bwa aln + samse
                shell(
                    r"""
                    mkdir -p 01aln
                    bwa aln -t {threads} {config[ref_fa]} {input.fastq} 2> {log} > 01aln/{wildcards.sample}.sai
                    bwa samse -r '{params.rg}' {config[ref_fa]} \
                        01aln/{wildcards.sample}.sai {input.fastq} 2>> {log} \
                        | samtools view -@ {threads} -b -o {output.bam}
                    rm 01aln/{wildcards.sample}.sai
                    """)

    # Step 2: Remove duplicates, filter, and sort (SE)
    rule process_bam_se:
        input:
            bam = "03aln/{sample}.unsorted.bam"
        output:
            bam = "03aln/{sample}.sorted.bam"
        threads: config.get("max_threads_align", 8)
        params:
            jobname = "{sample}"
        conda: "envs/alignment.yaml"
        resources:
            mem_mb=16000,
            runtime=360
        log:
            markdup = "00log/{sample}.markdup",
            sort = "00log/{sample}.sort"
        message: "Processing SE BAM for {wildcards.sample}: remove duplicates, filter, sort"
        shell:
            r"""
            samtools view -h {input.bam} \
                | samblaster -r --ignoreUnmated 2> {log.markdup} \
                | samtools view -@ {threads} -b -F 4 - \
                | samtools sort -m 2G -@ {threads} -T {output.bam}.tmp -o {output.bam} 2> {log.sort}
            """

    # Step 3: Index BAM (SE)
    rule index_bam_se:
        input:
            bam = "03aln/{sample}.sorted.bam"
        output:
            bai = "03aln/{sample}.sorted.bam.bai"
        threads: 2
        params:
            jobname = "{sample}"
        conda: "envs/alignment.yaml"
        resources:
            mem_mb=4000,
            runtime=60
        log: "00log/{sample}.index"
        message: "Indexing BAM for {wildcards.sample}"
        shell:
            """
            samtools index -@ {threads} {input.bam} 2> {log}
            """

## ============================================================================
## SHARED RULES (for both PE and SE)
## ============================================================================

# Check number of reads mapped by samtools flagstat, the output will be used for downsampling
rule flagstat_bam:
    input:  "03aln/{sample}.sorted.bam"
    output: "03aln/{sample}.sorted.bam.flagstat"
    log:    "00log/{sample}.flagstat_bam"
    threads: 1
    params: jobname = "{sample}"
    conda: "envs/alignment.yaml"
    resources:
        mem_mb=2000,
        runtime=30
    message: "flagstat_bam {input}: {threads} threads"
    shell:
        """
        samtools flagstat --threads {threads} {input} > {output} 2> {log}
        """

rule phantom_peak_qual:
    input:
        bam = "03aln/{sample}.sorted.bam",
        bai = "03aln/{sample}.sorted.bam.bai"
    output: "05phantompeakqual/{sample}_phantom.txt"
    log: "00log/{sample}_phantompeakqual.log"
    threads: config.get("max_threads_other", 4)
    params: jobname = "{sample}"
    conda: "envs/phantompeakqualtools.yaml"
    resources:
        mem_mb=8000,
        runtime=120
    message: "phantompeakqual for {input}"
    shell:
        """
        # Download run_spp.R if not present
        if [ ! -f run_spp.R ]; then
            curl -O https://raw.githubusercontent.com/kundajelab/phantompeakqualtools/master/run_spp.R
        fi

        # Install spp R package from GitHub if not already installed
        if ! Rscript -e "library(spp)" 2>/dev/null; then
            echo "Installing spp from GitHub..." 2>> {log}
            # Download from master branch (most reliable)
            curl -L -o spp.tar.gz https://github.com/hms-dbmi/spp/archive/master.tar.gz 2>> {log}
            # Extract and install from source
            tar -xzf spp.tar.gz 2>> {log}
            R CMD INSTALL spp-master 2>> {log}
            rm -rf spp.tar.gz spp-master
        fi

        # Run phantompeakqualtools
        Rscript run_spp.R -c={input.bam} -savp -rf -p={threads} -odir=05phantompeakqual -out={output} -tmpdir=05phantompeakqual 2>> {log}

        # Clean up temporary files created by run_spp.R (can be several GB)
        rm -f 05phantompeakqual/{wildcards.sample}.sorted.bam*[0-9a-f][0-9a-f][0-9a-f]*
        rm -f 05phantompeakqual/{wildcards.sample}.sorted.tagAlign*[0-9a-f][0-9a-f][0-9a-f]*
        """

rule down_sample:
    input:
        bam = "03aln/{sample}.sorted.bam",
        bai = "03aln/{sample}.sorted.bam.bai",
        flagstat = "03aln/{sample}.sorted.bam.flagstat"
    output:
        bam = "04aln_downsample/{sample}-downsample.sorted.bam",
        bai = "04aln_downsample/{sample}-downsample.sorted.bam.bai"
    log: "00log/{sample}_downsample.log"
    threads: config.get("max_threads_other", 4)
    params:
        jobname = "{sample}"
    conda: "envs/alignment.yaml"
    resources:
        mem_mb=16000,
        runtime=240
    message: "downsampling for {input}"
    run:
        import re
        with open (input.flagstat, "r") as f:
            # Read line 8 (index 7): primary mapped reads (excludes secondary/supplementary)
            line = f.readlines()[7]
            match_number = re.match(r'(\d+) \+', line)
            if not match_number:
                raise ValueError(f"Could not parse flagstat line: {line}")
            total_reads = int(match_number.group(1))
        if config["downsample"]:
            target_reads = config["target_reads"]
            if total_reads > target_reads:
                down_rate = target_reads/total_reads
            else:
                down_rate = 1
            shell("sambamba view -f bam -t {threads} --subsampling-seed=3 -s {rate} {inbam} | samtools sort -m 2G -@ {threads} -T {outbam}.tmp -o {outbam} 2> {log}".format(rate = down_rate, inbam = input.bam, outbam = output.bam, log = log, threads = threads))
            shell("samtools index -@ {threads} {outbam}".format(outbam = output.bam, threads = threads))
        else:
            shell("ln -sf $(pwd)/{inbam} {outbam}".format(inbam = input.bam, outbam = output.bam))
            shell("ln -sf $(pwd)/{inbai} {outbai}".format(inbai = input.bai, outbai = output.bai))

rule make_inputSubtract_bigwigs:
    input:
        control_bam = "04aln_downsample/{control}-downsample.sorted.bam",
        case_bam = "04aln_downsample/{case}-downsample.sorted.bam",
        control_bai = "04aln_downsample/{control}-downsample.sorted.bam.bai",
        case_bai = "04aln_downsample/{case}-downsample.sorted.bam.bai"
    output: "06bigwig_inputSubtract/{case}_subtract_{control}.bw"
    log: "00log/{case}_subtract_{control}.makebw"
    threads: config.get("max_threads_other", 4)
    params: jobname = "{case}_{control}"
    conda: "envs/bigwig.yaml"
    resources:
        mem_mb=16000,
        runtime=240
    message: "making input subtracted bigwig for {input}"
    shell:
        """
        bamCompare --bamfile1 {input.case_bam} --bamfile2 {input.control_bam} \
            --normalizeUsing RPKM --scaleFactorsMethod None --operation subtract --binSize 30 \
            --smoothLength 300 --numberOfProcessors {threads} \
            --extendReads 200 -o {output} 2> {log}
        """

rule make_bigwigs:
    input:
        bam = "04aln_downsample/{sample}-downsample.sorted.bam",
        bai = "04aln_downsample/{sample}-downsample.sorted.bam.bai"
    output: "07bigwig/{sample}.bw"
    log: "00log/{sample}.makebw"
    threads: config.get("max_threads_other", 4)
    params: jobname = "{sample}"
    conda: "envs/bigwig.yaml"
    resources:
        mem_mb=16000,
        runtime=240
    message: "making bigwig for {input}"
    shell:
        """
        bamCoverage -b {input.bam} --normalizeUsing RPKM --binSize 30 \
            --smoothLength 300 --numberOfProcessors {threads} \
            --extendReads 200 -o {output} 2> {log}
        """

rule call_peaks_macs3_narrow:
    input:
        control = "04aln_downsample/{control}-downsample.sorted.bam",
        case = "04aln_downsample/{case}-downsample.sorted.bam"
    output:
        peaks = "08peak_macs3/{case}_vs_{control}_narrow_peaks.narrowPeak",
        xls = "08peak_macs3/{case}_vs_{control}_narrow_peaks.xls"
    log: "00log/{case}_vs_{control}_call_peaks_macs3_narrow.log"
    threads: 1
    params:
        name = "{case}_vs_{control}_narrow",
        outdir = "08peak_macs3",
        jobname = "{case}"
    conda: "envs/peakcalling.yaml"
    resources:
        mem_mb=16000,
        runtime=240
    message: "call_peaks macs3 narrow {input}: {threads} threads"
    shell:
        """
        macs3 callpeak -t {input.case} \
            -c {input.control} --keep-dup all -f BAM -g {config[macs_g]} \
            --outdir {params.outdir} -n {params.name} \
            -q {config[macs_pvalue]} &> {log}
        """

rule call_peaks_macs3_broad:
    input:
        control = "04aln_downsample/{control}-downsample.sorted.bam",
        case = "04aln_downsample/{case}-downsample.sorted.bam"
    output:
        peaks = "09peak_macs3/{case}_vs_{control}_broad_peaks.broadPeak",
        xls = "09peak_macs3/{case}_vs_{control}_broad_peaks.xls"
    log: "00log/{case}_vs_{control}_call_peaks_macs3_broad.log"
    threads: 1
    params:
        name = "{case}_vs_{control}_broad",
        outdir = "09peak_macs3",
        jobname = "{case}"
    conda: "envs/peakcalling.yaml"
    resources:
        mem_mb=16000,
        runtime=240
    message: "call_peaks macs3 broad {input}: {threads} threads"
    shell:
        """
        macs3 callpeak -t {input.case} \
            -c {input.control} --keep-dup all -f BAM -g {config[macs2_g]} \
            --outdir {params.outdir} -n {params.name} \
            -q {config[macs2_pvalue]} --broad --broad-cutoff {config[macs2_pvalue_broad]} &> {log}
        """


if config["paired_end"]:
    rule multiQC:
        input :
            expand("00log/{sample}.bwa_align", sample = ALL_SAMPLES),
            expand("00log/{sample}.markdup", sample = ALL_SAMPLES),
            expand("03aln/{sample}.sorted.bam.flagstat", sample = ALL_SAMPLES),
            expand("02fqc/{sample}_{reads}_fastqc.zip", sample = ALL_SAMPLES, reads = ['R1', 'R2'])
        output: "10multiQC/multiQC_log.html"
        log: "00log/multiqc.log"
        threads: 1
        conda: "envs/qc.yaml"
        resources:
            mem_mb=4000,
            runtime=60
        message: "multiqc for all logs"
        shell:
            """
            multiqc 02fqc 03aln 00log -o 10multiQC -d -f -v -n multiQC_log 2> {log}
            """
else:
    rule multiQC:
        input :
            expand("00log/{sample}.bwa_align", sample = ALL_SAMPLES),
            expand("00log/{sample}.markdup", sample = ALL_SAMPLES),
            expand("03aln/{sample}.sorted.bam.flagstat", sample = ALL_SAMPLES),
            expand("02fqc/{sample}_fastqc.zip", sample = ALL_SAMPLES)
        output: "10multiQC/multiQC_log.html"
        log: "00log/multiqc.log"
        threads: 1
        conda: "envs/qc.yaml"
        resources:
            mem_mb=4000,
            runtime=60
        message: "multiqc for all logs"
        shell:
            """
            multiqc 02fqc 03aln 00log -o 10multiQC -d -f -v -n multiQC_log 2> {log}
            """

## ROSE2 for super-enhancer calling (Python 3 port)
## See: https://github.com/crazyhottommy/ROSE2
rule superEnhancer:
    input:
        control_bam = "04aln_downsample/{control}-downsample.sorted.bam",
        case_bam = "04aln_downsample/{case}-downsample.sorted.bam",
        control_bai = "04aln_downsample/{control}-downsample.sorted.bam.bai",
        case_bai = "04aln_downsample/{case}-downsample.sorted.bam.bai",
        peaks = "08peak_macs3/{case}_vs_{control}_narrow_peaks.narrowPeak"
    output: directory("11superEnhancer/{case}_vs_{control}-super/")
    log: "00log/{case}_vs_{control}_superEnhancer.log"
    threads: 4
    params:
        jobname = "{case}_vs_{control}"
    conda: "envs/rose.yaml"
    resources:
        mem_mb=16000,
        runtime=240
    message: "Calling super-enhancers with ROSE2 for {wildcards.case} vs {wildcards.control}"
    shell:
        """
        # ROSE2 - Python 3 version of ROSE for super-enhancer calling
        # See: https://github.com/crazyhottommy/ROSE2
        rose2 main \
            -i {input.peaks} \
            -r {input.case_bam} \
            -c {input.control_bam} \
            -g {config[rose_g]} \
            -o {output} 2> {log}
        """

################################################################################################################################

#                     Run ChromHMM for a subset of histone modifications                                                       #

#################################################################################################################################


rule bam2bed:
    input:
        "04aln_downsample/{sample}-downsample.sorted.bam"
    output:
        "12bed/{sample}.bed"
    params: jobname = "{sample}"
    log: "00log/{sample}_bam2bed.log"
    threads: 1
    conda: "envs/bedtools.yaml"
    resources:
        mem_mb=8000,
        runtime=120
    message: "converting bam to bed for {input}"
    shell:
        """
        bedtools bamtobed -i {input} > {output} 2> {log}
        """

if config["chromHMM"]:
    rule make_table:
        input : expand("12bed/{sample}.bed", sample = HISTONE_CASES + CONTROLS)
        output : "12bed/cellmarkfiletable.txt"
        log: "00log/make_table_chromHMM.log"
        message: "making a table for chromHMM"
        run:
            import os
            from os.path import join
            with open (output[0], "w") as f:
                for case in HISTONE_CASES:
                    sample = "_".join(case.split("_")[0:-1])
                    mark = case.split("_")[-1]
                    control = sample + "_" + CONTROL
                    case_bed = case + ".bed"
                    if os.path.exists(join("12bed", case_bed)):
                        f.write(sample + "\t" +  mark + "\t" + case + ".bed" + "\t" + control + ".bed" + "\n")

if config["chromHMM"]:
    rule chromHmm_binarize:
        input:
            beds = expand("12bed/{sample}.bed", sample = HISTONE_CASES + CONTROLS),
            table = "12bed/cellmarkfiletable.txt"
        output:
            directory("13chromHMM/binarizedData")
        log: "00log/chromhmm_bin.log"
        threads: 1
        conda: "envs/chromhmm.yaml"
        resources:
            mem_mb=16000,
            runtime=360
        shell:
            """
            # ChromHMM requires chromosome sizes file
            # Download if not present for your genome
            mkdir -p chromsizes
            if [ ! -f "chromsizes/{config[chromHmm_g]}.txt" ]; then
                echo "Warning: Chromosome sizes file not found for {config[chromHmm_g]}" > {log}
                echo "Please provide chromsizes/{config[chromHmm_g]}.txt" >> {log}
            fi

            ChromHMM.sh BinarizeBed -b {config[binsize]} \
                chromsizes/{config[chromHmm_g]}.txt \
                12bed/ {input.table} {output} 2> {log}
            """

if config["chromHMM"]:
    rule chromHmm_learn:
        input: "13chromHMM/binarizedData"
        output: directory("13chromHMM/MYOUTPUT")
        log: "00log/chromhmm_learn.log"
        threads: 10
        conda: "envs/chromhmm.yaml"
        resources:
            mem_mb=16000,
            runtime=720
        shell:
            """
            ChromHMM.sh LearnModel -p {threads} -b {config[binsize]} \
                {input} {output} {config[state]} {config[chromHmm_g]} 2> {log}
            """
