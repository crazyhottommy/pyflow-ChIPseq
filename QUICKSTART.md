# Quick Start Guide - pyflow-ChIPseq (Modernized 2025)

Get up and running with the modernized ChIP-seq pipeline in minutes!

## Prerequisites

- Conda or Mamba installed
- Snakemake >= 7.0
- FASTQ files from your ChIP-seq experiment

```bash
  # Create all conda environments without running the pipeline
  snakemake --use-conda --conda-create-envs-only --cores 1
```

## 5-Minute Setup

### 1. Download data 

Download the fastq files using fastq-dl

https://github.com/rpetit3/fastq-dl

```bash
conda install -c conda-forge -c bioconda fastq-dl

fastq-dl --accession SRR2518123
fastq-dl --accession SRR2518124 
fastq-dl --accession SRR2518125 
fastq-dl --accession SRR2518126 
```


### 2. Prepare Your Data

Create a metadata file `BRD4_meta.txt:

```bash
cat BRD4_meta.txt
sample_name fastq_name  factor  reads
MOLM-14_DMSO1_5 SRR2518123.fastq.gz BRD4    R1
MOLM-14_DMSO1_5 SRR2518124.fastq.gz Input   R1
MOLM-14_DMSO2_6 SRR2518125.fastq.gz BRD4    R1
MOLM-14_DMSO2_6 SRR2518126.fastq.gz Input   R1
```

This is a single-end example.
for paired-end data, you will have R1 and R2 for the same sample name 
in separate rows.


### 3. Generate Sample JSON

```bash
# see help
python sample2json.py -h 

# create the samples.json file 
python sample2json.py /folder/path/to/fastq/ BRD4_meta.txt
```

This creates `samples.json` that the pipeline uses.

```bash
cat samples.json 

{
    "MOLM-14_DMSO1_5": {
        "BRD4": {
            "R1": [
                "/Users/tommytang/githup_repo/pyflow-ChIPseq/data/SRR2518123.fastq.gz"
            ]
        },
        "Input": {
            "R1": [
                "/Users/tommytang/githup_repo/pyflow-ChIPseq/data/SRR2518124.fastq.gz"
            ]
        }
    },
    "MOLM-14_DMSO2_6": {
        "BRD4": {
            "R1": [
                "/Users/tommytang/githup_repo/pyflow-ChIPseq/data/SRR2518125.fastq.gz"
            ]
        },
        "Input": {
            "R1": [
                "/Users/tommytang/githup_repo/pyflow-ChIPseq/data/SRR2518126.fastq.gz"
            ]
        }
    }
}
```

It is a dictionary of dictionaries 

### 4. Configure the Pipeline

Edit `config.yaml`:

```yaml
# Essential settings
from_fastq: True
paired_end: False
long_reads: True

# Update this path!
ref_fa: /path/to/genome.fa

# Genome for MACS3
macs_g: mm  # or 'hs' for human

# Control name (must match meta.txt)
control: 'Input'
```

### 5. Run the Pipeline

#### Option A: Local Execution

```bash
# Dry run (check for errors)
snakemake -n --use-conda

# Run with 8 cores
snakemake --use-conda --cores 8 --keep-going
```

#### Option B: SLURM Cluster

```bash
# Create log directory
mkdir -p logs/slurm

# Dry run
snakemake --profile profiles/slurm -n

# Full run
snakemake --profile profiles/slurm
```

## Expected Outputs

After successful completion, you'll have:

```
00log/              # Log files
02fqc/              # FastQC reports
03aln/              # Aligned BAM files
04aln_downsample/   # Downsampled BAMs (50M reads)
05phantompeakqual/  # ChIP quality metrics
06bigwig_inputSubtract/  # Input-subtracted bigWigs
07bigwig/           # RPKM-normalized bigWigs
08peak_macs3/       # Narrow peaks
09peak_macs3/       # Broad peaks
10multiQC/          # Quality summary report
11superEnhancer/    # Super-enhancer calls (optional)
```

## Viewing Results

### 1. Quality Control

```bash
# Open MultiQC report in browser
open 10multiQC/multiQC_log.html
```

### 2. Peak Files

```bash
# Narrow peaks (for TFs, H3K4me3, etc.)
ls 08peak_macs3/*_narrow_peaks.narrowPeak

# Broad peaks (for H3K27me3, H3K9me3, etc.)
ls 09peak_macs3/*_broad_peaks.broadPeak
```

### 3. Visualization Tracks

```bash
# BigWig files for IGV/UCSC Genome Browser
ls 07bigwig/*.bw
ls 06bigwig_inputSubtract/*.bw
```

## Common Workflows

### Just QC and Alignment

```bash
snakemake --use-conda --cores 8 10multiQC/multiQC_log.html
```

### Specific Sample

```bash
snakemake --use-conda --cores 8 \
    08peak_macs3/Sample1_H3K27ac_vs_Sample1_Input_narrow_peaks.narrowPeak
```

### Re-run Failed Jobs

```bash
snakemake --use-conda --cores 8 --rerun-incomplete
```

## Troubleshooting

### "samples.json not found"
```bash
python sample2json.py /path/to/fastq meta.txt
```

### "Reference genome not indexed"
```bash
bwa index /path/to/genome.fa
```

### Jobs fail with memory errors
Edit `Snakefile` and increase `mem_mb` in the failing rule's `resources:` section.

### Conda is slow
Use mamba instead:
```bash
conda install -c conda-forge mamba
# Then use --conda-frontend mamba
```

## Next Steps

- **Quality Check**: Review MultiQC report
- **Peak Calling**: Adjust q-values in config.yaml if needed
- **Downstream Analysis**: Use peaks for motif discovery, annotation, etc.
- **Visualization**: Load bigWig files in IGV

## More Help

- **Full installation**: See [INSTALL.md](INSTALL.md)
- **All changes**: See [MODERNIZATION.md](MODERNIZATION.md)
- **Architecture**: See [CLAUDE.md](CLAUDE.md)

## One-Line Examples

```bash
# Dry run
snakemake -n --use-conda

# Run locally with 16 cores
snakemake --use-conda --cores 16 --keep-going

# Run on SLURM
snakemake --profile profiles/slurm

# Pre-create conda environments
snakemake --use-conda --conda-create-envs-only --cores 1

# Clean up and restart
rm -rf .snakemake/ 03aln/ 04aln_downsample/ 08peak_macs3/ 09peak_macs3/
snakemake --use-conda --cores 8
```

---

**Happy ChIP-seqing! 🧬🔬**
