# Installation Guide for pyflow-ChIPseq (Modernized 2025)

This guide will help you set up the modernized pyflow-ChIPseq pipeline on your system.

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Reference Genome Setup](#reference-genome-setup)
- [Configuration](#configuration)
- [Testing the Installation](#testing-the-installation)
- [Troubleshooting](#troubleshooting)

---

## Prerequisites

### Required Software

1. **Conda or Mamba** (for environment management)
   - Miniconda: https://docs.conda.io/en/latest/miniconda.html
   - Mamba (faster): `conda install -c conda-forge mamba`

2. **Snakemake >= 7.0**
   ```bash
   conda create -n snakemake -c bioconda -c conda-forge snakemake>=7.0
   conda activate snakemake
   ```

3. **Git** (for cloning the repository)
   ```bash
   # Usually pre-installed on Linux/Mac
   git --version
   ```

### System Requirements

- **Memory**: Minimum 32 GB RAM recommended (alignment and peak calling are memory-intensive)
- **Storage**: Depends on dataset size
  - Raw FASTQ: ~10-50 GB per sample
  - BAM files: ~5-20 GB per sample
  - Final outputs: ~2-10 GB per sample
- **CPU**: Multi-core processor (8+ cores recommended for parallel processing)

---

## Installation

### Step 1: Clone the Repository

```bash
# Clone from GitHub
git clone https://github.com/crazyhottommy/pyflow-ChIPseq.git
cd pyflow-ChIPseq

# Checkout the modernized branch
git checkout modernize-2025
```

### Step 2: Install Snakemake (if not already installed)

```bash
# Create a conda environment for Snakemake
conda create -n snakemake python=3.11 -y
conda activate snakemake

# Install Snakemake and mamba
conda install -c bioconda -c conda-forge snakemake mamba -y
```

### Step 3: Verify Installation

```bash
# Check Snakemake version (should be >= 7.0)
snakemake --version

# Check conda/mamba
conda --version
mamba --version
```

---

## Reference Genome Setup

The pipeline requires a reference genome in FASTA format and BWA indices.

### Option 1: Download Pre-built Genomes

#### Mouse (mm10)
```bash
# Create genome directory
mkdir -p genomes/mm10
cd genomes/mm10

# Download from UCSC
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
# or curl -O http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gunzip mm10.fa.gz

# Build BWA index
conda activate snakemake
mamba install -c bioconda bwa
bwa index mm10.fa

cd ../..
```

#### Human (hg38)
```bash
# Create genome directory
mkdir -p genomes/hg38
cd genomes/hg38

# Download from UCSC
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Build BWA index
mamba install bwa
bwa index hg38.fa

cd ../..
```

### Option 2: Use Existing Reference

If you already have a reference genome:

```bash
# Just create a symbolic link
ln -s /path/to/existing/genome.fa genomes/genome.fa
```

---

## Configuration

### Step 1: Prepare Your Metadata File

Create a tab-delimited file (`meta.txt`) describing your samples:

```
sample_name	fastq_name	factor	reads
Sample1	Sample1_Input_R1.fastq.gz	Input	R1
Sample1	Sample1_Input_R2.fastq.gz	Input	R2
Sample1	Sample1_H3K27ac_R1.fastq.gz	H3K27ac	R1
Sample1	Sample1_H3K27ac_R2.fastq.gz	H3K27ac	R2
Sample2	Sample2_Input_R1.fastq.gz	Input	R1
Sample2	Sample2_Input_R2.fastq.gz	Input	R2
Sample2	Sample2_H3K4me3_R1.fastq.gz	H3K4me3	R1
Sample2	Sample2_H3K4me3_R2.fastq.gz	H3K4me3	R2
```

**Column descriptions:**
- `sample_name`: Biological sample identifier
- `fastq_name`: FASTQ file name (must match actual files)
- `factor`: Chromatin mark or transcription factor (use "Input" for controls)
- `reads`: R1 for forward reads, R2 for reverse reads (paired-end)

For single-end data, only use R1.

### Step 2: Generate samples.json

```bash
# Run the metadata converter
python sample2json.py /path/to/fastq/directory meta.txt

# This creates samples.json which the pipeline uses
```

### Step 3: Edit config.yaml

Update the configuration file with your settings:

```yaml
# Adjust these settings for your experiment
from_fastq: True
paired_end: True
long_reads: True

# IMPORTANT: Update this path
ref_fa: genomes/mm10/mm10.fa

# Genome size for MACS3
macs_g: mm
macs2_g: mm

# Control sample name (must match your metadata)
control: 'Input'

# Peak calling stringency
macs_pvalue: 0.05
macs2_pvalue: 0.05
macs2_pvalue_broad: 0.1

# Downsampling
downsample: True
target_reads: 50000000

# Optional analyses
chromHMM: False
```

---

## Testing the Installation

### Step 1: Create Log Directories

```bash
mkdir -p logs/slurm
mkdir -p 00log
```

### Step 2: Dry Run (Local)

Test the workflow without actually running jobs:

```bash
# Activate snakemake environment
conda activate snakemake

# Dry run to check for errors
snakemake -n --use-conda
```

### Step 3: Small Test Run (Local)

Run on a small subset of data locally:

```bash
# Run with 4 cores
snakemake --use-conda --cores 4 --keep-going

# Or limit to specific targets
snakemake --use-conda --cores 4 02fqc/Sample1_Input_R1_fastqc.zip
```

### Step 4: Cluster Execution (SLURM)

If you have access to a SLURM cluster:

```bash
# Dry run with SLURM profile
snakemake --profile profiles/slurm -n

# Full run on cluster
snakemake --profile profiles/slurm
```

---

## Workflow Execution Options

### Local Execution (No Cluster)

```bash
# Use all available cores
snakemake --use-conda --cores all

# Limit to 8 cores
snakemake --use-conda --cores 8
```

### SLURM Cluster Execution

```bash
# Using the SLURM profile (recommended)
snakemake --profile profiles/slurm

# Manual SLURM submission (legacy)
snakemake --use-conda --cluster "sbatch --mem={resources.mem_mb} --time={resources.runtime}" --jobs 100
```

### Re-running Failed Jobs

```bash
# Rerun incomplete jobs
snakemake --use-conda --cores 8 --rerun-incomplete

# Force re-run specific rule
snakemake --use-conda --cores 8 --forcerun call_peaks_macs3_narrow
```

---

## Environment Management

The pipeline automatically creates conda environments for each tool group. To pre-create them:

```bash
# Create all environments in advance (recommended)
snakemake --use-conda --conda-create-envs-only --cores 1

# This creates environments in .snakemake/conda/
```

To clean up conda environments:

```bash
# Remove all conda environments
snakemake --use-conda --cleanup-conda
```

---

## Troubleshooting

### Issue: "samples.json not found"

**Solution**: Run `sample2json.py` to generate it:
```bash
python sample2json.py /path/to/fastq meta.txt
```

### Issue: "Reference genome not indexed"

**Solution**: Build BWA index:
```bash
bwa index genomes/mm10/mm10.fa
```

### Issue: Conda environment creation fails

**Solution**: Use mamba instead of conda:
```bash
conda install -c conda-forge mamba
# Then use --conda-frontend mamba with snakemake
```

### Issue: SLURM jobs fail with memory errors

**Solution**: Increase memory in the rule's `resources:` directive in Snakefile:
```python
resources:
    mem_mb=32000  # Increase this value
```

### Issue: "ChromHMM.jar not found"

**Solution**: ChromHMM is installed via conda environment. If using manually:
```bash
# Download ChromHMM
wget https://compbio.mit.edu/ChromHMM/ChromHMM.zip
unzip ChromHMM.zip
```

### Issue: Permissions error with phantompeakqualtools

**Solution**: The pipeline auto-downloads run_spp.R. Ensure write permissions:
```bash
chmod +w .
```

---

## Next Steps

After successful installation:

1. **Prepare your data**: Organize FASTQ files in a directory
2. **Create metadata**: Write meta.txt describing your samples
3. **Configure pipeline**: Edit config.yaml with your settings
4. **Generate samples.json**: Run sample2json.py
5. **Dry run**: Test with `snakemake -n`
6. **Execute**: Run the pipeline!

---

## Support

For issues and questions:
- **GitHub Issues**: https://github.com/crazyhottommy/pyflow-ChIPseq/issues
- **Original README**: See README.md for usage examples
- **Snakemake Documentation**: https://snakemake.readthedocs.io/

---

## Key Improvements in Modernized Version

- ✅ Python 3 compatible (MACS3 instead of MACS1/MACS2)
- ✅ Modern Snakemake syntax (v7+)
- ✅ Conda environments for reproducibility
- ✅ Updated tool versions (samtools 1.19+, deepTools 3.5+)
- ✅ SLURM profiles (replaces cluster.json)
- ✅ Resource specifications in Snakefile
- ✅ Improved error handling and logging
- ✅ Better documentation

Enjoy your modernized ChIP-seq pipeline!
