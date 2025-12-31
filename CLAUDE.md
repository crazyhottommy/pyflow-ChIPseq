# pyflow-ChIPseq: Repository Structure and Overview

**Analysis Date:** 2025-12-26
**Branch:** modernize-2025 (based on pairend)
**Purpose:** Documentation for understanding and modernizing the ChIP-seq pipeline

---

## Table of Contents

1. [Repository Purpose](#repository-purpose)
2. [Directory Structure](#directory-structure)
3. [Key Components](#key-components)
4. [Technology Stack](#technology-stack)
5. [Workflow Overview](#workflow-overview)
6. [Configuration System](#configuration-system)
7. [Modernization Priorities](#modernization-priorities)

---

## Repository Purpose

**pyflow-ChIPseq** is a Snakemake-based bioinformatics pipeline for processing and analyzing Chromatin Immunoprecipitation Sequencing (ChIP-seq) data. It automates the complete workflow from raw sequencing reads to peak calling and chromatin state analysis.

### Key Features

- Processes both public (GEO/SRA) and in-house ChIP-seq data
- Supports single-end and paired-end sequencing reads
- Handles both short (<70bp) and long (>70bp) reads
- Automated quality control and reporting
- Advanced downstream analysis (super-enhancers, chromatin states)
- Cluster deployment support (SLURM/DRMAA)

### Publication

This pipeline was published in the Journal of Visualized Experiments (JOVE):
*"An Integrated Platform for Genome-wide Mapping of Chromatin States Using High-throughput ChIP-sequencing in Tumor Tissues"*

---

## Directory Structure

```
pyflow-ChIPseq/
├── Snakefile                    # Main workflow definition (21,957 bytes)
├── config.yaml                  # Pipeline configuration and parameters
├── cluster.json                 # Cluster job submission resource specs
├── samples.json                 # Sample metadata (auto-generated)
│
├── sample2json.py               # Metadata converter (TSV → JSON)
├── sbatch_cluster.py            # SLURM job submission wrapper
├── pyflow-ChIPseq.sh            # Main execution script (SLURM)
├── pyflow-drmaa-ChIPseq.sh      # DRMAA-based execution (LSF)
├── jobscript.sh                 # Job script template
│
├── scripts/
│   └── sraDownload.R            # SRA data downloader
│
├── SRR.txt                      # Sample metadata template (GEO data)
├── meta.txt                     # Custom metadata template (in-house data)
│
├── README.md                    # User documentation
├── LICENSE                      # MIT License
├── rulegraph.png                # Workflow DAG visualization
├── GEO_rulegraph.png            # Alternative workflow visualization
└── TCGA_related files           # TCGA barcode documentation

### Output Directories (Created During Execution)

- `00log/` - Log files for all rules
- `01seq/` - Merged FASTQ files
- `02fqc/` - FastQC quality control reports
- `03aln/` - Aligned BAM files and indices
- `04aln_downsample/` - Downsampled BAM files
- `05phantompeakqual/` - Phantom peak quality metrics
- `06bigwig_inputSubtract/` - Input-subtracted bigWig tracks
- `07bigwig/` - RPKM-normalized bigWig tracks
- `08peak_macs1/` - MACS1 peak calls
- `09peak_macs2/` - MACS2 peak calls
- `10multiQC/` - MultiQC quality summary report
- `11superEnhancer/` - Super enhancer calls (ROSE)
- `12bed/` - BED format files
- `13chromHMM/` - Chromatin state predictions
```

---

## Key Components

### 1. Snakefile (Main Workflow)

The core orchestrator containing all pipeline rules, organized into:

#### Data Preparation Rules
- `merge_fastqs` - Combines multiple FASTQ files per sample/mark
- `fastqc` - Quality control on raw reads

#### Alignment Rules
- `align` - BWA mapping (supports both single/paired-end, short/long reads)
  - Short reads (<70bp): Uses `bwa aln + sampe/samse`
  - Long reads (>70bp): Uses `bwa mem`
  - Integrates `samblaster` for duplicate marking

#### BAM Processing Rules
- `flagstat_bam` - Counts mapped/unmapped reads
- `down_sample` - Normalizes read depth (default: 50M reads)

#### Quality Control Rules
- `phantom_peak_qual` - ChIP-seq quality assessment
- `multiQC` - Comprehensive quality report aggregation

#### Peak Calling Rules
- `call_peaks_macs1` - MACS v1.4.2 peak calling (narrow, nomodel)
- `call_peaks_macs2` - MACS v2.1.1 peak calling (broad peaks)

#### Visualization Rules
- `make_bigwigs` - RPKM-normalized bigWig generation
- `make_inputSubtract_bigwigs` - Input-subtracted bigWig tracks

#### Advanced Analysis Rules
- `superEnhancer` - Super-enhancer identification via ROSE
- `bam2bed` - BAM to BED conversion
- `chromHmm_binarize` - Binarize data for ChromHMM
- `chromHmm_learn` - Learn chromatin state models

### 2. sample2json.py (Python 3)

**Purpose:** Converts tab-delimited metadata into JSON format for Snakemake

**Input Format (meta.txt):**
```
sample_name     fastq_name              factor          reads
LKR10           Input_LKR10_1.fq.gz     Input           R1
LKR10           Input_LKR10_2.fq.gz     Input           R2
V6_5            Mll4_V6_5_1.fq.gz       Mll4            R1
```

**Output:** `samples.json` with nested structure:
```json
{
  "sample_name": {
    "factor": {
      "R1": ["path/to/file1.fq.gz"],
      "R2": ["path/to/file2.fq.gz"]
    }
  }
}
```

**Key Functions:**
- Walks directory to find all `.fq.gz` or `.fastq.gz` files
- Validates that declared files exist
- Groups files by sample, factor, and read direction

### 3. sbatch_cluster.py (Python 3)

**Purpose:** SLURM job submission wrapper

**Features:**
- Parses Snakemake job properties
- Reads cluster resource specs from `cluster.json`
- Constructs `sbatch` commands with proper resource allocation
- Manages memory, CPUs, time limits, job names, and log paths

### 4. pyflow-ChIPseq.sh (Bash)

**Purpose:** Main execution script

**Features:**
- Creates `sbatch_log/` directory for job logs
- Submits up to 99,999 jobs to SLURM cluster
- Uses 240-second latency wait for file system sync
- Supports dry-run (`-np`) and rerun (`-R`) options
- Integrates with `cluster.json` for resource management

### 5. scripts/sraDownload.R (R)

**Purpose:** Downloads SRA files from NCBI

**Features:**
- Uses SRAdb R package for metadata queries
- Integrates with Aspera Connect for fast downloads
- Converts SRA to FASTQ format

---

## Technology Stack

### Workflow Management
- **Snakemake** (Python 3) - Workflow orchestration

### Alignment & BAM Processing
- **BWA** - Burrows-Wheeler Aligner
- **samtools** v1.3.1 - BAM manipulation
- **samblaster** v0.1.22 - Duplicate marking
- **sambamba** - BAM subsampling

### Quality Control
- **FastQC** - Read quality assessment
- **MultiQC** - Aggregated QC reports
- **phantompeakqual** - ChIP-seq specific metrics

### Peak Calling
- **MACS1** v1.4.2 (Python 2.x) - Narrow peak calling
- **MACS2** v2.1.1 (Python 2.x) - Broad peak calling

### Visualization
- **deepTools** v2.3.3+ (`bamCoverage`, `bamCompare`) - BigWig generation

### Advanced Analysis
- **ROSE** - Super-enhancer detection
- **ChromHMM** - Chromatin state modeling
- **bedtools** - BED file operations

### Cluster Management
- **SLURM** - Job scheduler (sbatch)
- **DRMAA** - Alternative job control (LSF)

### Data Processing
- **Python 3** - Custom scripts
- **R** - SRA downloading
- **Conda** - Environment management

---

## Workflow Overview

### Primary Processing Pipeline

```
┌─────────────────┐
│ Raw FASTQ Files │ (R1/R2 for paired-end)
└────────┬────────┘
         ↓
┌────────────────┐
│  merge_fastqs  │ Combine multiple FASTQ files per sample
└────────┬───────┘
         ↓
┌────────────────┐
│     fastqc     │ Quality control assessment
└────────┬───────┘
         ↓
┌────────────────┐
│     align      │ BWA mapping (aln/mem based on read length)
└────────┬───────┘
         ↓
┌────────────────┐
│  flagstat_bam  │ Count aligned reads
└────────┬───────┘
         ↓
┌────────────────┐
│  down_sample   │ Normalize to target read depth (50M default)
└────────┬───────┘
         ↓
┌────────────────┐
│phantom_peak_qual│ Assess ChIP quality
└────────┬───────┘
         ↓
    ┌────┴────┬─────────────┬──────────────┬──────────────┐
    ↓         ↓             ↓              ↓              ↓
┌──────┐  ┌──────┐   ┌──────────┐  ┌──────────┐  ┌──────────┐
│bigwig│  │input │   │peaks_macs1│  │peaks_macs2│  │  super   │
│ RPKM │  │subtr │   │(narrow)   │  │(broad)    │  │ enhancer │
└──────┘  └──────┘   └──────────┘  └──────────┘  └──────────┘
    │         │             │              │              │
    └─────────┴─────────────┴──────────────┴──────────────┘
                            ↓
                     ┌──────────┐
                     │ multiQC  │ Final QC report
                     └──────────┘
```

### Optional ChromHMM Pathway

```
┌─────────────────┐
│ Downsampled BAM │
└────────┬────────┘
         ↓
┌────────────────┐
│    bam2bed     │ Convert BAM to BED
└────────┬───────┘
         ↓
┌────────────────┐
│  make_table    │ Create mark-sample table
└────────┬───────┘
         ↓
┌────────────────┐
│chromHmm_binarize│ Discretize read coverage
└────────┬───────┘
         ↓
┌────────────────┐
│ chromHmm_learn │ Train 15-state HMM model
└────────┬───────┘
         ↓
┌────────────────┐
│ Chromatin State│ Predictions
└────────────────┘
```

### Processing Modes

1. **From FASTQ** (`from_fastq: True`)
   - Start with raw sequencing reads
   - Complete alignment and processing

2. **From BAM** (`from_fastq: False`)
   - Start with pre-aligned BAM files
   - Skip to downsampling and downstream analysis

3. **GEO Data Processing**
   - Download SRA files via `sraDownload.R`
   - Convert SRA → FASTQ using `fastq-dump`
   - Process as standard FASTQ pipeline

4. **Custom In-house Data**
   - Provide metadata file (`meta.txt`)
   - Run `sample2json.py` to organize files
   - Process as standard FASTQ pipeline

### Sample Organization Logic

- **Sample Name:** Biological replicate or condition (e.g., "LKR10")
- **Factor:** Chromatin mark or TF (e.g., "H3K27ac", "Input")
- **Combined ID:** `{sample_name}_{factor}` creates output names
- **Control Matching:** Samples with `control` factor (default "Input") serve as background

---

## Configuration System

### config.yaml Parameters

**User Settings:**
- `email` - Notification email address

**Pipeline Mode:**
- `from_fastq` - Start from FASTQ (True) vs pre-aligned BAM (False)
- `paired_end` - Paired-end (True) vs single-end (False) reads
- `long_reads` - Long reads >70bp (True) vs short reads <70bp (False)

**Reference Genome:**
- `ref_fa` - Path to reference genome FASTA (e.g., `/path/to/mm9.fa`)
- `macs_g` - MACS1 genome size (e.g., "mm", "hs")
- `macs2_g` - MACS2 genome size (e.g., "mm", "hs")

**Control/Input:**
- `control` - Control sample factor name (default: "Input")

**Peak Calling Thresholds:**
- `macs_pvalue` - MACS1 p-value threshold (default: 1e-9)
- `macs2_pvalue` - MACS2 narrow peak p-value (default: 1e-9)
- `macs2_pvalue_broad` - MACS2 broad peak p-value (default: 1e-5)

**Tool Paths:**
- `phantom_path` - Path to phantompeakqualtools
- `rose_path` - Path to ROSE installation
- `chromHMM_path` - Path to ChromHMM JAR file

**Downsampling:**
- `downsample` - Enable downsampling (True/False)
- `target_reads` - Target read count (default: 50,000,000)

**ChromHMM Settings:**
- `chromHMM` - Enable chromatin state modeling (True/False)
- `histone_for_chromHMM` - Space-delimited histone marks (e.g., "H3K4me3 H3K27ac")

### cluster.json Resource Allocation

Defines per-rule cluster resources:
- **Memory:** 4GB - 32GB (depends on rule)
- **CPUs:** 1 - 10 cores
- **Time:** 20 minutes - 24 hours
- **Partition:** SLURM partition name

Example:
```json
{
  "align": {
    "__default__": {
      "mem": "32000",
      "time": "12:00:00",
      "cpus": "10"
    }
  }
}
```

---

## Modernization Priorities

### Critical Updates Needed

1. **Python 2 → Python 3 Migration**
   - MACS1 and MACS2 are Python 2.x based
   - Recommendation: Update to MACS3 (Python 3 compatible)

2. **Tool Version Updates**
   - samtools v1.3.1 → Latest (currently v1.19+)
   - deepTools v2.3.3 → Latest (currently v3.5+)
   - samblaster v0.1.22 → Latest (v0.1.26+)
   - BWA → Consider BWA-MEM2 for speed improvements

3. **Snakemake Best Practices**
   - Use `conda:` directive for per-rule environments
   - Migrate to `resources:` instead of cluster.json
   - Add `container:` directives for Docker/Singularity
   - Update to Snakemake v7+ syntax

4. **Conda Environment Management**
   - Create environment.yaml files
   - Pin tool versions for reproducibility
   - Separate environments per tool group

5. **Documentation Updates**
   - Add installation instructions
   - Create test dataset
   - Add troubleshooting guide
   - Update citations and references

6. **Code Quality**
   - Add type hints to Python scripts
   - Improve error handling
   - Add input validation
   - Add unit tests

7. **Performance Optimizations**
   - Enable Snakemake profiles
   - Use local rules for lightweight tasks
   - Implement shadow rules for temporary files
   - Add benchmark directives

8. **Security & Best Practices**
   - Remove hardcoded paths
   - Add input sanitization
   - Use environment variables for sensitive data
   - Update deprecated command-line flags

### Backward Compatibility Considerations

- Maintain support for existing config.yaml format
- Provide migration scripts for old metadata
- Keep output directory structure consistent
- Document breaking changes clearly

---

## Notes

- This document was generated on 2025-12-26 based on the `pairend` branch
- The `modernize-2025` branch is intended for updates while maintaining functionality
- Original pipeline design was robust for its time but requires updates for modern infrastructure
- Core workflow logic remains sound; focus modernization on tool versions and execution environment

---

**For questions or contributions, please see the main README.md or open an issue on GitHub.**
