# Modernization Summary - pyflow-ChIPseq (2025)

**Branch**: `modernize-2025`
**Base**: `pairend` branch
**Date**: December 2025
**Purpose**: Update the ChIP-seq pipeline to work with modern tools and infrastructure

---

## Overview

This document summarizes all changes made to modernize the pyflow-ChIPseq pipeline. The modernization ensures compatibility with current bioinformatics tools, Python 3, and modern cluster computing environments while maintaining the original functionality.

---

## Major Changes

### 1. Python 2 → Python 3 Migration

**Old Approach:**
- MACS1 v1.4.2 (Python 2.7)
- MACS2 v2.1.1 (Python 2.7)
- Required `source activate py27` in shell commands

**New Approach:**
- MACS3 v3.0.1+ (Python 3.9+)
- Native Python 3 compatibility throughout
- No environment activation in shell commands (handled by conda directives)

**Impact:**
- ✅ Compatible with modern Python ecosystems
- ✅ Better maintained and actively developed tools
- ✅ Improved performance and bug fixes

---

### 2. Tool Version Updates

| Tool | Old Version | New Version | Changes |
|------|-------------|-------------|---------|
| samtools | 1.3.1 | 1.19+ | Updated command syntax (`-Sb` → `-b`, added `-@` for threads) |
| deepTools | 2.3.3 | 3.5.5+ | `--ratio subtract` → `--operation subtract`, `--normalizeUsingRPKM` → `--normalizeUsing RPKM` |
| BWA | 0.7.x | 0.7.17+ | No major changes, but recommended update |
| samblaster | 0.1.22 | 0.1.26+ | Bug fixes and improvements |
| sambamba | < 1.0 | 1.0+ | Performance improvements |
| FastQC | 0.11.x | 0.12.1+ | Updated reports |
| MultiQC | 1.x | 1.19+ | Improved aggregation |
| bedtools | 2.x | 2.31.1+ | Bug fixes |
| ChromHMM | old | 1.25+ | Updated via conda |

---

### 3. Snakemake Modernization

#### Deprecated Features Removed

**Shell Prefix/Suffix** (Removed)
```python
# OLD
shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); exit $exitstat")

# NEW - removed, Snakemake handles this better now
```

**Module Load Commands** (Removed)
```python
# OLD
shell: """
    module load fastqc
    fastqc ...
"""

# NEW - use conda directives
conda: "envs/qc.yaml"
shell: """
    fastqc ...
"""
```

**cluster.json** (Replaced)
```python
# OLD
CLUSTER = json.load(open(config['CLUSTER_JSON']))
threads: CLUSTER["align"]["n"]

# NEW - use resources directive
threads: 10
resources:
    mem_mb=32000,
    runtime=720
```

#### New Features Added

**Conda Environment Integration**
```python
rule example:
    conda: "envs/alignment.yaml"  # Automatic environment management
    shell: "..."
```

**Resource Specifications**
```python
resources:
    mem_mb=16000,    # Memory in MB
    runtime=240      # Runtime in minutes
```

**Named Outputs**
```python
# OLD
output: "03aln/{sample}.sorted.bam", "03aln/{sample}.sorted.bam.bai"

# NEW
output:
    bam = "03aln/{sample}.sorted.bam",
    bai = "03aln/{sample}.sorted.bam.bai"
```

---

### 4. Conda Environment Management

Created 7 isolated conda environments for different tool groups:

1. **alignment.yaml** - BWA, samtools, samblaster, sambamba
2. **qc.yaml** - FastQC, MultiQC
3. **peakcalling.yaml** - MACS3
4. **bigwig.yaml** - deepTools (bamCoverage, bamCompare)
5. **phantompeakqualtools.yaml** - R, ChIPpeakAnno, required R packages
6. **bedtools.yaml** - bedtools
7. **chromhmm.yaml** - ChromHMM, OpenJDK

**Benefits:**
- ✅ Reproducible software versions
- ✅ Isolated dependencies (no conflicts)
- ✅ Easy installation (automatic via Snakemake)
- ✅ Portable across systems

---

### 5. Cluster Execution Modernization

**Old Approach (cluster.json):**
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

**New Approach (Snakemake Profiles):**
```yaml
# profiles/slurm/config.yaml
jobs: 100
use-conda: True
conda-frontend: "mamba"
cluster: "sbatch --parsable --mem={resources.mem_mb} ..."
default-resources:
  - mem_mb=8000
  - runtime=120
```

**Usage:**
```bash
# OLD
snakemake --cluster-config cluster.json --cluster '...' -j 99

# NEW
snakemake --profile profiles/slurm
```

---

### 6. Configuration File Improvements

**config.yaml** changes:

- ✅ Removed hardcoded user-specific paths
- ✅ Added comprehensive comments and documentation
- ✅ Organized into logical sections
- ✅ Changed MACS p-values to q-values (FDR) for MACS3
- ✅ Updated default parameters for modern standards
- ✅ Made paths configurable with examples

**Key Parameter Updates:**
```yaml
# OLD - p-value (not FDR corrected)
macs_pvalue: 1e-9

# NEW - q-value (FDR corrected, more appropriate)
macs_pvalue: 0.05
```

---

### 7. Peak Calling Updates

**MACS1/MACS2 → MACS3 Migration**

Changes:
- Unified command: `macs3 callpeak` (instead of `macs` and `macs2 callpeak`)
- Changed from `-p` (p-value) to `-q` (q-value/FDR) as primary threshold
- Updated output file extensions:
  - MACS1: `.bed` → MACS3: `.narrowPeak`
  - MACS2: `.xls` → MACS3: `.broadPeak` + `.xls`
- Removed separate "nomodel" calls (MACS3 handles this internally)

**Output Directories:**
- `08peak_macs1/` → `08peak_macs3/` (narrow peaks)
- `09peak_macs2/` → `09peak_macs3/` (broad peaks)

---

### 8. Improved Error Handling

**Updated Commands:**
```bash
# OLD - basic samtools without thread support
samtools sort -m 2G -@ 5 -T {output[0]}.tmp -o {output[0]}

# NEW - modern samtools with proper threading
samtools sort -m 2G -@ {threads} -T {output.bam}.tmp -o {output.bam}
samtools index -@ {threads} {output.bam}
```

**Improved deepTools:**
```bash
# OLD
bamCompare --bamfile1 {input[1]} --bamfile2 {input[0]} --normalizeUsingRPKM --ratio subtract -p 5

# NEW
bamCompare --bamfile1 {input.case_bam} --bamfile2 {input.control_bam} \
    --normalizeUsing RPKM --operation subtract \
    --numberOfProcessors {threads}
```

---

### 9. Documentation Enhancements

**New Files:**
1. **INSTALL.md** - Comprehensive installation guide
2. **MODERNIZATION.md** - This file, documenting all changes
3. **CLAUDE.md** - Repository structure analysis
4. **profiles/slurm/README.md** - SLURM profile documentation

**Updated Files:**
- `config.yaml` - Extensive inline documentation
- `README.md` - Updated for modern usage

---

### 10. Code Quality Improvements

**Named Inputs/Outputs:**
```python
# Easier to read and maintain
input:
    control_bam = "...",
    case_bam = "..."
output:
    peaks = "...",
    xls = "..."
```

**Consistent Threading:**
```python
# OLD - hardcoded
bwa mem -t 5 ...

# NEW - uses Snakemake threads
bwa mem -t {threads} ...
```

**Better Resource Management:**
- All rules now have explicit resource requirements
- Default resources defined in profile
- Easier to tune for different cluster environments

---

## Breaking Changes

### 1. Output File Locations and Names

| Old | New | Notes |
|-----|-----|-------|
| `08peak_macs1/{sample}_macs1_peaks.bed` | `08peak_macs3/{sample}_narrow_peaks.narrowPeak` | MACS3 narrow peaks |
| `08peak_macs1/{sample}_macs1_nomodel_peaks.bed` | Combined into above | MACS3 handles nomodel internally |
| `09peak_macs2/{sample}_macs2_peaks.xls` | `09peak_macs3/{sample}_broad_peaks.broadPeak` | MACS3 broad peaks |

### 2. Configuration Parameters

```yaml
# Parameter name changes
macs_pvalue: 1e-8  → macs_pvalue: 0.05  # Now expects q-value
phantom_path: /path/to/run_spp.R  → Removed # Auto-downloaded
chromHMM_path: /path/to/ChromHMM/ → Removed # Conda-managed
```

### 3. Execution Method

```bash
# OLD
./pyflow-ChIPseq.sh

# NEW
snakemake --profile profiles/slurm
# or
snakemake --use-conda --cores 8
```

---

## Backward Compatibility

### What Remains Compatible

- ✅ Sample metadata format (`meta.txt`)
- ✅ `sample2json.py` script (already Python 3)
- ✅ `sbatch_cluster.py` script (already Python 3)
- ✅ Directory structure for outputs
- ✅ Input file organization
- ✅ Core workflow logic

### What Requires Updates

- ❌ Old conda environments (need Python 3)
- ❌ cluster.json (replaced by profiles)
- ❌ Execution scripts (use Snakemake directly)
- ❌ Hardcoded paths in config.yaml

---

## Migration Guide

### For Existing Users

1. **Update config.yaml**
   ```bash
   # Review and update these settings
   - ref_fa: point to your genome
   - macs_pvalue: change from 1e-9 to 0.05 (or desired q-value)
   - Remove phantom_path, chromHMM_path (now conda-managed)
   ```

2. **Create Conda Environments**
   ```bash
   snakemake --use-conda --conda-create-envs-only --cores 1
   ```

3. **Test with Dry Run**
   ```bash
   snakemake -n --use-conda
   ```

4. **Run Pipeline**
   ```bash
   # Local
   snakemake --use-conda --cores 8

   # SLURM cluster
   snakemake --profile profiles/slurm
   ```

---

## Performance Improvements

### Conda/Mamba

- Using Mamba (faster conda alternative) reduces environment creation time
- Pre-creating environments recommended for cluster runs

### Threading

- Better thread utilization with parameterized `{threads}`
- Scales with available resources

### Resource Allocation

- More accurate resource requests reduce job failures
- Better cluster utilization

---

## Testing

### Validation Steps Performed

1. ✅ Snakefile syntax validation
2. ✅ Conda environment creation
3. ✅ Dry run completion
4. ✅ Tool version compatibility checks
5. ✅ Configuration parsing

### Recommended Testing

Before running on production data:

```bash
# 1. Dry run
snakemake -n --use-conda

# 2. Create environments
snakemake --use-conda --conda-create-envs-only --cores 1

# 3. Test on subset
# Limit to first few rules
snakemake --use-conda --cores 4 02fqc/test_sample_R1_fastqc.zip
```

---

## Future Enhancements (Recommended)

### High Priority

1. **Container Support** - Add Docker/Singularity containers for full reproducibility
2. **Automated Testing** - Add CI/CD with test datasets
3. **ROSE Alternative** - Replace Python 2.7 dependent ROSE with modern alternative

### Medium Priority

1. **BWA-MEM2** - Switch to faster BWA-MEM2 aligner
2. **Benchmark Directives** - Add runtime and memory benchmarking
3. **Report Generation** - Automated MultiQC-style reports for peak calls

### Low Priority

1. **Checkpoints** - Dynamic rule execution based on results
2. **Module System** - Modularize Snakefile for easier maintenance
3. **Parameter Optimization** - Auto-tune parameters based on data quality

---

## Acknowledgments

- Original pipeline: Ming Tang (https://github.com/crazyhottommy)
- Modernization: Claude Code (2025)
- Community: Snakemake developers and bioinformatics community

---

## Support and Issues

For questions or issues with the modernized version:

1. Check **INSTALL.md** for setup help
2. Review **CLAUDE.md** for architecture understanding
3. Open an issue on GitHub
4. Consult Snakemake documentation

---

## Summary of Files Changed

### Modified Files
- `Snakefile` - Complete modernization
- `config.yaml` - Reorganized and documented
- `.gitignore` - Add Snakemake artifacts (if needed)

### New Files
- `envs/*.yaml` - 7 conda environment files
- `profiles/slurm/config.yaml` - SLURM profile
- `profiles/slurm/README.md` - Profile documentation
- `INSTALL.md` - Installation guide
- `MODERNIZATION.md` - This document
- `CLAUDE.md` - Repository structure documentation

### Deprecated Files (Still Present)
- `cluster.json` - Use profiles instead
- `pyflow-ChIPseq.sh` - Use `snakemake --profile` instead
- `pyflow-drmaa-ChIPseq.sh` - Legacy DRMAA support

---

**The pipeline is now ready for modern ChIP-seq analysis!**
