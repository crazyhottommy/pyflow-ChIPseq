# SLURM Profile for pyflow-ChIPseq

This profile enables running the pipeline on SLURM clusters using modern Snakemake features.

## Usage

To run the pipeline with this profile:

```bash
# Dry run to check the workflow
snakemake --profile profiles/slurm -n

# Run the pipeline
snakemake --profile profiles/slurm
```

## Configuration

The `config.yaml` file in this directory contains cluster-specific settings:

- **jobs**: Maximum number of concurrent jobs
- **use-conda**: Automatically create and use conda environments
- **conda-frontend**: Use mamba (faster) or conda for environment management
- **cluster**: SLURM submission command template
- **default-resources**: Default memory, runtime, and thread allocations

## Customization

### Per-Rule Resources

Resources are now defined in the Snakefile using the `resources:` directive:

```python
rule example:
    input: ...
    output: ...
    threads: 8
    resources:
        mem_mb=16000,    # 16 GB memory
        runtime=240      # 4 hours in minutes
    conda: "envs/example.yaml"
    shell: "..."
```

### Cluster Partitions

To specify different SLURM partitions for different rules, you can modify the cluster command or create partition-specific profiles.

Example for adding partition:
```yaml
cluster: "sbatch --parsable --partition=general --mem={resources.mem_mb} ..."
```

### Log Files

SLURM output and error logs will be written to:
- `logs/slurm/{rule}_{wildcards}_{jobid}.out`
- `logs/slurm/{rule}_{wildcards}_{jobid}.err`

Make sure to create the logs directory:
```bash
mkdir -p logs/slurm
```

## Requirements

- Snakemake >= 7.0
- SLURM cluster environment
- Conda or Mamba installed

## Alternative Execution

### Local Execution (No Cluster)

To run locally without SLURM:
```bash
snakemake --use-conda --cores 8
```

### With Conda Only (No Cluster)

```bash
snakemake --use-conda --cores 16 --keep-going
```

## Troubleshooting

### Jobs Fail with Memory Errors

Increase `mem_mb` in the specific rule's resources directive.

### Jobs Timeout

Increase `runtime` (in minutes) in the specific rule's resources directive.

### Environment Creation Fails

Try using mamba instead of conda:
```bash
conda install -c conda-forge mamba
```

Then ensure `conda-frontend: "mamba"` is set in the profile config.
