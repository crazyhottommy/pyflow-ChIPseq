# GEOpyflow-ChIPseq
a snakemake pipeline to process ChIP-seq files from GEO

I have so many people asking me to process a public GEO ChIP-seq data set for them. I hate to repeat the same steps and decide to make a pipeline for it.

Now, this is working on LSF, I will have another branch for SLURM.

### work flow of the pipeline

![](./rule_diagram.png)


### Dependiencies

* [snakemake](https://bitbucket.org/snakemake/snakemake). snakemake is python3
* R > 3.3.0
you will need `optparse` package. `install.packages("optparse")`
`SRAdb`

```r
source("https://bioconductor.org/biocLite.R")
biocLite("SRAdb")
```

```
 Rscript sraDownload.R  -a ascp -QT -l 300m -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -t fastq SRR3144652
```

* aspera for downloading

check this blog post by MARK ZIEMANN http://genomespot.blogspot.com/2015/05/download-sra-data-with-aspera-command.html

```bash
sh <(curl -s aspera-connect-3.6.2.117442-linux-64.sh)
```

`sraDownload.R` is under the `scripts` folder from [Luke Zappia](https://github.com/lazappi):

```bash
## single quote your ascp command, otherwise R will be confused
Rscript sraDownload.R  -a 'ascp -QT -l 300m -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh' -t fastq SRR3144652

```

* [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* `bamCoverage` v2.3.3 from [deeptools](https://github.com/fidelram/deepTools) for making RPKM normalized and input subtracted bigwig files
* [bowtie1](http://bowtie-bio.sourceforge.net/index.shtml) for aligning short reads (< 50bp)
* [samtools](http://www.htslib.org/) v1.3.1
* [ROSE](http://younglab.wi.mit.edu/super_enhancer_code.html) for calling superEnhancer. ROSE has to be run inside the installation folder. now I hard coded the path in the Snakefile. (you will have to change that to the ROSE directory in your cluster). Todo: expose the path to the `config.yaml` file that one can change.
* [macs1](https://pypi.python.org/pypi/MACS/1.4.2) v1.4.2 and [macs2](https://github.com/taoliu/MACS) v2.1.1 for calling peaks (macs2 for broad peaks).
* [multiQC](http://multiqc.info/)

`macs1`, `macs2` and `ROSE` are python2.x, see this [Using Snakemake and Macs2](https://groups.google.com/forum/#!searchin/snakemake/macs%7Csort:relevance/snakemake/60txGSq81zE/NzCUTdJ_AQAJ) in the snakemake google group.

if you look inside the `Snakefile`, I did `source activate root` back to python2.x before running macs1 and macs2.

There will be [Integration of conda package management into Snakemake](https://bitbucket.org/snakemake/snakemake/pull-requests/92/wip-integration-of-conda-package/diff)


### How to distribute workflows

read [doc](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html)

```bash
ssh shark.mdanderson.org

# start a screen session
screen

# make a folder, name it yourself, I named it workdir for demon
mkdir /rsch2/genomic_med/krai/workdir/

cd /rsch2/genomic_med/krai/workdir/

git clone https://github.com/crazyhottommy/GEOpyflow-ChIPseq

cd GEOpyflow-ChIPseq

## go to downsampling branch. shark is LSF system
git checkout shark

## edit the config.yaml file as needed, e.g. set mouse or human for ref genome, p value cut off for peak calling, the number of reads you want to downsample to
nano config.yaml

## skip this if on Shark, samir has py351 set up for you. see below STEPS
conda create -n snakemake python=3 -c bioconda multiqc snakemake deeptools
source activate snakemake
```

## STEPS 

### Download the sra files

Prepare a txt file `SRR.txt` which has two columns: IP and Input:

e.g.

```bash
cat SRR.txt

IP      Input
SRR3144652      SRR3144654
SRR3144653      SRR3144654
SRR2518123      SRR2518124
SRR2518125      SRR2518126
SRR2518127      SRR2518128
SRR2518129      SRR2518130
SRR1616137      SRR1616139
SRR1616138      SRR1616139
SRR1616140      SRR1616142
SRR1616141      SRR1616142


```


```bash
cd GEOpyflow-ChIPseq
mkdir 01seq
cd 01seq
## because fastq-dump only convert the sra files to fastq to the current folder. one has to download the sra files to the 01seq folder
```

make a shell script:
`download.sh`

```bash
# /bin/bash
set -euo pipefail

## you will need to change the ascp command to get the right path
Rscript ../scripts/sraDownload.R -a 'ascp -QT -l 300m -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh' $1
```

`chmod u+x download.sh`

```bash
# inside the GEOpyflow-ChIPseq folder:
cat ../SRR.txt | sed '1d' | tr "\t" "\n" | sort | uniq > srr_unique.txt

## only have 4 jobs in parallel, good behavior on a cluster
cat srr_unique.txt | parallel -j 4 ./download.sh {}

# all the sra files will be downloaded in the current folder.
```

Now you have all `sra` files downloaded into `01seq` folder, proceed below:

### dry run to test 

```bash
## dry run
snakemake -np

## test for one sample
./pyflow-ChIPseq.sh  07bigwig/SRR2518125.bw

```

if no errors, preceed below.

### Using [DRMAA](https://www.drmaa.org/)

[job control through drmaa](http://drmaa-python.readthedocs.io/en/latest/tutorials.html#controlling-a-job)

DRMAA is only supported on `Shark` (LSF system).

```bash
module load drmma
./pyflow-drmaa-ChIPseq.sh
```

Using `drmaa` can `control + c` to stop the current run.

Dependent jobs are submitted one by one, if some jobs failed, the pipeline will stop. Good for initital testing.

### submit all jobs to the cluster

```bash
./pyflow-ChIPseq.sh 
```

All jobs will be submitted to the cluster on queue.  This is useful if you know your jobs will succeed for most of them and the jobs are on queue to gain priority.

### job control

To kill all of your pending jobs you can use the command:

```bash
bkill `bjobs -u krai |grep PEND |cut -f1 -d" "`
```

other useful commands: 

```
bjobs -pl
Display detailed information of all pending jobs of the invoker.

bjobs -ps
Display only pending and suspended jobs.

bjobs -u all -a
Display all jobs of all users.

bjobs -d -q short -m apple -u mtang1
Display all the recently finished jobs submitted by john to the
queue short, and executed on the host apple.
```

### rerun some of the jobs

```bash

# specify the name of the rule, all files that associated with that rule will be rerun. e.g. rerun macs2 calling peaks rule,
./pyflow-ChIPseq -R call_peaks_macs2

## rerun one sample, just specify the name of the target file

./pyflow-ChIPseq -R SRR3144652.sorted.bam

# check snakemake -f, -R, --until options
./pyflow-ChIPseq -f call_peaks_macs2
```

### checking results after run finish

```bash

snakemake --summary | sort -k1,1 | less -S

# or detailed summary will give you the commands used to generated the output and what input is used
snakemake --detailed-summary | sort -k1,1 > snakemake_run_summary.txt
```


### clean the folders

I use echo to see what will be removed first, then you can remove all later.

```
find . -maxdepth 1 -type d -name "[0-9]*" | xargs echo rm -rf
```


### Snakemake does not trigger re-runs if I add additional input files. What can I do?

Snakemake has a kind of “lazy” policy about added input files if their modification date is older than that of the output files. One reason is that information what to do cannot be inferred just from the input and output files. You need additional information about the last run to be stored. Since behaviour would be inconsistent between cases where that information is available and where it is not, this functionality has been encoded as an extra switch. To trigger updates for jobs with changed input files, you can use the command line argument –list-input-changes in the following way:

```bash
snakemake -n -R `snakemake --list-input-changes`

```

### How do I trigger re-runs for rules with updated code or parameters?

```bash
snakemake -n -R `snakemake --list-params-changes`
```

and

```bash
snakemake -n -R `snakemake --list-code-changes`
```


