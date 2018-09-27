# Legacy Exom sequencing workflow at QBiC written in Snakemake 
## Main Tools used

samtools Version 1.3  
bwa Version 0.7.10-r789  
Picard Version: 1.137  
tools from ngs-bits (https://github.com/imgag/ngs-bits) such as  
  SeqPurge Version (0.1-4-gaed0c94)
  BamClipOverlap Version (0.1-4-gaed0c94)  
Freebayes Version: v0.9.21-19-gc003c1e  
...

## Description how to start
The workflow can be downloaded and run on a cluster environment using the tool qproject provided on github: https://github.com/qbicsoftware/qproject.  
This workflow does QC, mapping and variant calling. It produces a multi-file vcf file (output of Freebayes) as final table that can be annotated in a next step with any tool such as Annovar or Snpeff.  
The workflow uses a module system to load the required software. Be sure to check the jobscript.sh file to see which software modules are required. The modules are loaded automatically when using qproject run to start the workflow, otherwise they have to be loaded manually.

1) One should use qproject to download the files. This also creates all folders necessary for the workflow.

```
qproject create -t . -w github:qbicsoftware/exomseq
```

2) Be sure to add a config file "params.json" in `etc` which should look like this:

```
{
    "indexed_genome": "/<folder_to_index>/BWAIndex/hg19/hg19"
}
```
where `indexed_genome` gives the full path to the genome fasta file. Here in this example in the folder /<folder_to_index>/BWAIndex/hg19/ live all files from the BWA index (created with the command: bwa index -p hg19 -a bwtsw hg19.fa) plus the genome fasta file hg19.fa. The basename of the files must be identical such as here hg19.fa and hg19.amb etc.  
You also need a samtools faidx created index file here (such as hg19.fa.fai) and a picard dict file (such as hg19.dict).   


3) Input Fastq files  
The input Fastq files (not fastq.gz, unzip them before) should be copied in the folder `data`.   
Correct filenames that would work are `sample1_R1.fastq` and `sample1_R2.fastq`, but not `sample1_L001_R1.fastq` and `sample1_L001_R2.fastq`. That means you want to avoid '_' anywhere else than before the R1 and R2 of the filename.


4) File grouping
The input files (see point 3)) need to be grouped in such a way that the workflow knows which forward (R1) and reverse (R2) file belongs together. This is done with a file named "design.csv" (content however, should be tab separated) that should be generated in `etc` that should look like this:

```
Identifier	SAMPLE TYPE
sample1	Q_TEST_SAMPLE
sample2	Q_TEST_SAMPLE
```
The column Identifier is used to group the fastq files in `data` into sample groups (here sample1 and sample2). The column SAMPLE TYPE is just set to a string (here Q_TEST_SAMPLE) to select the right


To run the workflow navigate to the `src` folder.
Using `snakemake -n` one can display the operations the workflow will perform.
Using the `--dag` parameter and piping it to `dot` one can create a .pdf version of the directed acyclic graph used by snakemake to inspect the behavious of the workflow on a local machine.

```
module load qbic/anaconda
cd src/
snakemake -n
snakemake --dag | dot -Tpdf > dag.pdf
```

To run the workflow:

```
qproject run -t ..
```

While running one can inspect the log files (e.g. in a different screen session) for the progress and errors generated by the workflow:

```
cd logs/
tail snake.err -f
```

And to check the jobs on the computing cluster one can use `qstat`.

Alternatively to using `qproject run` one could use `snakemake -j` to run the workflow, but then be sure to check the `jobscript.sh` to load the required modules manually and also note that this would also not use `qsub` to submit the jobs.
