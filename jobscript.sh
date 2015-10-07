#!/bin/sh
#PBS -q cfc
#PBS -A qbic
#PBS -l nodes=1:ppn=2:cfc
#PBS -l walltime=40:00:00
#PBS -e ../logs/jobscript.{job.rule.name}.e$PBS_JOBID
#PBS -o ../logs/jobscript.{job.rule.name}.o$PBS_JOBID
# properties = {properties}

set -e

module load bio/fastqc/0.10
module load qbic/anaconda
module load qbic/htseq/0.6.1p2
module load qbic/tophat
module load bio/samtools/1.2

{exec_job}
exit 0
