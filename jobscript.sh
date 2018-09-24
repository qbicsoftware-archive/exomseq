#!/bin/sh
#PBS -q cfc
#PBS -A qbic
#PBS -e ../logs/jobscript.{job.rule.name}.e$PBS_JOBID
#PBS -o ../logs/jobscript.{job.rule.name}.o$PBS_JOBID
#PBS -l nodes=1:ppn=10:cfc
#PBS -l walltime=30:00:00
#PBS -l mem=50g
# properties = {properties}

set -e

echo Running on machine $(hostname)

module load qbic/anaconda
module load devel/java_jdk/1.7.0u45
module load qbic/samtools
module load qbic/ngs-bits
module load qbic/bwa
module load qbic/stampy
module load qbic/picard/git
module load bio/gatk/3.3
module load qbic/annovar/0.1


{exec_job}
exit 0
