from os.path import join as pjoin
from os.path import exists as pexists
import pandas as pd
import os
import hashlib
import subprocess
import json

configfile: "config.json"

workdir: config["var"]

DATA = config['data']
RESULT = config['result']
LOGS = config['logs']
REF = config['ref']
INI_PATH = config['etc']
SNAKEDIR = config['src']

def data(path):
    return os.path.join(DATA, path)

def ref(path):
    return os.path.join(REF, path)

def log(path):
    return os.path.join(LOGS, path)

def result(path):
    return os.path.join(RESULT, path)

def etc(path):
    return os.path.join(INI_PATH, path)

try:
    with open(etc("params.json")) as f:
        parameters = json.load(f)


except OSError as e:
    print("Could not read parameter file: " + str(e), file=sys.stderr)
    sys.exit(1)
except ValueError as e:
    print("Invalid parameter file: " + str(e), file=sys.stderr)
    sys.exit(1)

default_params = {}
default_params.update(parameters)
parameters = default_params

parameters['indexed_genome'] = ref(parameters['indexed_genome'])
indexed_genome = parameters["indexed_genome"]
if not os.path.exists(indexed_genome + '.fa'):
    raise ValueError("Could not find indexed genome file %s" % indexed_genome)



def read_design():
    df = pd.read_csv(etc('design.csv'), sep='\t')

    test_samples = df[df['SAMPLE TYPE'] == 'Q_TEST_SAMPLE']
    #dna_samples = test_samples[test_samples['Q_SAMPLE_TYPE'] == 'DNA']

    pools = {}
    pool = {}
    pools['all'] = pool

    files = os.listdir(DATA)

    for barcode in test_samples.Identifier:
        sample_files = [file.split('.')[0] for file in files if (file.startswith(barcode) & file.endswith('.gz'))]
        if len(sample_files) == 2:
            pool[barcode] = sample_files
    return pools

pools = read_design()

print(pools)

def fastq_per_group(wildcards):
    pool = wildcards['pool']
    group = wildcards['group']
    print(group)
    if pool not in pools:
        return ['/i_do_not_exist']
    if group not in pools[pool]:
        return ['/i_do_not_exist']
    names = pools[pool][group]
    return expand(data("{name}.fastq.gz"), name=names)

all_files = []
for pool in pools.values():
    for group in pool.values():
        all_files.extend(group)


def groups_per_pool(format_string):
    def inner(wildcards):
        pool = wildcards['pool']
        if pool not in parameters['pools']:
            return ["/i_do_not_exist"]
        groups = parameters['pools'][pool]
        group_names = groups.keys()
        return expand(format_string, pool=[pool], group=group_names)
    return inner


def expand_pools_groups(format_string):
    files = []
    for pool in pools:
        for group in pools[pool]:
            files.append(expand(format_string, pool=pool, group=group))
    return files


def touch(fname):
    if os.path.exists(fname):
        os.utime(fname, None)
    else:
        open(fname, 'a').close()

ID="merged"


rule all:
    input:
        expand("checksums/{name}.txt", name=all_files),
        #expand("mapqc/{pool}.qcML", pool=pools.keys()),
        expand_pools_groups(result("{pool}_{group}_readqc.qcML")),
        expand("fastqc/{name}", name=all_files),
        expand_pools_groups("trim_fastqc/{pool}_{group}"),
        #expand_pools_groups("stampy/{pool}_{group}.sorted.bam.flagstat"),
        expand_pools_groups("map_bwa/{pool}_{group}.sorted.bam.flagstat"),
        expand_pools_groups("duplicates/{pool}_{group}.sorted.bam"),
        expand_pools_groups("clipped/{pool}_{group}.sorted.bam.flagstat"),
        expand_pools_groups("mapqc/{pool}_{group}.qcML"),
        expand_pools_groups("fastqc_clipped/{pool}_{group}"),
        #expand("fastqc_clipped/{pool}", pool=pools.keys())
        #expand("map_fastqc/{pool}.qcML", pool=pools.keys())
        #expand("result/{name}.qcML", name=all_files)
        expand("merged_bam/{id}_sorted.bam",id=ID),
        expand("merged_bam/{id}_sorted.bam.bai",id=ID),
        expand("annovar/annovar_in/{id}_sorted.vcf.avinput",id=ID),
        expand("annovar/annovar_out/{id}_sorted.vcf.avoutput",id=ID),
        expand("snpeff/{id}_sorted.snpeff.vcf",id=ID)



rule verify_checksums:
    input: filein=data("{name}.fastq.gz"),\
           cs=data("checksums/{name}.fastq.gz.sha256sum")
    output: "checksums/{name}.txt"
    run:
        try:
            sha=hashlib.sha256(open(str(input.filein), 'rb').read()).hexdigest()
            fo = open(str(output), "w")
            fo.write("%s"% sha)
            orig_sha = open(str(input.cs), 'r').read().split(" ")[0]
            fo.write("\n%s"% orig_sha)
            fo.close()
            if sha == orig_sha:
                touch(str(output)+"OK")
            else:
                touch(str(output)+"FAILED")
        except Exception:
            pass


rule flagstat:
    input: "{file}.sorted.bam"
    output: "{file}.sorted.bam.flagstat"
    shell: "samtools flagstat {input} > {output}"


rule read_qc:
    input: fastq_per_group
    output: result("{pool}_{group}_readqc.qcML")
    # TODO merge if != 2 elements per group / pairwise?
    shell: "ReadQC -in1 {input[0]} -in2 {input[1]} -out {output}"


rule fastqc:
    input: data("{name}.fastq.gz")
    output: "fastqc/{name}"
    run:
        try:
            os.mkdir(str(output))
        except Exception:
            pass
        shell("/lustre_cfc/software/qbic/FastQC/fastqc {input} -o {output}")


rule trim_fastqc:
    input: "trim/{pool}/{group}"
    output: "trim_fastqc/{pool}_{group}"
    run:
        try:
            os.mkdir(str(output))
            os.mkdir(os.path.join(str(output), 'forward'))
            os.mkdir(os.path.join(str(output), 'backward'))
        except Exception:
            pass
        shell("/lustre_cfc/software/qbic/FastQC/fastqc -o {output}/forward {input}/forward.fastq.gz")
        shell("/lustre_cfc/software/qbic/FastQC/fastqc -o {output}/backward {input}/backward.fastq.gz")


rule trim:
    input: fastq_per_group
    output: "trim/{pool}/{group}"
    params: match_perc="80", qcut="15"
    run:
        # TODO more than two elements per group
        # TODO get adapters from fastcq
        # TODO Default adapter is used
        try:
            os.mkdir(str(output))
        except Exception:
            pass
        shell("SeqPurge -in1 {input[0]} -in2 {input[1]} -out1 {output}/forward.fastq.gz -out2 {output}/backward.fastq.gz -match_perc 80 -qcut 15")


rule map_bwa:
    input: "trim/{pool}/{group}"
    output: "map_bwa/{pool}_{group}.unsorted.sam"
    threads: 10
    run:
        # TODO additional header LB (library)?
        options = r"-M -R '@RG\tID:{group}\tSM:{group}'".format(**wildcards)
        options = options + " " + indexed_genome
        shell("bwa mem -t {threads} %s {input}/forward.fastq.gz {input}/backward.fastq.gz > {output}" % options)


rule sam_to_bam:
    input: "{name}.unsorted.sam"
    output: "{name}.unsorted.bam"
    shell: "samtools view -Sb {input} -o {output}"


rule sam_sort:
    input: "{file}.unsorted.bam"
    output: "{file}.sorted.bam"
    shell: "samtools sort -O bam -T $(mktemp) -o {output} {input}"


rule sam_index:
    input: "{file}.sorted.bam"
    output: "{file}.sorted.bam.bai"
    run:
        shell("samtools index {input}")

rule sam_index_unsorted:
    input: "{file}.unsorted.bam"
    output: "{file}.unsorted.bam.bai"
    run:
        shell("samtools index {input}")



rule map_stampy:
    input: "map_bwa/{pool}_{group}.unsorted.bam"
    output: "stampy/{pool}_{group}.unsorted.sam"
    run:
        shell("stampy.py -g {fasta} -h {fasta} "
              "--readgroup=ID:{group} --bamkeepgoodreads -M "
              "{input} -o {output}".format(fasta=indexed_genome,
                                           input=input,
                                           output=output,
                                           **wildcards))


rule merge_groups:
    #input: groups_per_pool("stampy/{pool}_{group}.unsorted.sam")
    input: "map_bwa/{pool}_{group}.unsorted.bam"
    output: "merged/.unsorted.bam"
    run:
        # TODO len(group) != 2
        shell("java -Xmx2g {config['picard_tools']}/MergeSamFiles.jar "
              "I={input[0]} I={input[1]} O={output} AS=true")


rule duplicates:
    input: "map_bwa/{pool}_{group}.sorted.bam"
    output: bam="duplicates/{pool}_{group}.sorted.bam", \
            matrix="duplicates/{pool}_{group}.matrix"
    run:
        shell("picard MarkDuplicates "
              "I={input} O={output.bam} M={output.matrix} AS=true")


rule left_align:
    input: bam="duplicates/{name}.sorted.bam",\
           bai="duplicates/{name}.sorted.bam.bai"
    output: "left_align/{name}.unsorted.bam"
    run:
        shell("BamLeftAlign -in {input.bam} -out {output} -ref %s.fa"
               % indexed_genome)


rule re_align_intervals:
    input: bam="left_align/{name}.unsorted.bam", \
           bai="left_align/{name}.unsorted.bam.bai"
    output: "left_align/{name}.intervals"
    run:
        shell(("java -jar %s/GenomeAnalysisTK.jar "  # TODO Lite
               "-T RealignerTargetCreator "
               "-I {input.bam} "
               "-R %s.fa "
               "-o {output}")
              % (os.environ['GATK_JARS'], indexed_genome))
        #"--intervals {config['intervals']}")


rule re_align:
    input: bam="left_align/{name}.unsorted.bam", \
           bai="left_align/{name}.unsorted.bam.bai", \
           intervals="left_align/{name}.intervals"
    output: "realigned/{name}.unsorted.bam"
    run:
        shell(("java -jar %s/GenomeAnalysisTK.jar "  # TODO Lite
               "-T IndelRealigner "
               "-I {input.bam} "
               "-R %s.fa "
               "--targetIntervals {input.intervals} "
               "-o {output}")
              % (os.environ['GATK_JARS'], indexed_genome))


rule clip_overlap:
    input: bam="realigned/{name}.unsorted.bam", \
           bai="realigned/{name}.unsorted.bam.bai"
    output: "clipped/{name}.unsorted.bam"
    run:
        shell("BamClipOverlap -in {input.bam} -out {output}")


rule mapqc:
    input: bam="clipped/{name}.sorted.bam",\
           bai="clipped/{name}.sorted.bam.bai"
    output: "mapqc/{name}.qcML"
    shell: "MappingQC -in {input.bam} -out {output} -wgs hg19"

rule fastqc_clipped:
    input: "clipped/{pool}.sorted.bam"
    output: "fastqc_clipped/{pool}"
    run:
        try:
            os.mkdir(str(output))
        except Exception:
            pass
        shell("/lustre_cfc/software/qbic/FastQC/fastqc {input} -o {output}")



##create textfile of bam files to merge
rule create_txtfile:
    input: expand_pools_groups("clipped/{pool}_{group}.sorted.bam")
    output: "filenames.txt"
    run:
        shell("echo {input} | tr ' '  '\n' > {output}")


## merge bam files
rule merge_bam:
    input: "filenames.txt"
    output: "merged_bam/merged_sorted.bam"
    run:
        shell("samtools merge {output} -b {input}")


##index merged bam file
rule merge_sam_index:
    input:  "merged_bam/merged_sorted.bam"
    output: "merged_bam/merged_sorted.bam.bai"
    run:
        shell("samtools index {input}")


## create basic variant calls ,here with freebayes
rule freebayes:
    input: bam="merged_bam/merged_sorted.bam", \
           bai="merged_bam/merged_sorted.bam.bai"
    output: "vcf_base/merged_base.vcf"
    run:
        shell("/lustre_cfc/software/qbic/freebayes/bin/freebayes --min-base-quality 20 --min-alternate-qsum 90 -f %s.fa -b {input.bam} > {output}" % indexed_genome)


## filter variants according to variant quality>5 , alternate observations>=2, vcflib is needed and was installed with git clone --recursive git://github.com/ekg/vcflib.git, cd vcflib and make under qbic/software...
rule filter_variants:
    input: "vcf_base/merged_base.vcf"
    output: "vcf_filtered/merged_filtered.vcf"
    run:
        shell("/lustre_cfc/software/qbic/vcflib/bin/vcffilter -f 'QUAL > 5 & AO > 2' {input} > {output}")


## split complex variants to primitives
## this step has to be performed before vcfbreak multi - otherwise mulitallelic variants that contain both 'hom' and 'het' genotypes fail
rule split_complex_variants:
    input: "vcf_filtered/merged_filtered.vcf"
    output: "vcf_prim/merged_prim.vcf"
    run:
        shell("/lustre_cfc/software/qbic/vcflib/bin/vcfallelicprimitives -kg {input} > {output}")


## as the error ocurred Error : [E::bcf_calc_ac] Incorrect AN/AC counts at chr1:85062524 the next rule vcffixup had to be added. Why?:
## To follow up on this, the issue is due to having an improper AN allele output coming from splitting up FreeBayes MNP calls.
## It looks like the problem is that vcfallelicprimitives, which splits MNPs from FreeBayes output, creates invalid AN fields.
## These should be a single integer but in cases where you have multiple alleles and a MNP, it produces an allele based list:
## solution is to use vcffixup to correctly resolve it back to the right type:
rule vcffixup:
    input: "vcf_prim/merged_prim.vcf"
    output: "vcf_prim/merged_prim_fixed.vcf"
    run:
        shell("/lustre_cfc/software/qbic/vcflib/bin/vcffixup {input} > {output}")


## split multi-allelic variants: If multiple alleles are specified in a single record, break the record into multiple lines, preserving allele-specific INFO fields.
rule split_multiallelic_variants:
    input: "vcf_prim/merged_prim_fixed.vcf"
    output: "vcf_split/merged_split.vcf"
    run:
        shell("/lustre_cfc/software/qbic/vcflib/bin/vcfbreakmulti {input} > {output}")


## align INDELs to the left. Note that complex indels and multi-allelic deletions are not shifted!
#rule re_align_INDELS:
#    input: "vcf_split/merged_split.vcf"
#    output: "vcf_aligned/merged_aligned.vcf"
#    run:
#        shell("VcfLeftAlign -in {input} -out {output} -ref %s.fa"
#               % ref(config['params']['fasta']))

## align INDELs to the left. Note that complex indels and multi-allelic deletions are not shifted!
# tool VcfLeftAlign did not work thus I tried bcftools norm:
rule re_align_INDELS:
    input: "vcf_split/merged_split.vcf"
    output: "vcf_aligned/merged_aligned.vcf"
    run:
        shell("/lustre_cfc/software/qbic/bcftools-1.2/bin/bcftools norm -f %s.fa {input} -o {output}"
                % indexed_genome)


## sort variants by genomic position
#rule sort_variants:
#    input: "vcf_aligned/merged_aligned.vcf"
#    output: "vcf_sorted/merged_sorted.vcf"
#    run:
#        shell("VcfSort -in {input} -out {output}")

## sort variants by genomic position. as VcfSort did not work and lost the sample info over multiple files I decided to switch to vcf-sort (no tool in the faster bcftool kit was available for that:
rule sort_variants:
    input: "vcf_aligned/merged_aligned.vcf"
    output: "vcf_sorted/merged_sorted.vcf"
    run:
        shell("cat {input} | /lustre_cfc/software/qbic/vcftools_0.1.13/bin/vcf-sort > {output}")


#####
## some stats on base calling, filtered calling and final output which is the coordinate sorted vcf file
#####





### annotation

## 1) annovar annotation
## module is installed by Chris
## in /lustre_cfc/qbic/reference_genomes/ I downloaded database for annovar via annotate_variation.pl -buildver mm10 -downdb refGene mm10db_annovar/
## then check the log file in mm10db_annovar/
## #NOTICE: the FASTA file http://www.openbioinformatics.org/annovar/download mm10_refGeneMrna.fa.gz is not available to download but can be generated by the ANNOVAR software. PLEASE RUN THE FOLLOWING TWO COMMANDS CONSECUTIVELY TO GENERATE THE FASTA FILES:
## annotate_variation.pl --buildver mm10 --downdb seq mm10db_annovar/mm10_seq
## retrieve_seq_from_fasta.pl mm10db_annovar/mm10_refGene.txt -seqdir mm10db_annovar/mm10_seq -format refGene -outfile mm10db_annovar/mm10_refGeneMrna.fa
## check under /lustre_cfc/qbic/reference_genomes/mm10db_annovar/ how it eventually looks like

rule annovar_convert:
    input: "vcf_sorted/merged_sorted.vcf"
    output: "annovar/annovar_in/merged_sorted.vcf.avinput"
    run:
        try:
            shell("convert2annovar.pl -format vcf4 --allsample --withfreq --includeinfo {input} --outfile {output}")
        except subprocess.CalledProcessError:
            pass
        shell("convert2annovar.pl -format vcf4old --withfreq --includeinfo {input} --outfile {output}")
    #### the option --allsample is supported only if --format is 'vcf4'


rule annovar_annotate:
    input: "annovar/annovar_in/merged_sorted.vcf.avinput"
    output: "annovar/annovar_out/merged_sorted.vcf.avoutput"
    run:
        shell("annotate_variation.pl -buildver hg19 --otherinfo --infoasscore {input} > {output} --outfile {output} /lustre_cfc/qbic/reference_genomes/humandb_annovar/")




## 2) snpeff annotation
## installed snpeff under /qbic/software but not as module yet
## to check what databases are available you can do:
## java -Xmx4g -jar /lustre_cfc/software/qbic/snpEff/snpEff.jar databases | less
## #to find mouse specific stuff:
## java -Xmx4g -jar /lustre_cfc/software/qbic/snpEff/snpEff.jar databases | grep -i musculus
## Note: If you are running SnpEff from a directory different than the one it was installed, you will have to specify where the config file is. This is done using the '-c' command line option:
## in snpeff.config file changed data.dir: data.dir = /lustre_cfc/qbic/reference_genomes/snpeff
## run example
## java -Xmx4g -jar /lustre_cfc/software/qbic/snpEff/snpEff.jar eff -c /lustre_cfc/software/qbic/snpEff/snpEff.config -noStats -noLog -v GRCm38.79 sample_sorted.vcf > sample_sorted.snpeff.ann.vcf
## as pointed out in the command above direct to config with -c, set a version of the genome (should be found within databases available, see above command)
## also make sure a directory exists in data.dir named here GRCm38.79, then the db will be stored there.
## from Marc Sturm's commands leave out -spliceRegionIntronMax as not clear how to set this here
## #updated also the config file to add a line for Mitochondrial DNA for my genome: GRCm38.79.MT.codonTable : Vertebrate_Mitochondrial

rule snpeff:
    input: vcf="vcf_sorted/merged_sorted.vcf",\
           file="../etc/header_add.1.txt"
    output: "snpeff/merged_sorted.snpeff.vcf"
    run:
        shell("java -Xmx4g -jar /lustre_cfc/software/qbic/snpEff/snpEff.jar eff -c /lustre_cfc/software/qbic/snpEff/snpEff.config -v -cancer -cancerSamples {input.file} -noLog hg19 {input.vcf} > {output}")
