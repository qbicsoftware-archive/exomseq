configfile: "config.json"

params = config['params']

def fastq_per_group(wildcards):
    pool = wildcards['pool']
    group = wildcards['group']
    if pool not in params['pools']:
        return ['/i_do_not_exist']
    if group not in params['pools'][pool]:
        return ['/i_do_not_exist']
    names = params['pools'][pool][group]
    return expand("{name}.fastq.gz", name=names)


all_files = []
for pool in params['pools'].values():
    for group in pool.values():
        all_files.extend(group)


def groups_per_pool(format_string):
    def inner(wildcards):
        pool = wildcards['pool']
        if pool not in params['pools']:
            return ["/i_do_not_exist"]
        groups = params['pools'][pool]
        group_names = groups.keys()
        return expand(format_string, pool=[pool], group=group_names)
    return inner


rule all:
    input: expand("mapqc/{pool}.qcML", pool=params['pools'].keys()), \
           expand("map_fastqc/{pool}.qcML", pool=params['pools'].keys()), \
           expand("result/{name}.qcML", name=all_files)


rule read_qc:
    input: fastq_per_group
    output: "{group}_fastq.qcML"
    # TODO merge if != 2 elements per group / pairwise?
    shell: "ReadQC -in1 {input[0]} -in2 {input[1]} -out {output}"


rule fastqc:
    input: "{name}.fastq.gz"
    output: "result/{name}.qcML"
    shell: "fastqc {input} -out {output}"  # TODO check params


rule trim:
    input: fastq_per_group
    output: "trim/{pool}/{group}"
    params: match_perc="80", qcut="15"
    run:
        # TODO more than two elements per group
        # TODO get adapters from fastcq
        shell("SeqPurge -in1 {input[0]} -in2 {input[1]} "
              "-out1 {output[0]} -out2 {output[1]} "
              "-match_perc {match_perc} -qcut {qcut} "
              "-a1 {config[adapter[input[0]]]}")


rule map_bwa:
    input: "trim/{pool}/{group}"
    output: "map_bwa/{pool}_{group}.unsorted.sam"
    params: bwa_options="-M -R '@RG\tID:normal\tSM:normal\tLB:normal'"
    run:
        # TODO additional header LB (library)?
        read_group = "@RG\tID:{group}\tSM:{pool}"
        shell("bwa mem {fasta} {bwa_options} {input} > {output}")


rule sam_to_bam:
    input: "{name}.sam"
    output: "{name}.bam"
    shell: ""  # TODO


rule sam_sort:
    input: "{file}.unsorted.bam"
    output: "{file}.sorted.bam"
    shell: "samtools sort {input} {output}"


rule sam_index:
    input: "{name}.sorted.bam"
    output: "{name}.indexed.bam"
    run:
        shell("ln {input} {output}")
        shell("samtools index {output}")


rule map_stampy:
    input: "map_bwa/{pool}_{group}.indexed.bam"
    output: "stampy/{pool}_{group}.unsorted.sam"
    run:
        shell("stampy -g {config['fasta']} -h {config['fasta']} {params} "
              "--readgroup=ID:{group} --bamkeepgoodreads -M "
              "{input} -o {output}")


rule merge_groups:
    input: groups_per_pool("stampy/{pool}_{group}.unsorted.sam")
    output: "merged/{pool}.unsorted.bam"
    run:
        # TODO len(group) != 2
        shell("java -Xmx2g {config['picard_tools']}/MergeSamFiles.jar "
              "I={input[0]} I={input[1]} O={output} AS=true")


rule duplicates:
    input: "merged/{pool}.indexed.bam"
    output: bam="duplicates/{pool}.unsorted.bam", \
            matrix="duplicates/{pool}.matrix"
    run:
        shell("java -Xmx2g {config['picard_tools']}/MarkDuplicates.jar "
              "I={input} O={output['bam']} M={output['matrix']} AS=true")


rule left_align:
    input: "duplicates/{pool}.indexed.bam"
    output: "left_align/{pool}.unsorted.bam"
    run:
        shell("BamLeftAlign -in {input} -out {output} -ref {config['fasta']}")


rule re_align_intervals:
    input: "left_align/{pool}.indexed.bam"
    output: "left_align/{pool}.intervals"
    run:
        shell("java -Xmx4g -jar "
              "{config['GATK']}/GenomeAnalysisTKLite.jar "
              "-T RealignerTargetCreator "
              "-I {input} "
              "-R {config['fasta']} "
              "-o {output}")
        #"--intervals {config['intervals']}")


rule re_align:
    input: bam="left_align/{pool}.indexed.bam", \
           intervals="left_align/{pool}.intervals"
    output: "realigned/{pool}.unsorted.bam"
    run:
        shell("java -Xmx4g -jar "
              "{config['GATK']}/GenomeAnalysisTKLite.jar "
              "-T IndelRealigner "
              "-I {input['bam']} "
              "-R {config['fasta']} "
              "--targetIntervals {input['intervals']} "
              "-o {output}")


rule clip_overlap:
    input: "realigned/{pool}.indexed.bam"
    output: "clipped/{pool}.unsorted.bam"
    run:
        shell("BamClipOverlap -in {input} -out {output}")


rule mapqc:
    input: "clipped/{pool}.indexed.bam"
    output: "mapqc/{pool}.qcML"
    shell: "MappingQC -in {input} -out {output} -wgs hg19"


rule map_fastqc:
    input: "clipped/{pool}.indexed.bam"
    output: "map_fastqc/{pool}.qcML"
    shell: "fastqc {input}"  # TODO param for output dir?
