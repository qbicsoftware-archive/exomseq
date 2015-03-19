configfile: "config.json"


if "groups" not in config:
    raise ValueError("Config must contain 'groups'")

if len(config['groups']) != 2:
    raise ValueError("Only two groups are supported")

def get_fastq_per_group(wildcards):
    groups = config['groups']
    return expand("{name}.fastq.gz", name=groups[wildcards['group']])

GROUP_NAMES = list(config['groups'].keys())

rule all:
    input: "normal_realigned.indexed.bam", "normal_qc_fastq.qcML"

rule merge:
    input: get_fastq_per_group
    output: "{group}_merged.fastq.gz"
    shell: "zcat {input} | gzip -1 > {output}"

rule read_qc:
    input: expand("{group}_merged.fastq.gz", group=GROUP_NAMES)
    output: "normal_qc_fastq.qcML"
    shell: "ReadQC -in1 {input[0]} -in2 {input[1]} -out {output}"

rule trim:
    input: expand("{group}_merged.fastq.gz", group=GROUP_NAMES)
    output: expand("{group}_trimmed.fastq.gz", group=GROUP_NAMES)
    params: adapter1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTTA", \
            adapter2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTC", \
            match_perc="80", \
            qcut="15"
    run:
        shell("SeqPurge -in1 {input[0]} -in2 {input[1]} "
              "-out1 {output[0]} -out2 {output[1]} "
              "-match_perc {match_perc} -qcut {qcut}")

rule map_bwa:
    input: expand("{group}_trimmed.fastq.gz", group=GROUP_NAMES)
    output: "normal_bwa.sam"
    params: bwa_options="-M -R '@RG\tID:normal\tSM:normal\tLB:normal'"
    run:
        shell("bwa mem {fasta} {bwa_options} {input} > {output}")

rule sam_to_bam:
    input: "{name}.sam"
    output: "{name}.bam"
    shell: "samtools view -Sb {input} -o {output}"

rule sam_sort:
    input: "{file}.bam"
    output: "{file}.sorted.bam"
    shell: "samtools sort {input} {output}"

rule sam_index:
    input: "{name}.bam"
    output: "{name}.indexed.bam"
    run:
        shell("ln {input} {output}")
        shell("samtools index {output}")

rule map_stampy:
    input: "normal_bwa.bam"
    output: "normal_stampy.bam"
    params: "--readgroup=ID:normal --bamkeepgoodreads"
    run:
        shell("stampy -g {config['fasta']} -h {config['fasta']} {params} -M "
              "{input} -o {output}")

rule duplicates:
    input: "normal_stampy.sorted.bam"
    output: bam="normal_duplicates.bam", matrix="normal_duplicates.matrix"
    run:
        shell("java -Xmx2g {config['picard_tools']}/MarkDuplicates.jar "
              "I={input} O={output['bam']} M={output['matrix']} AS=true")

rule left_align:
    input: "normal_duplicates.indexed.bam"
    output: "normal_leftaligned.bam"
    run:
        shell("BamLeftAlign -in {input} -out {output} -ref {config['fasta']}")

rule re_align_intervals:
    input: "normal_leftaligned.indexed.bam"
    output: "normal_realigned.intervals"
    run:
        shell("java -Xmx4g -jar "
              "{config['GATK']}/GenomeAnalysisTKLite.jar "
              "-T RealignerTargetCreator "
              "-I {input} "
              "-R {config['fasta']} "
              "-o {output} "
              "--intervals {config['intervals']}")

rule re_align:
    input: bam="normal_leftaligned.indexed.bam", \
           intervals="normal_realigned.intervals"
    output: "normal_realigned.bam"
    run:
        shell("java -Xmx4g -jar "
              "{config['GATK']}/GenomeAnalysisTKLite.jar "
              "-T IndelRealigned "
              "-I {input['bam']} "
              "-R {config['fasta']} "
              "--targetIntervals {input['intervals']} "
              "-o {output}")
