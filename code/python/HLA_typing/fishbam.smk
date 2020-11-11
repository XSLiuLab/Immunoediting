rule all:
    input:
        fa = "path/to/fished/fastq"

rule clean_bam:
    input: 
        bam = "/path/to/input/bam"
    output: 
        cleaned_bam = temp("/path/to/clean/bam")
    conda:
        "HLA.yaml"
    params:
        samtools = "view -b -F 4"
    shell: 
        "samtools {params.samtools} {input.bam} > {output.cleaned_bam}"

rule bam2fq:
    input:
        bam = "/path/to/clean/bam"
    output:
        fa = temp("/path/to/clean/fastq")
    conda:
        "HLA.yaml"
    shell:
        "samtools bam2fq {input.bam} > {output.fa}"

rule fished_bam:
    input: 
        fa = "/path/to/clean/fastq",
        ref = "/public/home/liuxs/anaconda3/pkgs/optitype-1.3.2-py27_3/share/optitype-1.3.2-3/data/hla_reference_dna.fasta"
    output: 
        bam = temp("path/to/fished_bam")
    shell: 
        "/public/home/liuxs/anaconda3/envs/python3/bin/razers3 -tc 6 -i 95 -m 1 -dr 0 -o {output} {input.ref} {input.fa}"

rule bam2fq_fished:
    input: 
        bam = "path/to/fished_bam"
    output: 
        fa = "path/to/fished/fastq"
    shell: 
        "samtools bam2fq {input.bam} > {output.fa}"
