##########step1###########
#####sra to fastq
#PBS -N sra2fastq
#PBS -l nodes=1:ppn=16
#PBS -l walltime=24:00:00
#PBS -S /bin/bash
#PBS -q slst_pub
#PBS -j oe

if [ -f "/public/home/liuxs/anaconda3/etc/profile.d/conda.sh" ]; then
    . "/public/home/liuxs/anaconda3/etc/profile.d/conda.sh"
else
      export PATH="/public/home/liuxs/anaconda3/bin:$PATH"
fi
conda activate wes

fastq-dump --split-3 /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/RNA_seq/all_SRR/<srr> -O /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/RNA_seq/fastq/

#############step2############
#####quality control##########
####fastqc + remove adapt#########
#PBS -N fastqc
#PBS -l nodes=1:ppn=16
#PBS -l walltime=24:00:00
#PBS -S /bin/bash
#PBS -q slst_pub
#PBS -j oe

if [ -f "/public/home/liuxs/anaconda3/etc/profile.d/conda.sh" ]; then
    . "/public/home/liuxs/anaconda3/etc/profile.d/conda.sh"
else
      export PATH="/public/home/liuxs/anaconda3/bin:$PATH"
fi
conda activate wes

fastqc --outdir /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/RNA_seq/fastqc --threads 16 /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/RNA_seq/fastq/<fastq> >> /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/RNA_seq/fastqc/<fastq1>_fastqc.log 2>&1

#PBS -N multiqc
#PBS -l nodes=2:ppn=16
#PBS -l walltime=24:00:00
#PBS -S /bin/bash
#PBS -q slst_pub
#PBS -j oe

if [ -f "/public/home/wangshx/wt/miniconda3/etc/profile.d/conda.sh" ]; then
    . "/public/home/wangshx/wt/miniconda3/etc/profile.d/conda.sh"
 else
     export PATH="/public/home/wangshx/wt/miniconda3/bin:$PATH"
 fi

conda activate /public/home/wangshx/wt/miniconda3

multiqc /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/RNA_seq/fastqc/*zip -o /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/RNA_seq/multiqc

#PBS -N trim_galore
#PBS -l nodes=1:ppn=8
#PBS -l walltime=24:00:00
#PBS -S /bin/bash
#PBS -q slst_pub
#PBS -j oe
 
if [ -f "/public/home/wangshx/wt/miniconda3/etc/profile.d/conda.sh" ]; then
    . "/public/home/wangshx/wt/miniconda3/etc/profile.d/conda.sh"
 else
     export PATH="/public/home/wangshx/wt/miniconda3/bin:$PATH"
 fi

source activate /public/home/wangshx/wt/miniconda3

fq1=/public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/RNA_seq/fastq/<fastq1>.fastq
fq2=/public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/RNA_seq/fastq/<fastq2>.fastq
id=$(basename <fastq1> _1)

trim_galore --paired -q 20 --dont_gzip --phred33 --length 30 --stringency 3 -o /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/RNA_seq/clean_fastq/ $fq1 $fq2 >> /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/RNA_seq/clean_fastq/${id_trim}.log

################step3##############
########align+ sam to bam
#PBS -N hista2
#PBS -l nodes=1:ppn=7
#PBS -l walltime=24:00:00
#PBS -S /bin/bash
#PBS -q slst_pub
#PBS -j oe

name=$(basename <fastq1> _1_val_1.fq)

hisat2 --new-summary -p 6 -x /public/home/liuxs/biodata/reference/index/hisat/hg38/genome -1 /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/RNA_seq/clean_fastq/<fastq1> -2 /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/RNA_seq/clean_fastq/<fastq2> -S /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/RNA_seq/sam/${name}.sam

#PBS -N sam2bam
#PBS -l nodes=1:ppn=7
#PBS -l walltime=24:00:00
#PBS -S /bin/bash
#PBS -q slst_pub
#PBS -j oe

name=$(basename <sam> .sam)
samtools view -S /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/RNA_seq/sam/<sam> -b > /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/RNA_seq/bam/${name}.bam
samtools sort /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/RNA_seq/bam/${name}.bam -o /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/RNA_seq/sorted_bam/${name}_sorted.bam

####################step4#######################
###########calculate TPM
#PBS -N TPMCalculator
#PBS -l nodes=1:ppn=8
#PBS -l walltime=24:00:00
#PBS -S /bin/bash
#PBS -q slst_pub
#PBS -j oe
 
if [ -f "/public/home/wangshx/wt/miniconda3/etc/profile.d/conda.sh" ]; then
    . "/public/home/wangshx/wt/miniconda3/etc/profile.d/conda.sh"
 else
     export PATH="/public/home/wangshx/wt/miniconda3/bin:$PATH"
 fi

source activate /public/home/wangshx/wt/miniconda3

TPMCalculator -g /public/home/liuxs/biodata/reference/gtf/gencode/gencode.v29.annotation.gtf -b /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/RNA_seq/sorted_bam/<bam> -p -a -e 







