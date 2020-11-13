##fasta
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.gz

##bed
wget ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.current.txt
cat CCDS.current.txt | grep  "Public" | perl -alne '{/\[(.*?)\]/;next unless $1;$gene=$F[2];$exons=$1;$exons=~s/\s//g;$exons=~s/-/\t/g;print "$F[0]\t$_\t$gene" foreach split/,/,$exons;}'|sort -u |bedtools sort -i |awk '{if($3>$2) print "chr"$0"\t0\t+"}'  > hg38.exon.bed

##dict
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.dict

###step1################################
###bed to interval, use PreprocessIntervals to get interval lists
#PBS -N GATK_CNV_step1
#PBS -l nodes=2:ppn=16
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
GATK=/public/home/liuxs/anaconda3/envs/wes/share/gatk4-4.1.3.0-0/gatk-package-4.1.3.0-local.jar
ref=/public/home/liuxs/biodata/reference/genome/hg38/ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta
bed=/public/home/wangshx/wt/ref/hg38.exon.bed
dict=/public/home/liuxs/biodata/reference/genome/hg38/ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.dict

## bed to intervals_list
java -jar $GATK BedToIntervalList -I ${bed} -O /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/WES/cnv/step1/out/hg38.exon.interval_list -SD ${dict}

## Preprocess Intervals
java -jar $GATK PreprocessIntervals \
-L /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/WES/cnv/hg38.exon.interval_list \
--sequence-dictionary ${dict} \
--reference ${ref}  \
--padding 250 \
--bin-length 0 \
--interval-merging-rule OVERLAPPING_ONLY \
--output /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/WES/cnv/step1/out/targets.preprocessed.interval.list

###step2###############################
###CollectReadCounts counts reads that overlap the interval###########
#PBS -N GATK_CNV_step2
#PBS -l nodes=2:ppn=16
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

GATK=/public/home/liuxs/anaconda3/envs/wes/share/gatk4-4.1.3.0-0/gatk-package-4.1.3.0-local.jar
ref=/public/home/liuxs/biodata/reference/genome/hg38/ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta

sample=$(basename <bam> .marked.BQSR.bam)
java -jar -Xmx20G -Djava.io.tmpdir=./tmp $GATK CollectReadCounts \
-I /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/BQSR/<bam> \
-L /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/WES/cnv/step1/out/targets.preprocessed.interval.list \
--interval-merging-rule OVERLAPPING_ONLY \
-R ${ref} \
-O /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/WES/cnv/step2/out/${sample}.hdf5

##############step3##########################
######CreateReadCountPanelOfNormals
#!/bin/bash
GATK=/public/home/liuxs/anaconda3/envs/wes/share/gatk4-4.1.3.0-0/gatk-package-4.1.3.0-local.jar

java -jar -Xmx20G -Djava.io.tmpdir=./tmp $GATK CreateReadCountPanelOfNormals \
-I ../step2/out/SRR3083867.hdf5 \
-I ../step2/out/SRR3083840.hdf5 \
-I ../step2/out/SRR3083842.hdf5 \
-I ../step2/out/SRR3083838.hdf5 \
-I ../step2/out/SRR3083858.hdf5 \
-I ../step2/out/SRR3083864.hdf5 \
-I ../step2/out/SRR3083871.hdf5 \
-I ../step2/out/SRR4289714.hdf5 \
-I ../step2/out/SRR3083850.hdf5 \
-I ../step2/out/SRR4289716.hdf5 \
-I ../step2/out/SRR4289718.hdf5 \
-I ../step2/out/SRR3083883.hdf5 \
-I ../step2/out/SRR3083856.hdf5 \
-I ../step2/out/SRR3083869.hdf5 \
-I ../step2/out/SRR3083846.hdf5 \
-I ../step2/out/SRR4289720.hdf5 \
-I ../step2/out/SRR4289722.hdf5 \
-I ../step2/out/SRR3083860.hdf5 \
-I ../step2/out/SRR3083848.hdf5 \
-I ../step2/out/SRR3083873.hdf5 \
-I ../step2/out/SRR3083862.hdf5 \
-I ../step2/out/SRR3083879.hdf5 \
-I ../step2/out/SRR3083875.hdf5 \
-I ../step2/out/SRR3083854.hdf5 \
-I ../step2/out/SRR3083881.hdf5 \
-I ../step2/out/SRR3083844.hdf5 \
-I ../step2/out/SRR4289724.hdf5 \
-I ../step2/out/SRR3083852.hdf5 \
-I ../step2/out/SRR3083877.hdf5 \
-I ../step2/out/SRR4289727.hdf5 \
-I ../step2/out/SRR4289729.hdf5 \
-I ../step2/out/SRR4289731.hdf5 \
-I ../step2/out/SRR4289733.hdf5 \
-I ../step2/out/SRR4289735.hdf5 \
-I ../step2/out/SRR4289737.hdf5 \
-I ../step2/out/SRR4289739.hdf5 \
-I ../step2/out/SRR4289741.hdf5 \
-I ../step2/out/SRR4289743.hdf5 \
--minimum-interval-median-percentile 5.0 \
-O ./out/cnvponC.pon.hdf5

###########step4##################
########normalization and Denoise
#!/bin/bash
GATK=/public/home/liuxs/anaconda3/envs/wes/share/gatk4-4.1.3.0-0/gatk-package-4.1.3.0-local.jar

cat tumor_sample | while read id
do 
  java -jar -Xmx20G -Djava.io.tmpdir=./tmp $GATK DenoiseReadCounts \
  -I ../step2/out/${id}.hdf5 \
  --count-panel-of-normals ../step3/out/cnvponC.pon.hdf5 \
  --standardized-copy-ratios out/${id}.standardizedCR.tsv \
  --denoised-copy-ratios out/${id}.denoisedCR.tsv
done

###############step5##################
#PBS -N GATK_CNV_step5
#PBS -l nodes=1:ppn=8
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

GATK=/public/home/liuxs/anaconda3/envs/wes/share/gatk4-4.1.3.0-0/gatk-package-4.1.3.0-local.jar
ref=/public/home/liuxs/biodata/reference/genome/hg38/ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta

sample=$(basename <bam> .marked.BQSR.bam)
java -jar -Xmx20G -Djava.io.tmpdir=/public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/WES/cnv/step5/tmp $GATK CollectAllelicCounts \
-L /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/WES/cnv/step1/out/targets.preprocessed.interval.list \
-I /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/BQSR/<bam> \
-R ${ref} \
-O /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/WES/cnv/step5/out/${sample}.allelicCounts.tsv


#####step6###############
#####ModelSegments
#PBS -N GATK_CNV_step6
#PBS -l nodes=1:ppn=8
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

GATK=/public/home/liuxs/anaconda3/envs/wes/share/gatk4-4.1.3.0-0/gatk-package-4.1.3.0-local.jar

java -jar -Xmx20G -Djava.io.tmpdir=/public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/WES/cnv/step6/no_allelic_counts/tmp $GATK ModelSegments \
--denoised-copy-ratios /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/WES/cnv/step4/out/<file1>.denoisedCR.tsv \
--output /public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/sra/WES/cnv/step6/no_allelic_counts/out/ \
--output-prefix <file1>

############R
########do absolute
library(DoAbsolute)
library(data.table)
library(dplyr)

willy_seg <- readRDS("/public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/ABSOLUTE/willy_seg.rds")
willy_maf <- readRDS("/public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/ABSOLUTE/willy_maf.rds")
#/public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/absolute

DoAbsolute(Seg = willy_seg, Maf = willy_maf, platform = "Illumina_WES", copy.num.type = "total",
           results.dir = "/public/home/wangshx/wt/immune_therapy_dataset/willy_hugo_2016_cell/ABSOLUTE/output", nThread = 20, keepAllResult = TRUE, verbose = TRUE)



