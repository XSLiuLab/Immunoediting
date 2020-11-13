#PBS -N BQSR_<sample>
#PBS -l nodes=2:ppn=8
#PBS -l walltime=200:00:00
#PBS -S /bin/bash
#PBS -q normal_8
#PBS -j oe


if [ -f "/public/home/liuxs/anaconda3/etc/profile.d/conda.sh" ]; then
    . "/public/home/liuxs/anaconda3/etc/profile.d/conda.sh"
else
      export PATH="/public/home/liuxs/anaconda3/bin:$PATH"
fi


conda activate wes

dbsnp=/public/home/liuxs/biodata/reference/genome/hg38/ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz

dbsnp1000=/public/home/liuxs/biodata/reference/genome/hg38/ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

dbsnp1000G=/public/home/liuxs/biodata/reference/genome/hg38/ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz



java -jar -Xmx12G -Djava.io.tmpdir=/public/home/wangshx/wx/tmp  /public/home/liuxs/anaconda3/envs/wes/share/gatk4-4.1.3.0-0/gatk-package-4.1.3.0-local.jar BaseRecalibrator \
--reference /public/home/liuxs/biodata/reference/genome/hg38/hg38.fa \
--input /public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/worked_bam/mark_bam/<sample>.rmdup.bam \
--known-sites $dbsnp --known-sites $dbsnp1000 --known-sites $dbsnp1000G \
--output /public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/worked_bam/BQSR/<sample>_recal_data.table


java -jar -Xmx12G -Djava.io.tmpdir=/public/home/wangshx/wx/tmp  /public/home/liuxs/anaconda3/envs/wes/share/gatk4-4.1.3.0-0/gatk-package-4.1.3.0-local.jar ApplyBQSR \
-R /public/home/liuxs/biodata/reference/genome/hg38/hg38.fa \
-I /public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/worked_bam/mark_bam/<sample>.rmdup.bam \
--bqsr-recal-file /public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/worked_bam/BQSR/<sample>_recal_data.table \
-O /public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/worked_bam/BQSR/<sample>.marked.BQSR.bam
