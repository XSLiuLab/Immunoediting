#PBS -N mu_call_<sample>
#PBS -l nodes=3:ppn=8
#PBS -l walltime=50:00:00
#PBS -S /bin/bash
#PBS -q normal_8
#PBS -j oe

if [ -f "/public/home/liuxs/anaconda3/etc/profile.d/conda.sh" ]; then
    . "/public/home/liuxs/anaconda3/etc/profile.d/conda.sh"
else
      export PATH="/public/home/liuxs/anaconda3/bin:$PATH"
fi

conda activate wes


java -jar -Xmx12G -Djava.io.tmpdir=/public/home/wangshx/wx/tmp /public/home/liuxs/anaconda3/envs/wes/share/gatk4-4.1.3.0-0/gatk-package-4.1.3.0-local.jar GetPileupSummaries \
 -I /public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/worked_bam/BQSR/<sample1>.marked.BQSR.bam \
 -L /public/home/liuxs/biodata/reference/genome/hg38/mu_calling_need_file/whole_exome_illumina_coding_v1.Homo_sapiens_assembly38.targets.interval_list \
 -V /public/home/liuxs/biodata/reference/genome/hg38/mu_calling_need_file/af-only-gnomad.hg38.vcf.gz \
 -O /public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/vcf/pileups/<sample1>-tumor.pileups.table




java -jar -Xmx12G -Djava.io.tmpdir=/public/home/wangshx/wx/tmp /public/home/liuxs/anaconda3/envs/wes/share/gatk4-4.1.3.0-0/gatk-package-4.1.3.0-local.jar GetPileupSummaries \
 -I /public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/worked_bam/BQSR/<sample2>.marked.BQSR.bam \
 -L /public/home/liuxs/biodata/reference/genome/hg38/mu_calling_need_file/whole_exome_illumina_coding_v1.Homo_sapiens_assembly38.targets.interval_list \
 -V /public/home/liuxs/biodata/reference/genome/hg38/mu_calling_need_file/af-only-gnomad.hg38.vcf.gz \
 -O /public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/vcf/pileups/<sample2>-normal.pileups.table




java -jar -Xmx12G -Djava.io.tmpdir=/public/home/wangshx/wx/tmp /public/home/liuxs/anaconda3/envs/wes/share/gatk4-4.1.3.0-0/gatk-package-4.1.3.0-local.jar CalculateContamination \
 -I /public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/vcf/pileups/<sample1>-tumor.pileups.table \
 -matched /public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/vcf/pileups/<sample2>-normal.pileups.table \
 -O /public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/vcf/pileups/<sample>-contamination.table \
 -tumor-segmentation /public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/vcf/pileups/<sample>-segments.table
