#PBS -N mu_call_<sample>
#PBS -l nodes=1:ppn=8
#PBS -l walltime=40:00:00
#PBS -S /bin/bash
#PBS -q normal_8
#PBS -j oe

if [ -f "/public/home/liuxs/anaconda3/etc/profile.d/conda.sh" ]; then
    . "/public/home/liuxs/anaconda3/etc/profile.d/conda.sh"
else
      export PATH="/public/home/liuxs/anaconda3/bin:$PATH"
fi

conda activate wes


java -jar -Xmx12G -Djava.io.tmpdir=/public/home/wangshx/wx/tmp \
/public/home/liuxs/anaconda3/envs/wes/share/gatk4-4.1.3.0-0/gatk-package-4.1.3.0-local.jar Mutect2 \
 -R /public/home/liuxs/biodata/reference/genome/hg38/hg38.fa \
 -I /public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/worked_bam/BQSR/<sample1>.marked.BQSR.bam \
 -tumor <sample1> \
 -I /public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/worked_bam/BQSR/<sample2>.marked.BQSR.bam \
 -normal <sample2> \
 --germline-resource /public/home/liuxs/biodata/reference/genome/hg38/mu_calling_need_file/af-only-gnomad.hg38.vcf.gz \
 --panel-of-normals /public/home/liuxs/biodata/reference/genome/hg38/mu_calling_need_file/1000g_pon.hg38.vcf.gz \
 --af-of-alleles-not-in-resource 0.0000025 \
 --disable-read-filter  MateOnSameContigOrNoMappedMateReadFilter \
 --bam-output /public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/vcf/unfilter_vcf/<sample>.bam \
 --f1r2-tar-gz /public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/vcf/unfilter_vcf/<sample>.tar.gz \
 -O /public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/vcf/unfilter_vcf/<sample>.vcf
