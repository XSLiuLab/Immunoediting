#PBS -N fiter_<sample>
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
/public/home/liuxs/anaconda3/envs/wes/share/gatk4-4.1.3.0-0/gatk-package-4.1.3.0-local.jar FilterMutectCalls \
 -R /public/home/liuxs/biodata/reference/genome/hg38/hg38.fa \
 -V /public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/vcf/unfilter_vcf/<sample>.vcf \
 -O /public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/vcf/filter_vcf/<sample>.filter.vcf \
 --tumor-segmentation /public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/vcf/pileups/<sample>-segments.table \
 --contamination-table /public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/vcf/pileups/<sample>-contamination.table \
 --ob-priors /public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/vcf/artifacts/<sample>-tumor-artifact-prior.tar.gz
