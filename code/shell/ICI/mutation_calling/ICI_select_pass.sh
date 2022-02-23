#PBS -N filter_<sample>
#PBS -l nodes=1:ppn=8
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


bcftools filter -i ' FILTER=="PASS"' -o /public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/vcf/pass/<sample>_pass.vcf \
         -O v /public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/vcf/filter_vcf/<sample>.filter.vcf


