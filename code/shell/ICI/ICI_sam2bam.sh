#PBS -N GATK_<sample>
#PBS -l nodes=1:ppn=7
#PBS -l walltime=50:00:00
#PBS -l mem=10gb
#PBS -S /bin/bash
#PBS -q normal_8
#PBS -j oe

if [ -f "/public/home/liuxs/anaconda3/etc/profile.d/conda.sh" ]; then
    . "/public/home/liuxs/anaconda3/etc/profile.d/conda.sh"
else
      export PATH="/public/home/liuxs/anaconda3/bin:$PATH"
fi
conda activate wes

samtools view -bS /path/to/<sample>.sam>\
/path/to/<sample>.bam
