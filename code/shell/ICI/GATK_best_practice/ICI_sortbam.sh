#PBS -N GATK_<sample>
#PBS -l nodes=2:ppn=4
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

workdir=/path/to/bam


java -Djava.io.tmpdir=/public/home/wangshx/wx/tmp -jar /public/home/liuxs/anaconda3/envs/wes/share/gatk4-4.1.3.0-0/gatk-package-4.1.3.0-local.jar SortSam \
        -I=$workdir/<sample>.bam \
        -O=$workdir/<sample>.sort.bam \
        --SORT_ORDER=coordinate
