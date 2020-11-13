#PBS -N mark_<sample>
#PBS -l nodes=1:ppn=8
#PBS -l walltime=100:00:00
#PBS -S /bin/bash
#PBS -q normal_3
#PBS -j oe


if [ -f "/public/home/liuxs/anaconda3/etc/profile.d/conda.sh" ]; then
    . "/public/home/liuxs/anaconda3/etc/profile.d/conda.sh"
else
      export PATH="/public/home/liuxs/anaconda3/bin:$PATH"
fi

conda activate wes



java -jar /public/home/liuxs/anaconda3/envs/wes/share/gatk4-4.1.3.0-0/gatk-package-4.1.3.0-local.jar MarkDuplicates \
        -I=/public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/worked_bam/sort_bam/<sample>.sort.bam \
        -O=/public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/worked_bam/mark_bam/<sample>.rmdup.bam \
        --VALIDATION_STRINGENCY=LENIENT \
        --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
        -M=/public/home/wangshx/wt/immune_therapy_dataset/Nadeem_Riaz_2017_cell/files/worked_bam/mark_bam/<sample>.sort.addhead.rmdup.metric

