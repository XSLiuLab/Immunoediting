#PBS -N sra2fastq_<sample>
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

workdir=/path/to/sra
cd $workdir

fastq-dump --split-3 --gzip  <sample>.sra -O fastq
