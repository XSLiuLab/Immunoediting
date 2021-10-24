#PBS -N qc_<sample>
#PBS -l nodes=2:ppn=16
#PBS -l walltime=24:00:00
#PBS -S /bin/bash
#PBS -q normal_8
#PBS -j oe


if [ -f "/public/home/liuxs/anaconda3/etc/profile.d/conda.sh" ]; then
    . "/public/home/liuxs/anaconda3/etc/profile.d/conda.sh"
else
      export PATH="/public/home/liuxs/anaconda3/bin:$PATH"
fi
conda activate wes


input=/path/to/fastq
output=/path/to/out

fastp --thread=8 -i $input/<sample>_1.fastq.gz \
      -I $input/<sample>_2.fastq.gz \
      -o $output/<sample>_1_qc.fq.gz \
      -O $output/<sample>_2_qc.fq.gz 
