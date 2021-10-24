#PBS -N mem_<sample>
#PBS -l nodes=1:ppn=9
#PBS -l walltime=20:00:00
#PBS -S /bin/bash
#PBS -q normal_8
#PBS -j oe

if [ -f "/public/home/liuxs/anaconda3/etc/profile.d/conda.sh" ]; then
    . "/public/home/liuxs/anaconda3/etc/profile.d/conda.sh"
else
      export PATH="/public/home/liuxs/anaconda3/bin:$PATH"
fi

conda activate wes

workdir=/path/to/qc_file

cd $workdir
bwa mem -t 8 -M -R "@RG\tID:<sample>\t\
LM:<sample>\t\
SM:<sample>\t\
PL:illumina\tPU:<sample>"\
 /path/to/hg38.fa <sample>_1_qc.fq.gz\
 <sample>_2_qc.fq.gz\
> /path/to/<sample>.sam 
