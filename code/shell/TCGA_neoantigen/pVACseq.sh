#PBS -N predict_neo
#PBS -l nodes=1:ppn=20
#PBS -l walltime=24:00:00
#PBS -S /bin/bash
#PBS -q pub_fast
#PBS -j oe
 
if [ -f "/public/slst/home/wutao2/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/public/slst/home/wutao2/miniconda3/etc/profile.d/conda.sh"
else
        export PATH="/public/slst/home/wutao2/miniconda3/bin:$PATH"
fi 

conda activate pvactools
 
pvacseq run \
/public/slst/home/wutao2/anno_vcf/<samples_file> \
<samples_name> \
<HLA> \
MHCflurry MHCnuggetsI NetMHCpan \
/public/slst/home/wutao2/TCGA_pvacseq/out \
-e1 8,9,10 \
--iedb-install-directory /public/slst/home/wutao2/ \
-a sample_name --pass-only -t 20