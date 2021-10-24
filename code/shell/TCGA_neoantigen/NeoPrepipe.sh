##sample names
awk '{print $1}' panCancer_hla.tsv > file_names 

mkdir split
cd split/
##split
split -89 ../file_names
ls x* > split_files

cd ../batch
cat ../split/split_files | while read i;do mkdir $i ;done

###move
#!/bin/bash
cat ../split/split_files | while read i
do
  cat ../split/$i | while read j
  do
    mv ../nonsynonmous/$j.vcf ./$i
  done
done

cat ../split/split_files | while read i;do mkdir $i ;done
cat ../split/split_files  | while read i; do cp ../TCGA_HLA_typing.txt ./$i;done

##the pbs file do neoantigen prediction
#PBS -N run_NeoPredPipe
#PBS -l nodes=1:ppn=6
#PBS -l walltime=480:00:00
#PBS -S /bin/bash
#PBS -q pub_fast
#PBS -j oe

if [ -f "/public/slst/home/wutao2/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/public/slst/home/wutao2/miniconda3/etc/profile.d/conda.sh"
else
        export PATH="/public/slst/home/wutao2/miniconda3/bin:$PATH"
fi

conda activate python27
NeoPredPipe.py -I /public/slst/home/wutao2/TCGA_neopredpipe/batch/<dir> -H /public/slst/home/wutao2/TCGA_neopredpipe/batch_hla/<dir>/TCGA_HLA_typing.txt -o /public/slst/home/wutao2/TCGA_neopredpipe/batch_results/<dir> -n batch_<dir> -c 0 -E 8 9 10 -x /public/slst/home/wutao2/TCGA_neopredpipe/exp -a -m

###generate pbs files in batch
sync-qgen -m mapping -s split/split_files -f run_neopredpipe.pbs -d pbs
nohup ~/gosub ../pbs &

