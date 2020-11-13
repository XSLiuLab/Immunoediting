#PBS -N annotate_pancacner_mhc_i
#PBS -l nodes=1:ppn=10
#PBS -l walltime=24:00:00
#PBS -S /bin/bash
#PBS -q normal_8
#PBS -j oe

source activate neo

vcf_path=/public/home/wangshx/wx/neoantigens_delption/find_neoantigens_delotion_signal/vcf_file

out_path=/public/home/wangshx/wx/neoantigens_delption/find_neoantigens_delotion_signal/vcf_annoated/pancancer_maf/mhc_i

vep  --input_file $vcf_path/<sample>.vcf \
 --output_file $out_path/<sample>.annoated.vcf  \
--format vcf --vcf --symbol --term SO \
--hgvs \
--fasta  /public/home/liuxs/.vep/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--offline --cache --dir_cache /public/home/liuxs/.vep/ --cache_version 98 --assembly GRCh38 \
--plugin Downstream --plugin Wildtype \
--dir_plugin /public/home/wangshx/wt/VEP_plugins \
--transcript_version
