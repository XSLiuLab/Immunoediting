vcf=/public/home/wangshx/wx/neoantigens_delption/find_neoantigens_delotion_signal/vcf_annoated/pancancer/mhc-i/vcf

out_dir=/public/home/wangshx/wx/neoantigens_delption/find_neoantigens_delotion_signal/neo_pre/one_method/TCGA-QH-A870-01A

iedb_dir=/public/home/liuxs/biosoft/iedb

pvacseq run \
$vcf/TCGA-QH-A870-01A.annoated.vcf \
TCGA-QH-A870-01A \
HLA-A*02:01,HLA-A*02:01,HLA-B*13:02,HLA-B*51:01,HLA-C*04:01,HLA-C*14:02 \
NetMHCpan \
$out_dir \
-e 9 \
--iedb-install-directory $iedb_dir \
--fasta-size 50