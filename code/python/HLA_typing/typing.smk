rule typing:
    input: 
        fa = "/path/to/fish_fastq",
        config = "/public/home/wangshx/wx/neoantigens_delption/find_neoantigens_delotion_signal/demo/HLA_typing/snakemake/config.ini"
    params:
        op = "-d -o /output/dir"
    shell:
        "OptiTypePipeline.py -i {input.fa} {params.op} -c {input.config}" 