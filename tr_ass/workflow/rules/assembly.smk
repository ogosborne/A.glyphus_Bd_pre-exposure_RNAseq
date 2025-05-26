rule trinity:
    input:
        left = R1,
        right = R2,
    output:
        dir = temp(directory(output_dir + "trinity/" + asm_name)),
        fas = output_dir + "trinity/" + asm_name + ".Trinity.fasta",
        map = output_dir + "trinity/" + asm_name + ".Trinity.fasta.gene_trans_map",
    log:
        output_dir + 'logs/trinity/trinity.log',
    params:
        extra = "",
    threads: 16
    resources:
        mem_mb = get_mem_mb_trinity,
        runtime = "3d",
    wrapper:
        "v4.7.1/bio/trinity"
