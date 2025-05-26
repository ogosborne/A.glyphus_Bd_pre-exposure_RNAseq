rule busco:
    input:
        output_dir + "transdecoder/" + asm_name + ".Trinity.fasta.transdecoder.pep",
    output:
        short_json = output_dir + "busco/" + asm_name + "/short_summary.specific." + config["busco"]["lineage"] + "." + asm_name + ".json",
        short_txt = output_dir + "busco/" + asm_name + "/short_summary.specific." + config["busco"]["lineage"] + "." + asm_name + ".txt",
    log:
        output_dir + "logs/busco/busco.log",
    params:
        lineage = config["busco"]["lineage"],
        dataset_dir = config["busco"]["db_dir"],
        out_path = output_dir + "busco/",
        out_name = asm_name,
    threads: 8
    resources:
        runtime = "4h",
        mem_mb = 12000, 
    conda:
        "../envs/busco.yaml"
    shell:"""
    busco --cpu {threads} \
    --in {input} \
    --mode proteins \
    --offline \
    --download_path {params.dataset_dir} \
    --lineage {params.lineage} \
    --out {params.out_name} \
    --out_path {params.out_path} > {log}
    """

rule transrate:
    input:
        fasta = output_dir + "trinity/" + asm_name + ".Trinity.fasta",
    output:
        output_dir + "transrate/assemblies.csv",
        output_dir + "transrate/" + asm_name + ".Trinity/contigs.csv"
    log:
        output_dir + "logs/transrate/transrate.log",
    threads: 8
    resources:
        runtime = "12h",
        mem_mb = 160000,
    conda:
        "../envs/transrate.yaml"
    shell: """
    transrate --assembly={input.fasta} --threads={threads} --output={output} > {log}
    """

     