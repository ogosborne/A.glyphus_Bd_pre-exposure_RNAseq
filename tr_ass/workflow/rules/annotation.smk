# transdecoder
rule transdecoder_longorfs:
    input:
        fasta = output_dir + "trinity/" + asm_name + ".Trinity.fasta",
        gene_trans_map = output_dir + "trinity/" + asm_name + ".Trinity.fasta.gene_trans_map"
    output:
        temp(output_dir + "transdecoder/" + asm_name + ".Trinity.fasta.transdecoder_dir/longest_orfs.pep"),
    log:
        output_dir + "logs/transdecoder/longorfs.log"
    threads: 1
    resources:
        runtime = "4h",
        mem_mb = 12000,
    params:
        outdir = output_dir + "transdecoder/",
    conda:
        "../envs/transdecoder.yaml"
    shell: """
        TransDecoder.LongOrfs -t {input.fasta} --gene_trans_map {input.gene_trans_map} --output_dir {params.outdir}
    """

rule transdecoder_predict:
    input:
        fasta = output_dir + "trinity/" + asm_name + ".Trinity.fasta",
        longorfs = output_dir + "transdecoder/" + asm_name + ".Trinity.fasta.transdecoder_dir/longest_orfs.pep",
    output:
        output_dir + "transdecoder/" + asm_name + ".Trinity.fasta.transdecoder.bed",
        output_dir + "transdecoder/" + asm_name + ".Trinity.fasta.transdecoder.cds",
        output_dir + "transdecoder/" + asm_name + ".Trinity.fasta.transdecoder.pep",
        output_dir + "transdecoder/" + asm_name + ".Trinity.fasta.transdecoder.gff3",
    log:
        output_dir + "logs/transdecoder/predict.log"
    threads: 1
    resources:
        runtime = "4h",
        mem_mb = 12000,    
    params:
        outdir = output_dir + "transdecoder/",
        extra = "--single_best_only ",
        tmpdir =  output_dir + "transdecoder/" + asm_name + ".Trinity.fasta.transdecoder_dir/",
    conda:
        "../envs/transdecoder.yaml"
    shell: """
        TransDecoder.Predict -t {input.fasta} --output_dir {params.outdir} {params.extra}
        rm -r {params.tmpdir}
    """

## eggnog-mapper
rule emapper:
    input:
        # fasta
        fasta = output_dir + "transdecoder/" + asm_name + ".Trinity.fasta.transdecoder.pep",
    output:
        annots = output_dir + "emapper/" + asm_name + ".emapper.annotations",
        hits = output_dir + "emapper/" + asm_name + ".emapper.hits",
        outdir = directory(output_dir + "emapper/"),
    params: 
        db = config["emapper"]["db_dir"],
        asm_name = asm_name,
        tax_scope = str(config["emapper"]["tax_scope"]),
    log:
        output_dir + "logs/emapper/emapper.log"
    threads: 20
    resources:
        runtime = "4h",
        mem_mb = 160000,  
    conda: "../envs/emapper.yaml"
    shell: """
    emapper.py \
    -i {input.fasta} \
    --itype proteins \
    --tax_scope {params.tax_scope} \
    --output_dir {output.outdir} \
    --data_dir {params.db} \
    --output {params.asm_name} \
    -m diamond \
    --cpu {threads}
   
    """