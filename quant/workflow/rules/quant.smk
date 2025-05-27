rule salmon_index:
    input:
        sequences = config["transcriptome"],
    output:
        multiext(
            output_dir + "salmon/transcriptome_index/",
            "complete_ref_lens.bin",
            "ctable.bin",
            "ctg_offsets.bin",
            "duplicate_clusters.tsv",
            "info.json",
            "mphf.bin",
            "pos.bin",
            "pre_indexing.log",
            "rank.bin",
            "refAccumLengths.bin",
            "ref_indexing.log",
            "reflengths.bin",
            "refseq.bin",
            "seq.bin",
            "versionInfo.json",
        ),
        directory(output_dir + "salmon/transcriptome_index/"),
    log:
        output_dir + "logs/salmon/transcriptome_index.log",
    threads: 2
    params:
        # optional parameters
        extra="",
    resources:
        runtime = "2h",
        mem_mb = 24000, 
    wrapper:
        "v4.7.5/bio/salmon/index"

rule salmon_quant_reads:
    input:
        # If you have multiple fastq files for a single sample (e.g. technical replicates)
        # use a list for r1 and r2.
        r1 = get_fastq1,
        r2 = get_fastq2,
        index = output_dir + "salmon/transcriptome_index/",
    output:
        quant = output_dir + "salmon/{sample}/quant.sf",
        lib = output_dir + "salmon/{sample}/lib_format_counts.json",
    log:
        output_dir + "logs/salmon/{sample}.log",
    params:
        # optional parameters
        libtype = "A",
        extra = "",
    threads: 2
    resources:
        runtime = "2h",
        mem_mb = 24000, 
    wrapper:
        "v4.7.5/bio/salmon/quant"

rule tximport:
    input:
        quant = expand(output_dir + "salmon/{sample}/quant.sf", sample = samples),
        tx_to_gene = config["tx2gene"],
    output:
        txi = output_dir + "salmon/txi.RDS",
    params:
        extra = "type='salmon', txOut=TRUE",
    threads: 1
    resources:
        runtime = "1h",
        mem_mb = 12000, 
    wrapper:
        "v4.7.5/bio/tximport"

# summarise salmon mappings with multiqc 
#rule multiqc_salmon:
#    input:
#        expand(output_dir + "trimmed/report/{sample}.fastp.json", sample = samples),
#    output:
#        output_dir + "reports/read_trimming.html",
#        directory(output_dir + "reports/multiqc_data"),
#    params:
#        extra = ""
#    resources:
#        runtime = "20m",
#    log:
#        output_dir+"logs/multiqc/read_trimming.log"
#    wrapper:
#        "v4.7.1/bio/multiqc"