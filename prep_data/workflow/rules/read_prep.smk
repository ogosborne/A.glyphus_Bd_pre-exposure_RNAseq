######################################
# trim, filter and QC with fastp
rule fastp_pe:
    input:
        sample = get_fastq,
    output:
        trimmed = [output_dir + "trimmed/{sample}.R1.fastq.gz", output_dir + "trimmed/{sample}.R2.fastq.gz",],
        # Unpaired reads separately
        unpaired1 = temp(output_dir + "trimmed/{sample}_U1.fastq.gz"),
        unpaired2 = temp(output_dir + "trimmed/{sample}_U2.fastq.gz"),
        # Reports
        html = output_dir + "trimmed/report/{sample}.fastp.html",
        json = output_dir + "trimmed/report/{sample}.fastp.json"
    log:
        output_dir + "logs/fastp/{sample}.log"
    params:
        extra = "--detect_adapter_for_pe"
    threads: 6
    resources:
        runtime = "1h"
    wrapper:
        "v4.7.1/bio/fastp"

# summarise fastp outputs with multiqc 
rule multiqc_read_prep:
    input:
        expand(output_dir + "trimmed/report/{sample}.fastp.json", sample = samples),
    output:
        output_dir + "reports/read_trimming.html",
        directory(output_dir + "reports/multiqc_data"),
    params:
        extra = ""
    resources:
        runtime = "20m",
    log:
        output_dir+"logs/multiqc/read_trimming.log"
    wrapper:
        "v4.7.1/bio/multiqc"

## split reads with seal.sh
rule seal:
    input:
        # reads
        R1 = output_dir+"trimmed/{sample}.R1.fastq.gz",
        R2 = output_dir+"trimmed/{sample}.R2.fastq.gz",
        # fasta with rRNA sequences
        rrna = config["seal"]["rrna"],
        # fasta with known contaminant sequences, e.g. a pathogen genome
        known_contam = config["seal"]["known_contam"],
    output:
        # interleaved reads
        rrna_reads = temp(output_dir + "binned/rrna.{sample}.fq.gz"),
        contam_reads = temp(output_dir + "binned/contam.{sample}.fq.gz"),
        clean_reads = temp(output_dir + "binned/clean.{sample}.fq.gz"),
        # deinterleaved reads
        rrna_R1 = output_dir + "binned/rrna.{sample}.R1.fq.gz",
        rrna_R2 = output_dir + "binned/rrna.{sample}.R2.fq.gz",
        contam_R1 = output_dir + "binned/contam.{sample}.R1.fq.gz",
        contam_R2 = output_dir + "binned/contam.{sample}.R2.fq.gz",
        clean_R1 = output_dir + "binned/clean.{sample}.R1.fq.gz",
        clean_R2 = output_dir + "binned/clean.{sample}.R2.fq.gz",
        # binning stats
        refstats = output_dir + "binned/{sample}.refstats",
        # soft links to references so seal names outputs correctly
        rrna_ln = temp(output_dir + "binned/{sample}_temp/rrna.fasta"),
        contam_ln = temp(output_dir + "binned/{sample}_temp/contam.fasta"),
        tmpdir = temp(directory(output_dir + "binned/{sample}_temp/"))
    params: 
        speed = config["seal"]["speed"],
        pattern = output_dir + "binned/%.{sample}.fq.gz",
    log:
        output_dir + "logs/seal/{sample}.log"
    threads: 10
    resources:
        runtime = "2h",
        mem_mb = get_mem_mb_seal
    conda: "../envs/bbtools.yaml"
    shell: """
    
    # make renamed links for refs
    ln -s {input.rrna} {output.rrna_ln}
    ln -s {input.known_contam} {output.contam_ln}

    # run seal
    seal.sh in={input.R1} in2={input.R2} \
    ref={output.rrna_ln},{output.contam_ln} \
    ambiguous=first \
    refnames=t \
    pattern={params.pattern} \
    outu={output.clean_reads} \
    refstats={output.refstats} \
    threads={threads} \
    Xmx={resources.mem_mb}m \
    prealloc=t \
    speed={params.speed} &> {log}
    
    # deinterleave output reads
    reformat.sh in={output.rrna_reads} \
    out1={output.rrna_R1} \
    out2={output.rrna_R2}
    reformat.sh in={output.contam_reads} \
    out1={output.contam_R1} \
    out2={output.contam_R2}
    reformat.sh in={output.clean_reads} \
    out1={output.clean_R1} \
    out2={output.clean_R2}

    """
# summarise reference statistics and plot
rule refstats:
    input:
        expand(output_dir + "binned/{sample}.refstats", sample = samples),
    output: 
        tab = output_dir + "reports/bin_reads_by_ref.csv",
        fig = output_dir + "reports/bin_reads_by_ref.pdf"
    log:
        output_dir + "logs/seal/refstats.log"
    threads: 1
    resources:
        runtime = "10m",
    conda: "../envs/r.yaml"
    script: "../scripts/refstats.R"


