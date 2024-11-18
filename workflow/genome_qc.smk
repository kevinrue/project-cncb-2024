rule genome_reads_qc:
    input:
        "reads/genome/{fastqnoext}.fq.gz",
    output:
        html="qc/genome/fastqc/{fastqnoext}.html",
        zip="qc/genome/fastqc/{fastqnoext}_fastqc.zip"
    log:
        "logs/genome/fastqc/{fastqnoext}.log",
    threads: 8
    resources:
        runtime="1h",
        mem_mb = 1024,
    wrapper:
        "v5.1.0/bio/fastqc"

rule genome_multiqc:
    input:
        expand("qc/genome/fastqc/{fastqnoext}_fastqc.zip", fastqnoext=genome_fastqs_noext),
        "results/genome/alternate_reference.merged.fa.gz",
    output:
        "qc/genome/multiqc/multiqc_report.html",
    log:
        out="logs/genome/multiqc.out",
        err="logs/genome/multiqc.err",
    resources:
        runtime="10m",
        mem_mb = 1024,
    shell:
        "multiqc"
        " resources/genome"
        " results/genome"
        " qc/genome/fastqc"
        " -o qc/genome/multiqc"
        " --force"
        " > {log.out} 2> {log.err}"
