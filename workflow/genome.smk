def cmd_download_genome_fastas(config):
    tmp_dir = "tmp_genome"
    curl_template_cmd = "curl {source_url} > {tmp_dir}/{basename} && "
    cmd = f"mkdir -p {tmp_dir} && "
    for source_url in config["genome"]["fastas"]:
        basename = os.path.basename(source_url)
        cmd += curl_template_cmd.format(source_url=source_url, basename=basename, tmp_dir=tmp_dir)
    cmd += f"zcat {tmp_dir}/*.fa.gz > resources/genome/genome.fa && "
    cmd += "bgzip resources/genome/genome.fa && "
    cmd += f"rm -rf {tmp_dir}"
    return cmd

cmd_download_genome_fastas = cmd_download_genome_fastas(config)

rule prepare_genome_fasta:
    output:
        "resources/genome/genome.fa.gz",
    log:
        "logs/prepare_genome_fasta.log",
    shell:
        "{cmd_download_genome_fastas} 2> {log}"

rule prepare_genome_gatk_index:
    input:
        "resources/genome/genome.fa.gz",
    output:
        "resources/genome/genome.dict",
        "resources/genome/genome.fa.gz.fai",
        "resources/genome/genome.fa.gz.gzi",
    log:
        gatkout="logs/prepare_genome_gatk_index.gatk.out",
        gatkerr="logs/prepare_genome_gatk_index.gatk.err",
        samtoolkout="logs/prepare_genome_gatk_index.samtools.out",
        samtoolserr="logs/prepare_genome_gatk_index.samtools.err",
    shell:
        "gatk CreateSequenceDictionary -R {input} > {log.gatkout} 2> {log.gatkerr} && "
        "samtools faidx {input} > {log.samtoolkout} 2> {log.samtoolserr}"
