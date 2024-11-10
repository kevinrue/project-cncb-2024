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

cmd_download_genome_fastas_str = cmd_download_genome_fastas(config)

rule prepare_genome_fasta:
    output:
        "resources/genome/genome.fa.gz",
    log:
        "logs/prepare_genome_fasta.err",
    shell:
        "{cmd_download_genome_fastas_str} 2> {log}"

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

rule prepare_genome_bwa_index:
    input:
        "resources/genome/genome.fa.gz",
    output:
        "resources/genome/genome.fa.gz.amb",
        "resources/genome/genome.fa.gz.ann",
        "resources/genome/genome.fa.gz.bwt",
        "resources/genome/genome.fa.gz.pac",
        "resources/genome/genome.fa.gz.sa",
    log:
        out="logs/prepare_genome_bwa_index.out",
        err="logs/prepare_genome_bwa_index.err",
    shell:
        "bwa index {input} > {log.out} 2> {log.err}"

rule map_reads_to_genome:
    input:
        reads=[
            "/ceph/project/cncb/shared/proj140/analyses/novogene_sequencing/genome/download/X204SC24080649-Z01-F001/01.RawData/Gdna_1/Gdna_1_EKDN240047675-1A_22FVLJLT4_L5_1.fq.gz",
            "/ceph/project/cncb/shared/proj140/analyses/novogene_sequencing/genome/download/X204SC24080649-Z01-F001/01.RawData/Gdna_1/Gdna_1_EKDN240047675-1A_22FVLJLT4_L5_2.fq.gz"],
        idx=multiext("resources/genome/genome.fa.gz", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        "resources/genome_resequencing/mapped.bam",
    log:
        "logs/map_reads_to_genome.log",
    params:
        extra=r"-R '@RG\tID:Gdna_1\tSM:Gdna_1'",
        sorting="none",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
        sort_extra="",  # Extra args for samtools/picard.
    threads: 8
    resources:
        runtime="2h",
    wrapper:
        "v5.0.1/bio/bwa/mem"
