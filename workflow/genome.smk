# Workflow readme <https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/master/PairedEndSingleSampleWf.md>
# Workflow wdl <https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/master/PairedEndSingleSampleWf.wdl>

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
        "({cmd_download_genome_fastas_str}) 2> {log}"

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
        "gatk CreateSequenceDictionary -R {input} > {log.gatkout} 2> {log.gatkerr} &&"
        " samtools faidx {input} > {log.samtoolkout} 2> {log.samtoolserr}"

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

## The rule below uses a wrapper that seems to require 'conda' installed in the container
# rule map_reads_to_genome:
#     input:
#         reads=[
#             "/ceph/project/cncb/shared/proj140/analyses/novogene_sequencing/genome/download/X204SC24080649-Z01-F001/01.RawData/Gdna_1/Gdna_1_EKDN240047675-1A_22FVLJLT4_L5_1.fq.gz",
#             "/ceph/project/cncb/shared/proj140/analyses/novogene_sequencing/genome/download/X204SC24080649-Z01-F001/01.RawData/Gdna_1/Gdna_1_EKDN240047675-1A_22FVLJLT4_L5_2.fq.gz"],
#         idx=multiext("resources/genome/genome.fa.gz", ".amb", ".ann", ".bwt", ".pac", ".sa"),
#     output:
#         "resources/genome_resequencing/mapped.bam",
#     log:
#         "logs/map_reads_to_genome.log",
#     params:
#         extra=r"-R '@RG\tID:Gdna_1\tSM:Gdna_1'",
#         sorting="none",  # Can be 'none', 'samtools' or 'picard'.
#         sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
#         sort_extra="",  # Extra args for samtools/picard.
#     threads: 8
#     resources:
#         runtime="2h",
#     wrapper:
#         "v5.0.1/bio/bwa/mem"

rule map_reads_to_genome:
    input:
        reads=[
            "/ceph/project/cncb/shared/proj140/analyses/novogene_sequencing/genome/download/X204SC24080649-Z01-F001/01.RawData/Gdna_1/Gdna_1_EKDN240047675-1A_22FVLJLT4_L5_1.fq.gz",
            "/ceph/project/cncb/shared/proj140/analyses/novogene_sequencing/genome/download/X204SC24080649-Z01-F001/01.RawData/Gdna_1/Gdna_1_EKDN240047675-1A_22FVLJLT4_L5_2.fq.gz"],
        genome="resources/genome/genome.fa.gz",
        amb="resources/genome/genome.fa.gz.amb",
        ann="resources/genome/genome.fa.gz.ann",
        bwt="resources/genome/genome.fa.gz.bwt",
        pac="resources/genome/genome.fa.gz.pac",
        sa="resources/genome/genome.fa.gz.sa",
    output:
        "resources/genome_resequencing/mapped.bam",
    log:
        "logs/map_reads_to_genome.log",
    threads: 8
    resources:
        runtime="2h",
    shell:
        "(bwa mem"
        " -t {threads}"
        " -R '@RG\\tID:Gdna_1\\tSM:Gdna_1'"
        " {input.genome}"
        " {input.reads}"
        " | samtools view -@ {threads} -b - > {output}) 2> {log}"

# Command <https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/master/PairedEndSingleSampleWf-fc-hg38.wdl#L1007>
# ASSUME_SORT_ORDER="queryname" <https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/master/PairedEndSingleSampleWf-fc-hg38.wdl#L1022>
rule mark_duplicates:
    input:
        "resources/genome_resequencing/mapped.bam",
    output:
        bam="resources/genome_resequencing/mapped.dedup.bam",
        metrics="resources/genome_resequencing/mapped.dedup.metrics.txt",
    log:
        out="logs/mark_duplicates.out",
        err="logs/mark_duplicates.err",
    resources:
        runtime="1h",
    shell:
        "picard MarkDuplicates"
        " --INPUT {input}"
        " --OUTPUT {output.bam}"
        " --METRICS_FILE {output.metrics}"
        " --VALIDATION_STRINGENCY SILENT"
        " --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500"
        " --ASSUME_SORT_ORDER \"queryname\""
        " --CLEAR_DT \"false\""
        " --ADD_PG_TAG_TO_READS false"
        " > {log.out} 2> {log.err}"

rule sort_bam:
    input:
        "resources/genome_resequencing/mapped.dedup.bam",
    output:
        "resources/genome_resequencing/mapped.dedup.sorted.bam",
    log:
        out="logs/sort_bam.out",
        err="logs/sort_bam.err",
    resources:
        runtime="2h",
    shell:
        "picard SortSam"
        " --INPUT {input}"
        " --OUTPUT {output}"
        " --SORT_ORDER coordinate"
        " --CREATE_INDEX true"
        " --CREATE_MD5_FILE true"
        " --MAX_RECORDS_IN_RAM 300000"
        " > {log.out} 2> {log.err}"
    
rule haplotype_caller:
    input:
        genome="resources/genome/genome.fa.gz",
        bam="resources/genome_resequencing/mapped.dedup.sorted.bam",
    output:
        gvcf="resources/genome_sequencing/mapped.dedup.sorted.{interval}.gvcf",
        bamout="resources/genome_sequencing/mapped.dedup.sorted.{interval}.bamout.bam",
    log:
        out="logs/haplotype_caller.{interval}.out",
        err="logs/haplotype_caller.{interval}.err",
    resources:
        runtime="4h",
    shell:
        "gatk HaplotypeCaller"
        " -R {input.genome}"
        " -I {input.bam}"
        " -L {wildcards.interval}"
        " -O {output.gvcf}"
        " -contamination 0"
        " -G StandardAnnotation"
        " -G StandardHCAnnotation"
        " -G AS_StandardAnnotation"
        " -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90"
        " -ERC GVCF"
        " -bamout {output.bamout}"
        " > {log.out} 2> {log.err}"

def merge_gcvfs_inputs(config):
    cmd = ""
    for interval in config['genome']['intervals']:
        cmd += f" -I resources/genome_resequencing/mapped.dedup.sorted.{interval}.gvcf"
    return cmd

merge_gcvfs_inputs_str = merge_gcvfs_inputs(config)

rule merge_gcvfs:
    input:
        expand("resources/genome_resequencing/mapped.dedup.sorted.{interval}.gvcf", interval=config['genome']['intervals']),
    output:
        "resources/genome_resequencing/mapped.dedup.sorted.merged.g.vcf.gz",
    log:
        out="logs/merge_gcvfs.out",
        err="logs/merge_gcvfs.err",
    resources:
        runtime="10m",
    shell:
        "gatk SortVcf"
        " {merge_gcvfs_inputs_str}"
        " -O {output}"
        " > {log.out} 2> {log.err}"
