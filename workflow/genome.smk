# Source 1
# Workflow readme <https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/master/PairedEndSingleSampleWf.md>
# Workflow wdl <https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/master/PairedEndSingleSampleWf.wdl>

# Source 2
# Workflow readme <https://github.com/broadinstitute/warp/tree/develop/pipelines/broad/dna_seq/germline/single_sample/wgs>
# Workflow wdl <https://github.com/broadinstitute/warp/blob/develop/pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl>

# Function called by rule: prepare_reference_genome_fasta
def prepare_reference_genome_fasta_cmd(config):
    # constants
    tmp_dir = "tmp_genome"
    curl_template_cmd = "curl {source_url} > {tmp_dir}/{basename} && "
    # concatenate commands into a string
    cmd = f"mkdir -p {tmp_dir} && "
    for source_url in config["genome"]["fastas"]:
        basename = os.path.basename(source_url)
        cmd += curl_template_cmd.format(source_url=source_url, basename=basename, tmp_dir=tmp_dir)
    cmd += f"zcat {tmp_dir}/*.fa.gz | bgzip > resources/genome/reference.fa.gz && "
    cmd += f"rm -rf {tmp_dir}"
    return cmd

prepare_reference_genome_fasta_cmd_str = prepare_reference_genome_fasta_cmd(config)

# Download a list of genome FASTA files (see config/config.yaml).
# Concatenate them into a single FASTA file.
# Compress it.
# See above for a function that generates the command (cmd_download_genome_fastas).
rule genome_prepare_reference_fasta:
    output:
        "resources/genome/reference.fa.gz",
    log:
        "logs/genome/prepare_reference_fasta.log",
    resources:
        mem="8G",
        runtime="5m",
    threads: 1
    shell:
        "({prepare_reference_genome_fasta_cmd_str}) 2> {log}"

# Build GATK index for the combined genome FASTA file.
# Necessary for GATK tools.
rule genome_index_reference_for_gatk:
    input:
        "resources/genome/reference.fa.gz",
    output:
        "resources/genome/reference.dict",
        "resources/genome/reference.fa.gz.fai",
        "resources/genome/reference.fa.gz.gzi",
    log:
        gatkout="logs/genome/index_reference_for_gatk.gatk.out",
        gatkerr="logs/genome/index_reference_for_gatk.gatk.err",
        samtoolout="logs/genome/index_reference_for_gatk.samtools.out",
        samtoolserr="logs/genome/index_reference_for_gatk.samtools.err",
    resources:
        mem="8G",
        runtime="5m",
    threads: 1
    shell:
        "gatk CreateSequenceDictionary -R {input} > {log.gatkout} 2> {log.gatkerr} &&"
        " samtools faidx {input} > {log.samtoolout} 2> {log.samtoolserr}"

# Build BWA index for the combined genome FASTA file.
# Necessary for mapping reads using BWA.
rule genome_index_reference_for_bwa:
    input:
        "resources/genome/reference.fa.gz",
    output:
        idx=multiext("resources/genome/reference.fa.gz", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/genome/index_reference_for_bwa.log",
    resources:
        mem="8G",
        runtime="5m",
    threads: 1
    wrapper:
        "v5.1.0/bio/bwa/index"

# Map DNA-resequencing reads to the genome using BWA.
rule genome_map_reads:
    input:
        reads=expand("reads/genome/{fastq}", fastq=config['genome']['fastqs']),
        idx=multiext("resources/genome/reference.fa.gz", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        "results/genome/mapped.bam",
    log:
        "logs/genome/map_reads.log",
    params:
        extra=r"-R '@RG\tID:Gdna_1\tSM:Gdna_1'",
        sorting="none",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
        sort_extra="",  # Extra args for samtools/picard.
    resources:
        mem="8G",
        runtime="2h",
    threads: 8
    wrapper:
        "v5.0.1/bio/bwa/mem"

# Command <https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/master/PairedEndSingleSampleWf-fc-hg38.wdl#L1007>
# ASSUME_SORT_ORDER="queryname" <https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/master/PairedEndSingleSampleWf-fc-hg38.wdl#L1022>
rule genome_mark_duplicates:
    input:
        bams="results/genome/mapped.bam",
    output:
        bam="results/genome/mapped.dedup.bam",
        metrics="results/genome/mapped.dedup.metrics.txt",
    log:
        "logs/genome/mark_duplicates.log",
    params:
        extra=
            " --VALIDATION_STRINGENCY SILENT"
            " --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500"
            " --ASSUME_SORT_ORDER queryname"
            " --CLEAR_DT false"
            " --ADD_PG_TAG_TO_READS false"
    resources:
        mem="8G",
        runtime="1h",
    threads: 1
    wrapper:
        "v5.1.0/bio/picard/markduplicates"

rule genome_sort_bam:
    input:
        "results/genome/mapped.dedup.bam",
    output:
        "results/genome/mapped.dedup.sorted.bam",
    log:
        log="logs/genome/sort_bam.log",
    resources:
        mem="8G",
        runtime="30m",
    threads: 8
    wrapper:
        "v5.1.0/bio/samtools/sort"

rule genome_index_bam:
    input:
        "results/genome/mapped.dedup.sorted.bam",
    output:
        "results/genome/mapped.dedup.sorted.bam.bai",
    log:
        "logs/genome/index_bam.log",
    params:
        extra="",
    resources:
        mem="8G",
        runtime="30m",
    threads: 8
    wrapper:
        "v5.1.0/bio/samtools/index"
    
rule genome_haplotype_caller:
    input:
        genome="resources/genome/reference.fa.gz",
        fai="resources/genome/reference.fa.gz.fai",
        bam="results/genome/mapped.dedup.sorted.bam",
        bai="results/genome/mapped.dedup.sorted.bam.bai",
    output:
        gvcf="results/genome/intervals/haplotype_caller.{interval}.gvcf",
        bamout="results/genome/intervals/haplotype_caller.{interval}.bamout.bam",
    log:
        out="logs/genome/haplotype_caller/HaplotypeCaller.{interval}.out",
        err="logs/genome/haplotype_caller/HaplotypeCaller.{interval}.err",
    resources:
        mem="16G",
        runtime="4h",
    threads: 2
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
        cmd += f" -I results/genome/intervals/haplotype_caller.{interval}.gvcf"
    return cmd

merge_gcvfs_inputs_str = merge_gcvfs_inputs(config)

rule genome_merge_gcvfs:
    input:
        expand("results/genome/intervals/haplotype_caller.{interval}.gvcf", interval=config['genome']['intervals']),
    output:
        "results/genome/haplotype_caller.merged.g.vcf.gz",
    log:
        out="logs/genome/merge_gcvfs.out",
        err="logs/genome/merge_gcvfs.err",
    resources:
        mem="8G",
        runtime="15m",
    threads: 2
    shell:
        "gatk SortVcf"
        " {merge_gcvfs_inputs_str}"
        " -O {output}"
        " > {log.out} 2> {log.err}"

rule genome_genotype_gvcfs:
    input:
        vcf="results/genome/haplotype_caller.merged.g.vcf.gz",
        genome="resources/genome/reference.fa.gz",
    output:
        "results/genome/genotype.merged.vcf.gz",
    log:
        out="logs/genome/genotype_gvcfs.out",
        err="logs/genome/genotype_gvcfs.err",
    resources:
        mem="8G",
        runtime="10m",
    threads: 2
    shell:
        "gatk GenotypeGVCFs"
        " -R {input.genome}"
        " -V {input.vcf}"
        " -O {output}"
        " > {log.out} 2> {log.err}"

rule genome_make_alternate_reference:
    input:
        genome="resources/genome/reference.fa.gz",
        vcf="results/genome/genotype.merged.vcf.gz",
    output:
        fasta="results/genome/intervals/alternate_reference.{interval}.fa.gz",
    log:
        out="logs/genome/make_alternate_reference/FastaAlternateReferenceMaker.{interval}.out",
        err="logs/genome/make_alternate_reference/FastaAlternateReferenceMaker.{interval}.err",
    resources:
        mem="8G",
        runtime="10m",
    threads: 2
    shell:
        "gatk FastaAlternateReferenceMaker"
        " -R {input.genome}"
        " -O {output.fasta}"
        " -L {wildcards.interval}"
        " -V {input.vcf}"
        " > {log.out} 2> {log.err}"

rule genome_merge_alternate_reference_fastas:
    input:
        expand("results/genome/intervals/alternate_reference.{interval}.fa.gz", interval=config['genome']['intervals']),
    output:
        "results/genome/alternate_reference.merged.fa.gz",
    log:
        "logs/genome/merge_alternate_reference_fastas.log",
    resources:
        mem="2G",
        runtime="10m",
    threads: 1
    shell:
        "zcat {input} | bgzip > {output} 2> {log}"
