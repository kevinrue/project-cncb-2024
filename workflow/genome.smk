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
rule prepare_reference_genome_fasta:
    output:
        "resources/genome/reference.fa.gz",
    log:
        "logs/prepare_reference_genome_fasta.log",
    shell:
        "({prepare_reference_genome_fasta_cmd_str}) 2> {log}"

# Build GATK index for the combined genome FASTA file.
# Necessary for GATK tools.
rule index_reference_genome_for_gatk:
    input:
        "resources/genome/reference.fa.gz",
    output:
        "resources/genome/reference.dict",
        "resources/genome/reference.fa.gz.fai",
        "resources/genome/reference.fa.gz.gzi",
    log:
        gatkout="logs/index_reference_genome_for_gatk.gatk.out",
        gatkerr="logs/index_reference_genome_for_gatk.gatk.err",
        samtoolkout="logs/index_reference_genome_for_gatk.samtools.out",
        samtoolserr="logs/index_reference_genome_for_gatk.samtools.err",
    shell:
        "gatk CreateSequenceDictionary -R {input} > {log.gatkout} 2> {log.gatkerr} &&"
        " samtools faidx {input} > {log.samtoolkout} 2> {log.samtoolserr}"

# Build BWA index for the combined genome FASTA file.
# Necessary for mapping reads using BWA.
rule index_reference_genome_for_bwa:
    input:
        "resources/genome/reference.fa.gz",
    output:
        idx=multiext("resources/genome/reference.fa.gz", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/index_reference_genome_for_bwa.log",
    wrapper:
        "v5.1.0/bio/bwa/index"

# Map DNA-resequencing reads to the genome using BWA.
rule map_dna_resequencing_reads:
    input:
        reads=expand("reads/genome/{fastq}", fastq=config['genome']['fastqs']),
        idx=multiext("resources/genome/reference.fa.gz", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        "results/genome/mapped.bam",
    log:
        "logs/map_dna_resequencing_reads.log",
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

# Command <https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/master/PairedEndSingleSampleWf-fc-hg38.wdl#L1007>
# ASSUME_SORT_ORDER="queryname" <https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/master/PairedEndSingleSampleWf-fc-hg38.wdl#L1022>
rule mark_duplicates:
    input:
        bams="results/genome/mapped.bam",
    output:
        bam="results/genome/mapped.dedup.bam",
        metrics="results/genome/mapped.dedup.metrics.txt",
    log:
        "logs/mark_duplicates.log",
    params:
        extra=
            " --VALIDATION_STRINGENCY SILENT"
            " --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500"
            " --ASSUME_SORT_ORDER queryname"
            " --CLEAR_DT false"
            " --ADD_PG_TAG_TO_READS false"
    resources:
        mem_mb=1024,
        runtime="1h",
    wrapper:
        "v5.1.0/bio/picard/markduplicates"

# rule sort_bam:
#     input:
#         "results/genome/mapped.dedup.bam",
#     output:
#         "results/genome/mapped.dedup.sorted.bam",
#     log:
#         out="logs/sort_bam.out",
#         err="logs/sort_bam.err",
#     resources:
#         runtime="2h",
#     shell:
#         "picard SortSam"
#         " --INPUT {input}"
#         " --OUTPUT {output}"
#         " --SORT_ORDER coordinate"
#         " --CREATE_INDEX true"
#         " --CREATE_MD5_FILE true"
#         " --MAX_RECORDS_IN_RAM 300000"
#         " > {log.out} 2> {log.err}"

rule sort_dna_resequencing_bam:
    input:
        "results/genome/mapped.dedup.bam",
    output:
        "results/genome/mapped.dedup.sorted.bam",
    log:
        log="logs/sort_dna_resequencing_bam.log",
    resources:
        runtime="30m",
    threads: 8
    wrapper:
        "v5.1.0/bio/samtools/sort"
    
rule haplotype_caller:
    input:
        genome="resources/genome/reference.fa.gz",
        bam="results/genome/mapped.dedup.sorted.bam",
    output:
        gvcf="results/genome/intervals/mapped.dedup.sorted.{interval}.gvcf",
        bamout="results/genome/intervals/mapped.dedup.sorted.{interval}.bamout.bam",
    log:
        out="logs/haplotype_caller/{interval}.out",
        err="logs/haplotype_caller/{interval}.err",
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
        cmd += f" -I results/genome/mapped.dedup.sorted.{interval}.gvcf"
    return cmd

merge_gcvfs_inputs_str = merge_gcvfs_inputs(config)

rule merge_gcvfs:
    input:
        expand("results/genome/intervals/mapped.dedup.sorted.{interval}.gvcf", interval=config['genome']['intervals']),
    output:
        "results/genome/mapped.dedup.sorted.merged.g.vcf.gz",
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

rule genotype_gvcfs:
    input:
        vcf="results/genome/mapped.dedup.sorted.merged.g.vcf.gz",
        genome="resources/genome/reference.fa.gz",
    output:
        "results/genome/mapped.dedup.sorted.merged.vcf.gz",
    log:
        out="logs/genotype_gvcfs.out",
        err="logs/genotype_gvcfs.err",
    resources:
        runtime="10m",
    shell:
        "gatk GenotypeGVCFs"
        " -R {input.genome}"
        " -V {input.vcf}"
        " -O {output}"
        " > {log.out} 2> {log.err}"

rule make_alternate_reference:
    input:
        genome="resources/genome/reference.fa.gz",
        vcf="results/genome/mapped.dedup.sorted.merged.vcf.gz",
    output:
        fasta="results/genome/intervals/mapped.dedup.sorted.{interval}.fa.gz",
    log:
        out="logs/make_alternate_reference.{interval}.out",
        err="logs/make_alternate_reference.{interval}.err",
    resources:
        runtime="4h",
    shell:
        "gatk FastaAlternateReferenceMaker"
        " -R {input.genome}"
        " -O {output.fasta}"
        " -L {wildcards.interval}"
        " -V {input.vcf}"
        " > {log.out} 2> {log.err}"

rule merge_alternate_fastas:
    input:
        expand("results/genome/intervals/mapped.dedup.sorted.{interval}.fa.gz", interval=config['genome']['intervals']),
    output:
        "results/genome/mapped.dedup.sorted.merged.fa.gz",
    log:
        err="logs/merge_alternate_fastas.err",
    resources:
        runtime="10m",
    shell:
        "zcat {input} | bgzip > {output} 2> {log.err}"
