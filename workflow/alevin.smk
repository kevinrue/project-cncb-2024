# Source <https://combine-lab.github.io/alevin-fry-tutorials/2023/simpleaf-piscem/>
# define env variable
#export ALEVIN_FRY_HOME="$AF_SAMPLE_DIR/af_home"
# simpleaf configuration
#simpleaf set-paths
#ulimit -n 2048

# Function called by rule: prepare_reference_genome_fasta
def download_gtf_cmd(config):
    source_url = config["genome"]["gtf"]
    basename = os.path.basename(source_url)
    cmd = f"curl {source_url} > resources/genome/reference.gtf.gz"
    return cmd

download_gtf_cmd_str = download_gtf_cmd(config)

rule download_gtf:
    output:
        "resources/genome/reference.gtf.gz",
    log:
        "logs/genome/download_gtf.log",
    resources:
        mem="2G",
        runtime="5m",
    shell:
        "({download_gtf_cmd_str}) 2> {log}"

rule alevin_build_reference_index:
    input:
        genome="resources/genome/reference.fa.gz",
        gtf="resources/genome/reference.gtf.gz",
    output:
        index=directory("resources/genome/index/alevin"),
    log:
        out="logs/genome/alevin/build_reference_index.out",
        err="logs/genome/alevin/build_reference_index.err",
    threads: 16
    resources:
        mem="8G",
        runtime="4h",
    shell:
        "jobdir=$(pwd) &&"
        " cd $TMPDIR &&"
        " export ALEVIN_FRY_HOME=af_home &&"
        " simpleaf set-paths &&"
        " gunzip -c $jobdir/{input.genome} > tmp_reference.fa  &&"
        " simpleaf index"
        " --output tmp_alevin_index"
        " --fasta tmp_reference.fa"
        " --gtf $jobdir/{input.gtf}"
        " --rlen 150"
        " --threads {threads}"
        " --use-piscem"
        " > $jobdir/{log.out} 2> $jobdir/{log.err} &&"
        " mv tmp_alevin_index $jobdir/{output.index} &&"
        " rm tmp_reference.fa"

# Function used by rule: alevin_quant
def get_alevin_quant_input_fastqs(wildcards):
    sample_fastqs_info=SAMPLES[SAMPLES['sample_name'] == 'WPPm048hrs_rep1'].filter(items=['directory', 'R1', 'R2'])
    fq1=[os.path.join(os.path.realpath(sample_fastqs_info['directory'][i]), sample_fastqs_info['R1'][i]) for i in range(sample_fastqs_info.shape[0])]
    fq2=[os.path.join(os.path.realpath(sample_fastqs_info['directory'][i]), sample_fastqs_info['R2'][i]) for i in range(sample_fastqs_info.shape[0])]
    return {
        'fq1': fq1,
        'fq2': fq2,
    }


rule alevin_quant:
    input:
        unpack(get_alevin_quant_input_fastqs),
        index="resources/genome/index/alevin",
    output:
        "results/alevin/{sample}",
    params:
        reads1=lambda wildcards, input: ','.join(input.fq1),
        reads2=lambda wildcards, input: ','.join(input.fq2),
    log:
        out="logs/genome/alevin/quant/{sample}.out",
        err="logs/genome/alevin/quant/{sample}.err",
    threads: 16
    resources:
        mem="8G",
        runtime="4h",
    shell:
        "jobdir=$(pwd) &&"
        " cd $TMPDIR &&"
        " export ALEVIN_FRY_HOME=af_home &&"
        " simpleaf set-paths &&"
        " simpleaf quant"
        " --reads1 {params.reads1}"
        " --reads2 {params.reads2}"
        " --index $jobdir/{input.index}/index"
        " --chemistry 10xv3 --resolution cr-like --expected-ori fw --unfiltered-pl" # from tutorial, to be confirmed
        " --t2g-map $jobdir/{input.index}/index/t2g_3col.tsv"
        " --threads {threads}"
        " --output alevin_quant_{wildcards.sample}"
        " > $jobdir/{log.out} 2> $jobdir/{log.err} &&"
        " mv alevin_quant_{wildcards.sample} $jobdir/{output}"
