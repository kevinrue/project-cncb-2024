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
        runtime="6h",
    shell:
        "export ALEVIN_FRY_HOME=af_home &&"
        " simpleaf set-paths &&"
        " gunzip -c {input.genome} > tmp_reference.fa  &&"
        " simpleaf index"
        " --output {output.index}"
        " --fasta tmp_reference.fa"
        " --gtf {input.gtf}"
        " --rlen 150"
        " --threads 16"
        " --use-piscem"
        " > {log.out} 2> {log.err} &&"
        " rm tmp_reference.fa"
