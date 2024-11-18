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
        "logs/genome/alevin/build_reference_index.log",
    threads: 16
    resources:
        runtime="1h",
        mem_mb=8192,
    shell:
        "simpleaf index"
        " --output {output.index}"
        " --fasta {input.genome}"
        " --gtf {input.gtf}"
        " --rlen 150"
        " --threads 16"
        " --use-piscem"
