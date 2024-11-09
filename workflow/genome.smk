rule prepare_genome_fasta:
    output:
        "resources/genome/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.fa.gz",
    log:
        "logs/prepare_genome_fasta.log",
    shell:
        "mkdir -p tmp_genome && "
        "curl https://ftp.ensembl.org/pub/release-113/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.2L.fa.gz > tmp_genome/2L.fa.gz 2> {log} && "
        "curl https://ftp.ensembl.org/pub/release-113/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.2R.fa.gz > tmp_genome/2R.fa.gz 2> {log} && "
        "curl https://ftp.ensembl.org/pub/release-113/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.3L.fa.gz > tmp_genome/3L.fa.gz 2> {log} && "
        "curl https://ftp.ensembl.org/pub/release-113/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.3R.fa.gz > tmp_genome/3R.fa.gz 2> {log} && "
        "curl https://ftp.ensembl.org/pub/release-113/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.4.fa.gz > tmp_genome/4.fa.gz 2> {log} && "
        "curl https://ftp.ensembl.org/pub/release-113/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.X.fa.gz > tmp_genome/X.fa.gz 2> {log} && "
        "curl https://ftp.ensembl.org/pub/release-113/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.Y.fa.gz > tmp_genome/Y.fa.gz 2> {log} && "
        "curl https://ftp.ensembl.org/pub/release-113/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.mitochondrion_genome.fa.gz 2> {log} > tmp_genome/mitochondrion_genome.fa.gz && "
        "zcat tmp_genome/*.fa.gz > resources/genome/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.fa 2> {log} && "
        "bgzip resources/genome/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.fa 2> {log} && "
        "rm -rf tmp_genome"

rule prepare_genome_gatk_index:
    input:
        "resources/genome/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.fa.gz",
    output:
        "resources/genome/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.dict",
    log:
        gatkout="logs/prepare_genome_gatk_index.gatk.out",
        gatkerr="logs/prepare_genome_gatk_index.gatk.err",
        samtoolkout="logs/prepare_genome_gatk_index.samtools.out",
        samtoolserr="logs/prepare_genome_gatk_index.samtools.err",
    shell:
        "gatk CreateSequenceDictionary -R {input} > {log.gatkout} 2> {log.gatkout} && "
        "samtools faidx {input} > {log.samtoolkout} 2> {log.samtoolserr}"
