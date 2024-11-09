rule prepare_genome_fasta:
    output:
        "resources/genome/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.main.fa.gz"
    shell:
        "mkdir -p tmp_genome &&"
        "curl https://ftp.ensembl.org/pub/release-113/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.2L.fa.gz > tmp_genome/2L.fa.gz &&"
        "curl https://ftp.ensembl.org/pub/release-113/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.2R.fa.gz > tmp_genome/2R.fa.gz &&"
        "curl https://ftp.ensembl.org/pub/release-113/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.3L.fa.gz > tmp_genome/3L.fa.gz &&"
        "curl https://ftp.ensembl.org/pub/release-113/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.3R.fa.gz > tmp_genome/3R.fa.gz &&"
        "curl https://ftp.ensembl.org/pub/release-113/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.4.fa.gz > tmp_genome/4.fa.gz &&"
        "curl https://ftp.ensembl.org/pub/release-113/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.X.fa.gz > tmp_genome/X.fa.gz &&"
        "curl https://ftp.ensembl.org/pub/release-113/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.Y.fa.gz > tmp_genome/Y.fa.gz &&"
        "curl https://ftp.ensembl.org/pub/release-113/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.mitochondrion_genome.fa.gz > tmp_genome/mitochondrion_genome.fa.gz &&"
        "zcat tmp_genome/*.fa.gz > resources/genome/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.main.fa.gz &&"
        "rm -rf tmp_genome"

