def get_final_output():
    intervals = config['genome']['intervals']
    final_output = expand(
        "genome_sequencing/aln-pe.rg.dedup.sorted.${interval}.gvcf",
        interval=intervals,
    )
    final_output.append(
        expand(
            "genome_sequencing/aln-pe.rg.dedup.sorted.${interval}.bamout.bam",
            interval=intervals,
        )
    )
    return final_output