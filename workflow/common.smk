def get_final_output():
    intervals = config['genome']['intervals']
    final_output = expand(
        "resources/genome_sequencing/mapped.dedup.sorted.{interval}.gvcf",
        interval=intervals,
    )
    final_output.append(
        expand(
            "resources/genome_sequencing/mapped.dedup.sorted.{interval}.bamout.bam",
            interval=intervals,
        )
    )
    return final_output
