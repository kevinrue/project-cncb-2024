library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(stringr)

# Raw input ----

# Manually typed list based on email from Aaron (29-Oct-2024)
samples_list <- list(
  list(
    sample_name="WPPm048hrs_rep1",
    fastqs="/project/cncb/shared/proj140/analyses/novogene_sequencing/single_cell/20240925_run_2/fastq/X204SC24073110-Z01-F003_002/01.RawData/Sample_1/",
    sequencing_sample="Sample_1-SCI7T003-SCI5T003_22FVLJLT4,Sample_1-SCI7T003-SCI5T003_22FYLKLT4,Sample_1-SCI7T003-SCI5T003_22G25CLT4,Sample_1-SCI7T003-SCI5T003_22GF2KLT4"
  ),
  list(
    sample_name="WPPm024hrs_rep1",
    fastqs="/project/cncb/shared/proj140/analyses/novogene_sequencing/single_cell/20240925_run_2/fastq/X204SC24073110-Z01-F003_010/01.RawData/Sample_2/",
    sequencing_sample="Sample_2-SCI7T015-SCI5T015_22FVLJLT4,Sample_2-SCI7T015-SCI5T015_22GF2KLT4"
  ),
  list(
    sample_name="WPPp000hrs_rep1",
    fastqs="/project/cncb/shared/proj140/analyses/novogene_sequencing/single_cell/20240925_run_2/fastq/X204SC24073110-Z01-F003_011/01.RawData/Sample_3/",
    sequencing_sample="Sample_3-SCI7T027-SCI5T027_22FVLJLT4,Sample_3-SCI7T027-SCI5T027_22GF2KLT4"
  ),
  list(
    sample_name="WPPp024hrs_rep1",
    fastqs="/project/cncb/shared/proj140/analyses/novogene_sequencing/single_cell/20240925_run_2/fastq/X204SC24073110-Z01-F003_012/01.RawData/Sample_4/",
    sequencing_sample="Sample_4-SCI7T039-SCI5T039_22FVLJLT4,Sample_4-SCI7T039-SCI5T039_22GF2KLT4"
  ),
  list(
    sample_name="WPPp048hrs_rep1",
    fastqs="/project/cncb/shared/proj140/analyses/novogene_sequencing/single_cell/20240925_run_2/fastq/X204SC24073110-Z01-F003_013/01.RawData/Sample_5/",
    sequencing_sample="Sample_5-SCI7T051-SCI5T051_22FVLJLT4,Sample_5-SCI7T051-SCI5T051_22FYLKLT4,Sample_5-SCI7T051-SCI5T051_22G25CLT4,Sample_5-SCI7T051-SCI5T051_22GF2KLT4"
  ),
  list(
    sample_name="WPPp072hrs_rep1",
    fastqs="/project/cncb/shared/proj140/analyses/novogene_sequencing/single_cell/20240925_run_2/fastq/X204SC24073110-Z01-F003_014/01.RawData/Sample_6/",
    sequencing_sample="Sample_6-SCI7T063-SCI5T063_22FVLJLT4,Sample_6-SCI7T063-SCI5T063_22GF2KLT4"
  ),
  list(
    sample_name="WPPp096hrs_rep1",
    fastqs="/project/cncb/shared/proj140/analyses/novogene_sequencing/single_cell/20240925_run_2/fastq/X204SC24073110-Z01-F003_015/01.RawData/Sample_7/",
    sequencing_sample="Sample_7-SCI7T075-SCI5T075_22FVLJLT4,Sample_7-SCI7T075-SCI5T075_22GF2KLT4"
  ),
  list(
    sample_name="WPPp120hrs_rep1",
    fastqs="/project/cncb/shared/proj140/analyses/novogene_sequencing/single_cell/20240925_run_2/fastq/X204SC24073110-Z01-F003_016/01.RawData/Sample_8/",
    sequencing_sample="Sample_8-SCI7T087-SCI5T087_22FVLJLT4,Sample_8-SCI7T087-SCI5T087_22GF2KLT4"
  ),
  list(
    sample_name="WPPm048hrs_rep2",
    fastqs="/project/cncb/shared/proj140/analyses/novogene_sequencing/single_cell/20240925_run_2/fastq/X204SC24073110-Z01-F003_017/01.RawData/Sample_9/",
    sequencing_sample="Sample_9-SCI7T004-SCI5T004_22FVLJLT4,Sample_9-SCI7T004-SCI5T004_22GF2KLT4"
  ),
  list(
    sample_name="WPPm024hrs_rep2",
    fastqs="/project/cncb/shared/proj140/analyses/novogene_sequencing/single_cell/20240925_run_2/fastq/X204SC24073110-Z01-F003_003/01.RawData/Sample_10/",
    sequencing_sample="Sample_10-SCI7T016-SCI5T016_22FVLJLT4,Sample_10-SCI7T016-SCI5T016_22FYLKLT4,Sample_10-SCI7T016-SCI5T016_22GF2KLT4"
  ),
  list(
    sample_name="WPPp000hrs_rep2",
    fastqs="/project/cncb/shared/proj140/analyses/novogene_sequencing/single_cell/20240925_run_2/fastq/X204SC24073110-Z01-F003_004/01.RawData/Sample_11/",
    sequencing_sample="Sample_11-SCI7T028-SCI5T028_22FVLJLT4,Sample_11-SCI7T028-SCI5T028_22FYLKLT4,Sample_11-SCI7T028-SCI5T028_22G25CLT4,Sample_11-SCI7T028-SCI5T028_22GF2KLT4"
  ),
  list(
    sample_name="WPPp024hrs_rep2",
    fastqs="/project/cncb/shared/proj140/analyses/novogene_sequencing/single_cell/20240925_run_2/fastq/X204SC24073110-Z01-F003_005/01.RawData/Sample_12/",
    sequencing_sample="Sample_12-SCI7T040-SCI5T040_22FVLJLT4,Sample_12-SCI7T040-SCI5T040_22FYLKLT4,Sample_12-SCI7T040-SCI5T040_22GF2KLT4"
  ),
  list(
    sample_name="WPPp048hrs_rep2",
    fastqs="/project/cncb/shared/proj140/analyses/novogene_sequencing/single_cell/20240925_run_2/fastq/X204SC24073110-Z01-F003_006/01.RawData/Sample_13/",
    sequencing_sample="Sample_13-SCI7T052-SCI5T052_22FVLJLT4,Sample_13-SCI7T052-SCI5T052_22FYLKLT4,Sample_13-SCI7T052-SCI5T052_22G25CLT4,Sample_13-SCI7T052-SCI5T052_22GF2KLT4"
  ),
  list(
    sample_name="WPPp072hrs_rep2",
    fastqs="/project/cncb/shared/proj140/analyses/novogene_sequencing/single_cell/20240925_run_2/fastq/X204SC24073110-Z01-F003_007/01.RawData/Sample_14/",
    sequencing_sample="Sample_14-SCI7T064-SCI5T064_22FVLJLT4,Sample_14-SCI7T064-SCI5T064_22GF2KLT4"
  ),
  list(
    sample_name="WPPp096hrs_rep2",
    fastqs="/project/cncb/shared/proj140/analyses/novogene_sequencing/single_cell/20240925_run_2/fastq/X204SC24073110-Z01-F003_008/01.RawData/Sample_15/",
    sequencing_sample="Sample_15-SCI7T076-SCI5T076_22FVLJLT4,Sample_15-SCI7T076-SCI5T076_22FYLKLT4,Sample_15-SCI7T076-SCI5T076_22GF2KLT4"
  ),
  list(
    sample_name="WPPp120hrs_rep2",
    fastqs="/project/cncb/shared/proj140/analyses/novogene_sequencing/single_cell/20240925_run_2/fastq/X204SC24073110-Z01-F003_009/01.RawData/Sample_16/",
    sequencing_sample="Sample_16-SCI7T088-SCI5T088_22FVLJLT4,Sample_16-SCI7T088-SCI5T088_22FYLKLT4,Sample_16-SCI7T088-SCI5T088_22G25CLT4,Sample_16-SCI7T088-SCI5T088_22GF2KLT4"
  )
)

# Functions ----

process_sample <- function(sample_info) {
  fastqs_files <- tibble(
    basename = list.files(path = sample_info[["fastqs"]], pattern = ".fastq.gz")
  ) %>%
    mutate(
      directory = sample_info[["fastqs"]],
      sample_name = sample_info[["sample_name"]],
      sequencing_sample = str_extract(basename, "(Sample_[[:digit:]]+-[[:alnum:]]{8}-[[:alnum:]]{8}_[[:alnum:]]{9})", group = 1),
      type = str_extract(basename, "_(I1|I2|R1|R2)_", group = 1),
      lane = str_extract(basename, "_(L[[:digit:]]{3})_", group = 1)
    )
  stopifnot(all(fastqs_files$sequencing_sample %in% str_split(string = sample_info[["sequencing_sample"]], pattern = ",")[[1]]))
  fastqs_files <- fastqs_files %>%
    pivot_wider(values_from = basename, names_from = type) %>%
    write_tsv(file = "sample1.tsv")
  return(fastqs_files)
}

# Execution ----

sample_fastq_table <- do.call("rbind", lapply(X = samples_list, FUN = process_sample))

# Output ----

write_tsv(x = sample_fastq_table, file = "../snakemake-cncb-2024/config/samples.tsv")
