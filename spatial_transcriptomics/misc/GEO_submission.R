# R script to consolidate sample FASTQ information
library(dplyr)
library(stringr)
library(readr)
library(here)
library(glue)

source(here("misc/paths.R"))

id <- "01-spaceranger/data/sample_id.txt" %>%
    here() %>%
    read_csv()

id <- id %>% select(-c(subproject, index))

info_32 <- "01-spaceranger/projects/BCLLATLAS_32/info.txt" %>%
    here() %>%
    read_tsv()

info_51 <- "01-spaceranger/projects/BCLLATLAS_51/info.txt" %>%
    here() %>%
    read_tsv()

path <- "/scratch/project/production/fastq"

info <- bind_rows(info_32, info_51) %>%
    left_join(id, by = c("library" = "library_id"))

info_long <- info %>%
    mutate(
        P1 = glue("{path}/{info$flowcell}/1/fastq/{info$flowcell}_{info$lane}_{info$index}_{info$lane}.fastq.gz"),
        P2 = glue("{path}/{info$flowcell}/2/fastq/{info$flowcell}_{info$lane}_{info$index}_{info$lane}.fastq.gz")) %>% 
    tidyr::pivot_longer(c(P1, P2),
        names_to = "pair_id",
        values_to = "fastq_path")

info_subset <- info_long %>%
    dplyr::select(c(
        subproject, gem_id, library, type, donor_id, index, slide, area,
        lane, pair_id, fastq_path)) %>%
    rename(
        library_id = library,
        read = lane
    )

here("misc/fastq_paths.tsv") %>%
    write_tsv(x = info_subset, file = .)


## Images
img_path <- "/scratch/devel/melosua/phd/projects/BCLLatlas/tonsil_atlas/spatial_transcriptomics/01-spaceranger/img"
img_df <- lapply(unique(info_subset$gem_id), function(i) {
    slide <- info_subset %>% filter(gem_id == i) %>% pull(slide) %>% unique()
    area <- info_subset %>% filter(gem_id == i) %>% pull(area) %>% unique()
    p <- glue("{img_path}/{i}_{slide}_{area}.jpg")
    data.frame(gem_id = i, img_path = p)
}) %>%
    bind_rows()

here("misc/img_paths.tsv") %>%
    write_tsv(x = img_df, file = .)


