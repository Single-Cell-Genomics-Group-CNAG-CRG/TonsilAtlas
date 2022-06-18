# This script plots the results of the label transfer for both cases (102, 413)


# Load packages
library(tidyverse)


# Read data
df_102 <- readRDS(here::here("MCL/results/R_objects/dataframe_CD4_T_label_transfer_102.rds"))
df_413 <- readRDS(here::here("MCL/results/R_objects/dataframe_CD4_T_label_transfer_413.rds"))


# Plot UMAP
annotation_df <- bind_rows(df_102, df_413)
annotation_df$type <- case_when(
  annotation_df$type == "reference" ~ "reference",
  annotation_df$type == "query" & annotation_df$donor_id == "M102" ~ "M102",
  annotation_df$type == "query" & annotation_df$donor_id == "M413" ~ "M413",  
)
reordered_levels <- c("Naive", "CM PreTfh", "CM Pre-non-Tfh", "Tfh T:B border",
                      "Tfh-LZ-GC", "GC-Tfh-0X40", "GC-Tfh-SAP", "Tfh-Mem",
                      "T-Trans-Mem", "T-Eff-Mem", "T-helper", "GC-Tf-regs",
                      "non-GC-Tf-regs", "Eff-Tregs")
annotation_df$label <- factor(annotation_df$label, levels = reordered_levels) 
umap_annotation <- annotation_df %>%
  mutate(type = factor(type, levels = c("reference", "M102", "M413"))) %>%
  ggplot(aes(UMAP1, UMAP2, color = label)) +
    geom_point(size = 0.35) +
    facet_wrap(~type) +
    scale_color_manual(values = pals::glasbey()) +
    theme_classic() +
    theme(legend.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 2)))


# Plot stacked barplot
adult_ids <- c("BCLL-2-T", "BCLL-6-T", "BCLL-14-T", "BCLL-15-T")
annotation_df$type2 <- case_when(
  annotation_df$type == "M102" ~ "M102",
  annotation_df$type == "M413" ~ "M413",
  annotation_df$type == "reference" & annotation_df$donor_id %in% adult_ids ~ "adult",
  annotation_df$type == "reference" & !(annotation_df$donor_id %in% adult_ids) ~ "kid"
)
annotation_df$type2 <- factor(
  annotation_df$type2,
  levels = c("kid", "adult", "M102", "M413")
)
proportion_df <- annotation_df %>%
  group_by(type2, label) %>%
  summarize(n_cells = n()) %>%
  group_by(type2) %>%
  mutate(total_cells = sum(n_cells), pct_cells = n_cells / total_cells * 100)
proportion_gg <- proportion_df %>%
  ggplot(aes(type2, pct_cells, fill = label)) +
    geom_col() +
    scale_fill_manual(values = pals::glasbey()) +
    ylab("Percentage of cells (%)") +
    theme_classic() +
    theme(axis.title.x = element_blank(), legend.title = element_blank())


# Save
ggsave(
  plot = umap_annotation,
  filename = here::here("MCL/results/plots/umap_CD4_T_label_transfer.png"),
  width = 21,
  height = 11.8,
  units = "cm"
)
ggsave(
  plot = proportion_gg,
  filename = here::here("MCL/results/plots/stacked_barplot_CD4_T_label_transfer.png"),
  width = 21,
  height = 11.8,
  units = "cm"
)
