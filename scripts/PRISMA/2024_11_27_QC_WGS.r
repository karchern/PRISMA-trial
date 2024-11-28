library(here)
library(tidyverse)
library(ggembl)

# depth
p <- read_tsv(here("data/fastqc_sequence_counts_plot.tsv"), col_names = TRUE) %>%
    mutate(batch_raw = str_split_fixed(Sample, " [|] ", n = 6)[, 5]) %>%
    mutate(
        batch = str_replace(batch_raw, "_for_Q", ""),
        batch = str_replace(batch,  ".*_PRISMA_", "")
        ) %>% 
        select(-batch_raw, -Sample) %>%
        mutate(num_bases = `Unique Reads` * 150) %>%
        mutate(giga_bases = num_bases / 1e9) %>%
        select(batch, giga_bases) %>%
        mutate(batch = factor(batch, levels = sort(unique(batch)))) %>%
    ggplot(aes(x = batch, y = giga_bases)) +
    geom_boxplot() +
    theme_presentation() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ylab("Sequencing depth (giga-base pairs)") +
    xlab("WGS Sequencing batch")

ggsave(
    plot = p,
    filename = here("plots/KLGPG_221206/WGS_QC_read_depth.pdf"),
    width = 6, height = 4
)

# Read quality
data <- read_tsv(here("data/fastqc_per_base_sequence_quality_plot.tsv"), col_names = FALSE)
# Separate X and Y data
x_data <- data %>% filter(X2 == "X") %>% select(-X2) %>% rename(s = X1)
y_data <- data %>% filter(X2 == "Y") %>% select(-X2) %>% rename(s = X1)

x_sc <- as.data.frame(x_data[1, ])[1, ][2:40]
x_sc <- x_sc %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column('x') %>%
    rename(x_actual = `1`)

longo <- y_data %>%
  pivot_longer(cols = starts_with("X"), names_to = "x", values_to = "value") %>%
  left_join(x_sc) %>%
  rename(Sample = s) %>%
  mutate(batch_raw = str_split_fixed(Sample, " [|] ", n = 6)[, 5]) %>%
    mutate(
        batch = str_replace(batch_raw, "_for_Q", ""),
        batch = str_replace(batch,  ".*_PRISMA_", "")
        ) %>% 
        select(-batch_raw)
  

p <- ggplot(
    data = longo
) +
    geom_line(aes(x = x_actual, y = value, group = Sample), alpha = 0.2) +
    theme_publication() +
    facet_wrap(.~batch,ncol = 1) +
    xlab("Read position") +
    ylab("Mean PHRED score")
ggsave(
    plot = p,
    filename = here("plots/KLGPG_221206/WGS_QC_read_qual.pdf"),
    width = 2.5, height = 4.5
)
