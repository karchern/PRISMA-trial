library(tidyverse)
library(ggplot2)
library(patchwork)
library(vegan)
library(ggrepel)
library(RColorBrewer)
library(eeptools)
library(readxl)
library(here)
library(ggembl)
library(ggquantileplot)
# source('/home/karcher/utils/utils.r')
source(here('scripts/utils.r'))

# Load data
obj_path <- here('objects/PRISMA.rdata')
load_data(obj_path)


# tmp <- pairwiseDistances %>%
#     as.matrix()
# tmp[upper.tri(tmp)] <- NA
# tmp <- tmp %>%
#     as.data.frame() %>%
#     rownames_to_column('tmp') %>%
#     pivot_longer(-tmp) %>%
#     mutate(sampleID_1 = str_split_fixed(tmp, "___", n = 2)[, 1]) %>%
#     mutate(batch_1 = str_split_fixed(tmp, "___", n = 2)[, 2]) %>%
#     mutate(sampleID_2 = str_split_fixed(name, "___", n = 2)[, 1]) %>%
#     mutate(batch_2 = str_split_fixed(name, "___", n = 2)[, 2]) %>%
#     rename(distance = value) %>%
#     filter(!is.na(distance)) %>%
#     left_join(meta %>%
#         select(sampleID, PSN, visit), by = c("sampleID_1" = "sampleID")) %>%
#     left_join(meta %>%
#         select(sampleID, PSN, visit), by = c("sampleID_2" = "sampleID"), suffix = c("_1", "_2")) %>%
#     mutate(visit_1 = as.factor(visit_1)) %>%
#     mutate(visit_2 = as.factor(visit_2)) %>%
#     filter(PSN_1 == PSN_2) %>%
#     rename(PSN = PSN_1) %>%
#     select(-PSN_2) %>%
#     filter(visit_1 != visit_2) %>%
#     as_tibble()

for (taxLevel in c(
    # "phylum",
    # "class",
    # "order",
    "family"
    # "genus"
)
) {
    set.seed(2)
    taxL <- taxLevel
    numTaxa <- 10
    levelsToShow <- profiles %>%
        mutate(relAb = 10^relAb - pseudoCount) %>%
        group_by(sampleID, .data[[taxL]]) %>%
        summarize(relAb = sum(relAb)) %>%
        group_by(.data[[taxL]]) %>%
        summarize(m = mean(relAb)) %>%
        arrange(desc(m)) %>%
        head(numTaxa) %>%
        pull(.data[[taxL]])
    getPalette <- colorRampPalette(brewer.pal(9, "Paired"))
    colors <- sample(getPalette(numTaxa))
    colors <- c(colors, "#808080", "#D3D3D3")

    ggsave(
        plot = get_family_level_barplot_for_all_samples_WIDE(ggplot(),
            profiles %>%
                # filter(visit != 3) %>%
                mutate(visit = as.numeric(as.character(visit))),
            taxL,
            levelsToShow = levelsToShow) +
            geom_rect(
                data = outcomeInformation %>%
                    filter(CD > 2) %>%
                    rename(PSN = patientID, visit = visitNumber) %>%
                    mutate(PSN = factor(PSN, levels = orderDFPatientID$PSN)) %>%
                    mutate(PSN = as.numeric(PSN)),
                aes(xmin = PSN - 0.5, xmax = PSN + 0.5, ymin = 0, ymax = 1), color = 'black', fill = alpha("grey", 0),
                inherit.aes = F) +
            scale_fill_manual(values = colors) +
            NULL,
        # ggtitle("red dots indicate clinical complication\nblack asterisks indicate ABx intake (other than prophylactic)"),
        filename = str_c(here("plots/KLGPG_221206/taxa_barplots_visit_v3_visit_3_removed_treatment_dichotomized_all__WIDE_VERSION_"), taxL, ".pdf", collapse = ""),
        width = 8,
        height = 10
    )
}
