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
# source('/home/karcher/utils/utils.r')
source(here('scripts/utils.r'))

# obj_path <- here('objects/PRISMA_idtaxa.rdata')
obj_path <- here('objects/PRISMA.rdata')
load_data(obj_path)

pcoa_plot <- ggplot() +
    geom_point(data = pcoa, aes(x = V1, y = V2, color = visit)) +
    theme_presentation() +
    xlab("PCo 1") +
    ylab("PCo 2") +
    scale_color_manual(values = time_point_colors)

ggsave(plot = pcoa_plot, filename = here("plots/KLGPG_221206/pcoa_visit_v1.pdf"), width = 5.5, height = 3.5)

pcoa_plot <- ggplot() +
    geom_point(data =
        pcoa %>%
            filter(visit != 3) %>%
            mutate(visit = factor(map_chr(visit, \(x) labelLink[x]), levels = labelLink)), aes(x = V1, y = V2, color = visit)) +
    theme_presentation() +
    xlab("PCo 1") +
    ylab("PCo 2") +
    scale_color_manual(values = time_point_colors[c(1:2, 4:7)])

ggsave(plot = pcoa_plot, filename = here("plots/KLGPG_221206/pcoa_visit_v2.pdf"), width = 8, height = 3.5)



# Compute a PERMANOVA in distance space using batch
stopifnot(all(pairwiseDistances %>% as.matrix() %>% rownames() %>% map(\(x) str_split(x, "___")[[1]][1]) %>% unlist() == pcoa$sampleID))
pcoa_perm <- adonis2(pairwiseDistances ~ batch, data = pcoa, permutations = 999)
pcoa_plot <- ggplot() +
    geom_point(data = pcoa, aes(x = V1, y = V2, color = batch)) +
    theme_presentation() +
    xlab("PCo 1") +
    ylab("PCo 2") +
    scale_color_manual(values = batch_colors) +
    annotate('text', x = min(pcoa$V1) + (0.0 * abs(max(pcoa$V2) - min(pcoa$V2))), y = min(pcoa$V2) + (0.025 * abs(max(pcoa$V1) - min(pcoa$V1))), label = str_c("PERMANOVA p-value: ", round(pcoa_perm$`Pr(>F)`[which(rownames(pcoa_perm) == "batch")], 3), "\nPERMANOVA R2: ", round(pcoa_perm$`R2`[which(rownames(pcoa_perm) == "batch")], 3)), hjust = 0)

ggsave(plot = pcoa_plot, filename = here("plots/KLGPG_221206/pcoa_batch.pdf"), width = 5.5, height = 3.5)

pcoa_plot <- ggplot() +
    geom_point(data = pcoa %>%
        left_join(clinicalMetadata %>%
            select(patientID, studyCenter) %>% distinct(), by = c("PSN" = 'patientID')) %>%
        mutate(studyCenter = case_when(
            studyCenter == "k" ~ "Kinderklinik",
            studyCenter == "m" ~ "Muenster",
            studyCenter == "n" ~ "Nierenzentrum",
        )), aes(x = V1, y = V2, color = studyCenter)) +
    theme_presentation() +
    xlab("PCo 1") +
    ylab("PCo 2") +
    NULL

ggsave(plot = pcoa_plot, filename = here("plots/KLGPG_221206/pcoa_visit_study_center.pdf"), width = 8, height = 5)

pcoa_plot <- ggplot() +
    geom_point(data = pcoa %>%
        left_join(clinicalMetadata %>%
            select(patientID, age) %>% distinct(), by = c("PSN" = 'patientID')) %>%
        mutate(age = as.numeric(age)) %>%
        mutate(age = case_when(
            age < 6 ~ "< 6 years",
            age >= 6 & age < 12 ~ ">= 6, < 12 years",
            age >= 12 ~ "> 12 years",
        )) %>%
        mutate(age = factor(age, levels = c("< 6 years", ">= 6, < 12 years", "> 12 years"))), aes(x = V1, y = V2, color = age)) +
    theme_presentation() +
    xlab("PCo 1") +
    ylab("PCo 2") +
    NULL

ggsave(plot = pcoa_plot, filename = here("plots/KLGPG_221206/pcoa_visit_age.pdf"), width = 8, height = 5)
for (indName in unique(pcoa$PSN)) {
    print(indName)
    # for (indName in c('RiGr-NZHD-8')) {
    ggsave(
        plot = add_line_with_ind(ggplot() + geom_point(data = pcoa, aes(x = V1, y = V2, color = visit)), pcoa, indName) +
            theme_presentation() +
            # theme(legend.title = element_text(size = 0), legend.text = element_text(size = 12)) +
            scale_color_manual(values = c(
                "#5b92fc",
                "#5b92fc",
                "#3579fc",
                "#ff3636",
                "#f75e5e",
                "#facaca",
                "#fcebeb")) +
            xlab("PCo 1") +
            ylab("PCo 2"),
        filename = str_c(here("plots/KLGPG_221206/single_ind_pcoa_plots/pcoa_visit_v3_visit_3_removed_treatment_dichotomized_"), indName, ".pdf", collapse = ""), width = 5.5, height = 3.5)
}

nice_looking_individuals <- c(
    # These are from interims cohort
    "AbWa-KKHD-12",
    "RaCu-NZHD-17",
    "BaEr-NZHD-39",
    "DeZi-KKHD-2",
    "AnCz-KKHD-11",
    'MaBa-NZHD-10',
    # ... and these are from the modelling cohort
    'AlGoe-NZHD-28',
    "DaBe-NZMU-21"
)

get_family_level_barplot <- function(pObj, dataB, indName, taxLevel = 'family', levelsToShow = NULL) {
    dataB <- dataB %>%
        filter(PSN == indName) %>%
        mutate(relAb = (10^relAb) - pseudoCount) %>%
        mutate(visit = as.factor(visit)) %>%
        mutate(taxa = .data[[taxLevel]]) %>%
        mutate(taxa = as.character(taxa)) %>%
        mutate(taxa = ifelse(taxa %in% levelsToShow, taxa, "other")) %>%
        # mutate(taxa = factor(taxa, levels = c(levelsToShow, "other", 'unclassified'))) %>%
        group_by(taxa, sampleID, visit) %>%
        summarize(relAb = sum(relAb))
    dataC <- dataB %>%
        ungroup() %>%
        group_by(sampleID, visit) %>%
        summarize(relAb = 1 - sum(relAb)) %>%
        mutate(taxa = "unclassified")
    dataB <- rbind(dataB, dataC) %>%
        mutate(taxa = factor(as.character(taxa), levels = c(levelsToShow, "other", 'unclassified')))

    pObj <- pObj +
        geom_bar(data = dataB,
            aes(x = visit, y = relAb, fill = taxa), position = 'stack', stat = 'identity') +
        theme_presentation()
    return(pObj)
}

for (taxL in c(
    # "phylum",
    # "class",
    # "order",
    "family"
    # "genus"
)
) {
    set.seed(2)
    # taxL <- "family"
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

    # for (indName in unique(pcoa$PSN)) {
    # for (indName in c("AbWa-KKHD-12")) {
    for (indName in nice_looking_individuals) {
        if (indName %in% profiles$PSN) {
            print(str_c("Producing family-level barplot over time for ", indName))
            ggsave(
                plot = get_family_level_barplot(ggplot(),
                    profiles %>%
                        filter(visit != 3) %>%
                        mutate(visit = as.numeric(as.character(visit))),
                    indName,
                    taxL,
                    levelsToShow = levelsToShow) +
                    scale_fill_manual(values = colors) + scale_x_discrete_prisma() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                    guides(fill = guide_legend(ncol = 2)),
                filename = str_c(here("plots/KLGPG_221206/single_ind_tax_barplots/taxa_barplots_visit_v3_visit_3_removed_treatment_dichotomized_"), indName, "__", taxL, ".pdf", collapse = ""),
                width = 5.75,
                height = 3
            )
        } else {
            print(str_c("Skipping family-level barplot over time for ", indName, " cause individual not found in profiles..."))
        }

    }
}
