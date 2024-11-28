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


#####
# taxonomy_annot is not being accessed if WGS data is loaded!
# If you indeed use WGS data, make sure to set this to ""
taxonomy_annot <- ""
# taxonomy_annot <- "ncbi_mapseq"
# taxonomy_annot <- "gtdb_idtaxa"

if (!taxonomy_annot %in% c("ncbi_mapseq", "gtdb_idtaxa", "")) {
    if (taxonomy_annot == "") {
        print("No taxonomy annotation specified, this means youre working with WGS data")
    }
    stop("Unknown taxonomy annotation")
}

#obj_path <- here(str_c('objects/PRISMA_', taxonomy_annot, '.rdata'))
obj_path <- here(str_c('objects/PRISMA_WGS', taxonomy_annot, '.rdata'))
load_data(obj_path)
obj_path <- here(str_c('objects/PRISMA_', "ncbi_mapseq", '.rdata'))
load_data(obj_path)

# Get combined distance matrix between WGS and 16S
wgs <- profiles_wgs_genus %>%
    select(sampleID, visit, genus, relAb) %>%
    mutate(type = 'wgs') %>%
    inner_join(meta %>% select(PSN, visit, batch), by = c("PSN", 'visit'))
sixteens <- profiles %>%
    select(sampleID, visit, genus, relAb) %>%
    mutate(type = '16S') %>%
    inner_join(meta %>% select(PSN, visit, batch), by = c("PSN", 'visit'))
combined <- bind_rows(wgs, sixteens)

shared_genera <- combined %>%
        group_by(genus, type) %>% 
        summarize(m = mean(relAb)) %>%
        pivot_wider(names_from = type, values_from = m) %>%
        filter(!(is.na(`16S`))) %>%
        filter(!is.na(wgs)) %>%
        select(genus)

# Get PCOA
pairwiseDistancesGenusAll <- combined %>%
    mutate(s = str_c(PSN, visit, sampleID, type, batch, sep = "___")) %>%
    ungroup() %>%
    select(s, genus, relAb) %>%
    inner_join(shared_genera) %>%
    pivot_wider(names_from = genus, values_from = relAb) %>%
    as.data.frame() %>%
    column_to_rownames('s') %>%
    dist(method = 'euclidean')


# Get Permanova
pcoa <- cmdscale(pairwiseDistancesGenusAll, k = 2)
pcoa <- pcoa %>%
    as.data.frame() %>%
    mutate(
        PSN = str_split(rownames(.), "___") %>% map(\(x) x[1]) %>% unlist(),
        visit = str_split(rownames(.), "___") %>% map(\(x) x[2]) %>% unlist(),
        sampleID = str_split(rownames(.), "___") %>% map(\(x) x[3]) %>% unlist(),
        type = str_split(rownames(.), "___") %>% map(\(x) x[4]) %>% unlist(),
        batch = str_split(rownames(.), "___") %>% map(\(x) x[5]) %>% unlist()
        ) %>%
        as.tibble()

pcoa_plot <- ggplot() +
    geom_point(data = pcoa, aes(x = V1, y = V2, color = type)) +
    theme_presentation() +
    xlab("PCo 1") +
    ylab("PCo 2")

ggsave(
    plot = pcoa_plot,
    filename =  here("plots/KLGPG_221206/pcoa_WGS_vs_16S_v1.pdf"),
    width = 6, height = 6
)
