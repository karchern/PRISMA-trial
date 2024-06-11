library(tidyverse)

#setwd("/Volumes/zeller/karcher/PRISMA/")
setwd("/g/scb/zeller/karcher/PRISMA")

################# CAREFUL ##################
# When comparing WGS and 16S reads, make sure to not join by sampleID!
# You will find overlap but the matching is incorrect!
# Instead always join by c(PSN, Visit) to ensure correct matching!!
############################################

wgs_meta <- read_csv('data/WGS_metadata/230307_PRISMA_Batch1_Metadata.csv', skip = 1) %>%
  mutate(Visit = map_chr(SampleName, function(x) str_split(x, "_")[[1]][length(str_split(x, "_")[[1]])])) %>%
  select(-SampleName)
wgs_number_reads_raw <- readRDS('profiles/WGS/AAAY7KJHV_230111_MG_Batch1/res_libsize_raw_before_QC.rds') %>%
  rename(sampleID = Sample_ID) %>%
  mutate(sampleID = map_chr(sampleID, function(x) str_replace(x, ".raw", "")))
wgs_number_reads_qc <- readRDS('profiles/WGS/AAAY7KJHV_230111_MG_Batch1/res_libsize.rds') %>%
  rename(sampleID = Sample_ID)
wgs_motus <- readRDS('profiles/WGS/AAAY7KJHV_230111_MG_Batch1/res_mOTUs.rds') %>%
  as.data.frame() %>%
  rownames_to_column('taxa') %>%
  pivot_longer(-taxa) %>%
  rename(sampleID = name,
         counts = value) %>%
  group_by(sampleID) %>% 
  mutate(relAb = counts/sum(counts)) %>%
  ungroup()

motu_counts <- wgs_motus %>%
  group_by(sampleID) %>%
  summarize(counts = sum(counts))

wgs_motus <- wgs_motus %>%
  filter(taxa != 'unassigned')

amplicon_counts <- readRDS('profiles/16S/process_KLGPG_221206_PRISMA_MG/res_libsize.rds') %>%
  rename(sampleID = ID,
         mapseq_counts_HQ = libsize) %>%
  arrange(mapseq_counts_HQ) %>%
  mutate(sampleID2 = str_replace(sampleID, ".*MG", "MG_")) %>%
  mutate(isWater = str_detect(sampleID, "water"))

amplicon_meta <- read_csv('data/16S_metadata/221227_PRISMA_Concentrations_Sequencing.csv') %>%
  select(`PSN`, `Visit`, `Sample_ID`, ID) %>%
  rename(sampleID = Sample_ID)
  #rename(PSN = `PSN...1`,
  #       Visit = `Visit...2`,
  #       sampleID = `Sample_ID...4`)

amplicon_profiles <- readRDS('profiles/16S/process_KLGPG_221206_PRISMA_MG/res_mapseq.rds') %>%
  as.data.frame() %>%
  rownames_to_column('genus') %>%
  pivot_longer(-genus) %>% 
  rename(sampleID = name,
       counts = value) %>%
  filter(genus != 'Bacteria') %>%
    group_by(sampleID) %>%
    mutate(relAb = counts / sum(counts)) %>%
    select(-counts) %>%
    filter(genus != 'not_resolved')

# 1. Num. good reads as a function of the concentration values
#motu_counts %>%
  #full_join(wgs_number_reads_raw) %>%
  #full_join(wgs_number_reads_qc) %>% 
p <-  wgs_number_reads_qc %>%
  mutate(sampleID = str_replace(sampleID, ".*MG", "MG_")) %>%
  pivot_longer(-sampleID) %>%
  rename(depth_measures = name,
         depth = value) %>%
  left_join(wgs_meta %>% 
              rename(sampleID = Sample_ID)) %>% print(n=30) %>%
  select(sampleID, depth_measures, depth, `Conc [ng/uL]`, `GeneCore HS Qubit`) %>%
  pivot_longer(-c(sampleID, depth_measures, depth)) %>%
  rename(dna_concentration_value = name,
         dna_concentration = value) %>%
  mutate(depthBases = depth * 100) %>%
  mutate(depthGigaBases = depthBases / 1E9) %>%
  ggplot(aes(y = depthGigaBases, x = dna_concentration)) +
  geom_point() +
  theme_bw() +
  facet_grid(~dna_concentration_value) +
    #scale_x_log10() +
    scale_y_log10() +
    xlab("DNA Concentration") +
    ylab("Sequencing depth\nafter QC (giga bases)")

ggsave(plot = p, filename = 'plots/depth_vs_DNA_concentration.pdf', height = 2.25, width = 4.5)

# 2. 16S seq depth boxplots
p <- amplicon_counts %>%
  ggplot(aes(x = isWater, y = mapseq_counts_HQ)) +
  geom_boxplot() +
  theme_bw() +
  scale_y_log10()
ggsave(plot = p, filename = 'plots/16S_depth.pdf', height = 2.25, width = 3.5)

# 3. Correlation between genus-level profiles of matched smaples between 16S and WGS.
p <- amplicon_profiles %>%
  mutate(sampleID = str_replace(sampleID, ".*MG", "MG_")) %>%
  rename(relab_16S = relAb) %>%
  # Water samples are lost
  inner_join(amplicon_meta, by = 'sampleID') %>%
  inner_join(wgs_meta %>% 
               select(PSN, Visit) %>%
               mutate(Visit = as.numeric(Visit)), by = c("PSN", "Visit")) %>%
  group_by(genus) %>%
  mutate(kick = all(relab_16S == 0)) %>%
  filter(!kick) %>%
  select(-kick) %>%
  ungroup() %>%
  select(PSN, Visit, genus, relab_16S) %>%
  full_join(wgs_motus %>%
               mutate(genus = str_split_fixed(taxa, "[|]", n = 7)[, 6]) %>%
               mutate(genus = str_replace(genus, "g__", "")) %>%
               select(sampleID, genus, relAb) %>%
               mutate(sampleID = str_replace(sampleID, ".*MG", "MG_")) %>%
               left_join(wgs_meta %>% 
                           select(`Sample_ID`, PSN, Visit) %>%
                           rename(sampleID = `Sample_ID`), by = 'sampleID') %>%
               rename(relAb_WGS = relAb) %>%
              select(PSN, Visit, genus, relAb_WGS) %>%
              mutate(Visit = as.numeric(Visit)) %>%
              group_by(Visit, PSN, genus) %>%
              summarize(relAb_WGS = sum(relAb_WGS)),
             by = c('PSN', 'Visit', "genus")) %>%
  rename(relAb_16S = relab_16S) %>%
  mutate(relAb_16S = ifelse(is.na(relAb_16S), 0, relAb_16S)) %>%
  mutate(relAb_WGS = ifelse(is.na(relAb_WGS), 0, relAb_WGS)) %>%
  ggplot(aes(x = relAb_16S+1E-6, y = relAb_WGS+1E-6)) +
  geom_point(alpha = 0.2) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw()

ggsave(plot = p, filename = "plots/scatter_wgs_16s.pdf", width = 6, height = 6)

