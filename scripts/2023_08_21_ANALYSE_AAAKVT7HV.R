library(tidyverse)
library(ggplot2)
library(patchwork)
library(readxl)
library(vegan)
library(RColorBrewer)
library(taxonomizr)
library(ape)
library(ggunileg)
library(ggrepel)
library(ggembl)

reloadHealthyProfiles <- TRUE


setwd('/g/scb/zeller/karcher/PRISMA')

theme_paper <- ggembl::theme_presentation() +
    theme(
        axis.title = element_text(face = "bold", size = 12), # Bold axis titles
        panel.border = element_rect(fill=NA, colour='black', size=1.5),
        #axis.line = element_line(size = 1.5), # Change line thickness of axes
        axis.text = element_text(face = "bold", size = 12),
    )

leg_inside <-
    theme(
        legend.title = element_blank(),
        legend.position = c(-0, 0),
        # legend.background = element_rect(color = "black", size = 0.5),
        legend.justification = c(0, 0),
        legend.margin = margin(6, 6, 6, 6))

leg_inside_tl <-
    theme(
        legend.title = element_blank(),
        legend.position = c(-0, 1),
        # legend.background = element_rect(color = "black", size = 0.5),
        legend.justification = c(0, 1),
        legend.margin = margin(6, 6, 6, 6))
leg_inside_tr <-
    theme(
        legend.title = element_blank(),
        legend.position = c(1, 1),
        # legend.background = element_rect(color = "black", size = 0.5),
        legend.justification = c(1, 1),
        legend.margin = margin(6, 6, 6, 6))
leg_inside_br <-
    theme(
        legend.title = element_blank(),
        legend.position = c(1, 0),
        # legend.background = element_rect(color = "black", size = 0.5),
        legend.justification = c(1, 0),
        legend.margin = margin(6, 6, 6, 6))

##################################################################
## Data loading
##################################################################

#######################################
### healthy profiles for reference..
#######################################

# Loading the reference profiles can take around 5 minutes and doesn't need to be redone all the time.
if (reloadHealthyProfiles) {
    # load metadata.
    set.seed(1321311)
    metaOther <- read_csv('/g/scb2/zeller/SHARED/DATA/metadata/Western_non_Western_metadata_meta_analysis.csv') %>%
        rename(sampleID = sample_id) %>%
        filter(non_westernized == "Westernized")

    # load profiles. Load all and then subset based on metadata
    referenceProfiles <- map(list.files('/g/scb2/zeller/karcher/CAZY_project_v2/data/profiles/')[!str_detect(list.files('/g/scb2/zeller/karcher/CAZY_project_v2/data/profiles/'), 'plot')], function(x) readRDS(str_c('/g/scb2/zeller/karcher/CAZY_project_v2/data/profiles/', x)))
    referenceProfiles <- map(referenceProfiles, function(x) x %>%
        as.data.frame() %>%
        rownames_to_column('taxon') %>%
        pivot_longer(-taxon) %>%
        rename(sampleID = name,
            count = value))

    referenceProfiles <- do.call('rbind', referenceProfiles)

    referenceProfiles <- referenceProfiles %>%
        mutate(sampleID = str_replace(sampleID, "bgi-", "")) %>%
        mutate(sampleID = str_replace(sampleID, "[.]singles", "")) %>%
        # distinct(sampleID) %>%
        group_by(sampleID, taxon) %>%
        summarize(count = sum(count)) %>%
        left_join(read_tsv('/g/scb2/zeller/karcher/CAZY_project_v2/data/Liu_2016_sample_info') %>%
            select(run_accession,
                sample_alias) %>%
            rename(newSampleID = run_accession,
                oldSampleID = sample_alias),
        by = c('sampleID' = 'oldSampleID')) %>%
        mutate(sampleID = ifelse(is.na(newSampleID), sampleID, newSampleID)) %>%
        select(-newSampleID) %>%
        left_join(read_csv('/g/scb2/zeller/karcher/CAZY_project_v2/data/Schirmer_2016_SRA_Metadata.txt') %>%
            filter(str_detect(`Library Name`, "_pe")) %>%
            select(Run, `Sample Name`) %>%
            rename(newSampleID = Run,
                oldSampleID = `Sample Name`),
        by = c('sampleID' = 'oldSampleID')) %>%
        mutate(sampleID = ifelse(is.na(newSampleID), sampleID, newSampleID)) %>%
        select(-newSampleID)

    referenceProfiles <- referenceProfiles %>%
        pivot_wider(id_cols = taxon, names_from = sampleID, values_from = count, values_fill = 0) %>%
        pivot_longer(-taxon) %>%
        rename(sampleID = name, count = value)

    referenceProfiles <- referenceProfiles %>%
        inner_join(metaOther, by = 'sampleID') %>%
        select(-family)

    tax <- str_split_fixed(referenceProfiles$taxon, '[|]', n = 8)
    colnames(tax) <- c('kingdom', 'phylum', "class", "order", "family", "genus", 'species', 'mOTU_ID')
    referenceProfiles <- cbind(tax, referenceProfiles)

    # Calculate relAbs
    # SUBSET PROFILES LATER
    referenceProfiles <- referenceProfiles %>%
        group_by(sampleID) %>%
        mutate(relAb = count / sum(count)) %>%
        # select(-count) %>%
        # inner_join(prev_mOTUs %>%
        #              select(taxon),
        #            by = 'taxon') %>%
        filter(taxon != 'unassigned')


    referenceProfiles <- referenceProfiles %>%
        select(
            kingdom,
            phylum,
            class,
            order,
            family,
            genus,
            species,
            mOTU_ID,
            taxon,
            sampleID,
            count,
            study_name,
            relAb
        )
}

referenceProfiles %>% 
pivot_wider(id_cols = c(kingdom, phylum,  class, order, family, genus, species, mOTU_ID, taxon), names_from = sampleID, values_from = relAb) %>%
write_tsv("/g/scb/zeller/karcher/PRISMA/profiles/referenceProfilesWide.tsv")

load(url("https://github.com/AlessioMilanese/motus_taxonomy/blob/master/data/motus_taxonomy_3.0.1.Rdata?raw=true"))

species_taxID <- motus3.0_taxonomy %>%
    select(Kingdom, Phylum, Class, Order, Family, Genus, Species, profiled) %>%
    mutate(NCBI_tax_ID = str_split_fixed(Species, " ", n = 2)[, 1]) %>%
    mutate(mOTUs_ID = str_replace(profiled, ".*\\[", "")) %>%
    mutate(mOTUs_ID = str_replace(mOTUs_ID, "]", "")) %>%
    mutate(mOTUs_ID = str_replace(mOTUs_ID, "_v3_", "_v31_")) %>%
    as_tibble()

profile <- readRDS('/g/scb/zeller/fspringe/Projects/metaG_mOTUs3.1_profiles/211213_MB002A_WGS/Results/collated/res_mOTUs.rds')
colnames(profile) <- str_replace(colnames(profile), ".*_lane1", "")
profile <- profile %>%
    as.data.frame() %>%
    rownames_to_column('taxon') %>%
    pivot_longer(-taxon) %>%
    as_tibble() %>%
    rename(sampleIDGenecore = name, count = value)
meta <- read_tsv('/g/scb/zeller/karcher/PRISMA/data/WGS_metadata/MB002A_Samples_for_WGS.tsv',
    col_types = readr::cols(`Sample ID (Genecore)` = col_character())) %>%
    select(Samples, `Sample ID (Genecore)`, `Anaerob (AA)/Microaerob (MA)`) %>%
    rename(
        sampleIDGenecore = `Sample ID (Genecore)`,
        sampleID = Samples,
        oxygen = `Anaerob (AA)/Microaerob (MA)`
    )
meta %>% write_tsv('/g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/metadata/overnight_cultures_WGS_meta.tsv')
dimBefore <- dim(profile)[1]
profile <- profile %>%
    inner_join(meta)
stopifnot(dimBefore == dim(profile)[1])


profile <- profile %>%
    # as.data.frame() %>%
    anti_join(data.frame(sampleID = c("MB006S", "MB010 160321"))) %>%
    mutate(sampleID = ifelse(sampleID == "MB006A", "MB006", sampleID)) %>%
    mutate(sampleID = ifelse(sampleID == "MB010 290720", "MB010", sampleID))

sampleID_donorID_map <- read.delim(text = "    sampleID	ExpID	stoolDonor
MB001	1	A-01
MB002	2	A-02
MB003	3	A-03
MB005	4	A-04
MB010	5	A-05
MB006	6	A-06
MB007	7	A-07
MB008	8	A-08
MB009	9	A-09
MB021	21	A-10
MB015	15	P-01
MB016	16	P-02
MB017	17	P-03
MB018	18	P-04
MB013	13	T-01
MB014	14	T-02
MB019	19	T-03
MB020	20	T-04
MB011	22	T-05", sep = "\t") %>%
    as_tibble()

dimBefore <- dim(profile)[1]
profile <- profile %>%
    inner_join(sampleID_donorID_map, by = 'sampleID')
stopifnot(dimBefore == dim(profile)[1])

profile <- profile %>%
    mutate(phylum = str_split_fixed(taxon, "\\|", n = 7)[, 2]) %>%
    mutate(class = str_split_fixed(taxon, "\\|", n = 7)[, 3]) %>%
    mutate(order = str_split_fixed(taxon, "\\|", n = 7)[, 4]) %>%
    mutate(family = str_split_fixed(taxon, "\\|", n = 7)[, 5]) %>%
    mutate(genus = str_split_fixed(taxon, "\\|", n = 7)[, 6]) %>%
    mutate(species = str_split_fixed(taxon, "\\|", n = 7)[, 7]) %>%
    mutate(mOTUs_ID = str_split_fixed(taxon, "\\|", n = 8)[, 8])

profile <- profile %>%
    group_by(stoolDonor, oxygen, sampleIDGenecore) %>%
    mutate(relAb = count / sum(count)) %>%
    filter(taxon != 'unassigned')

profile %>%
ungroup() %>%
select(sampleID, stoolDonor, oxygen, phylum, class, order, family, genus, genus, species, mOTUs_ID, relAb) %>%
write_tsv('/g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/profiles/overnight_cultures_WGS_profiles_long.tsv')

library(vegan)
library(ggsignif)
library(lmerTest)
library(lme4)

tmp <- profilesTogetherLong %>%
    select(sampleID, oxygen, profileType, Genus, relAb) %>%
    group_by(sampleID, oxygen, profileType, Genus) %>%
    mutate(relAb = 10^relAb - pseudoCount) %>%
    summarize(relAb = sum(relAb)) %>%
    group_by(sampleID, oxygen, profileType) %>%
    summarize(richness = sum(relAb > pseudoCount), `Shannon` = diversity(relAb, 'shannon')) %>%
    pivot_longer(-c(sampleID, oxygen, profileType))

tmp <- tmp %>%
    rename(`Diversity index` = name) %>%
    mutate(donorType = case_when(
        str_detect(sampleID, "A-") ~ "healthy adults",
        str_detect(sampleID, "T-") ~ "transplant patients",
        str_detect(sampleID, "P-") ~ "healthy children",
    #)) %>%
    .default = "reference")) %>%
    mutate(donorType = factor(donorType, levels = c(
        "healthy adults",
        "healthy children",
        "transplant patients",
        "reference"
    )))



comparisons = list(
    c("healthy adults", "healthy children"),
    c("healthy adults", "transplant patients"),
    c("healthy children", "transplant patients"),
    c("healthy adults", "reference"),
    c("healthy children", "reference"),
    c("transplant patients", "reference")            
)

diffResults <- c()
for (compPair in comparisons) {
    a <- compPair[[1]]
    b <- compPair[[2]]

    tmp2 <- tmp %>% filter(donorType %in% c(a,b))
    lmmResultsShannon <- summary(lmer(value ~ donorType + (1|oxygen), data = tmp2 %>% filter(`Diversity index` == "Shannon")))
    lmmResultsRichness <- summary(lmer(value ~ donorType + (1|oxygen), data = tmp2 %>% filter(`Diversity index` == "richness")))

    diffResults[[length(diffResults) + 1]] <- lmmResultsShannon$coefficients[rownames(lmmResultsShannon$coefficients)!='(Intercept)', 'Pr(>|t|)']
    names(diffResults)[length(diffResults)] <- str_c(str_c(compPair[[1]], compPair[[2]], sep = "___"), "Shannon", sep = "___")

    diffResults[[length(diffResults) + 1]] <- lmmResultsRichness$coefficients[rownames(lmmResultsShannon$coefficients)!='(Intercept)', 'Pr(>|t|)']
    names(diffResults)[length(diffResults)] <- str_c(str_c(compPair[[1]], compPair[[2]], sep = "___"), "richness", sep = "___")    
}

bla <- tmp %>%
group_by(`Diversity index`) %>%
summarize(maxVal = max(value))

maxVals <- bla$maxVal
names(maxVals)<- bla$`Diversity index`

annotations <- data.frame(diffResults, check.names = FALSE) %>%
t() %>%
data.frame(check.names = FALSE) %>%
rename(values = '.') %>%
rownames_to_column('raw') %>%
as_tibble() %>%
mutate(comparisons = map(raw, \(x){
    tmp <- str_split_fixed(x, "___", n = 3)[, 1:2]
    return(c(tmp[[1]], tmp[[2]]))
})) %>%
mutate(`Diversity index` = map_chr(raw, \(x) str_split_fixed(x, "___", n = 3)[, 3])) %>%
mutate(comparisonsLong = map_chr(comparisons, \(x) str_c(x[[1]], x[[2]], sep = "___"))) %>%
mutate(label = round(values, 4))  %>%
mutate(start = map_chr(comparisons, \(x) x[[1]])) %>%
mutate(end = map_chr(comparisons, \(x) x[[2]])) %>%
group_by(`Diversity index`) %>%
nest() %>%
mutate(data = map2(data, `Diversity index`, \(x, di) {
    desiredOrder <- map_chr(comparisons, \(x) str_c(x[[1]], x[[2]], sep = "___"))
    x <- x[match(desiredOrder, x$comparisonsLong),]
    #maxVal  <- max(x$values)
    maxVal <- maxVals[[di]]
    x$y <- map_dbl((1:dim(x)[1])/dim(x)[1], \(x) {
        return(maxVal + maxVal * x * 0.75)
    })
    return(x)
})) %>%
unnest()


(
    #mutate(SampleType = ifelse(is.na(donorType), "reference", "Stool donor")) %>%    
    ggplot(data = tmp) +
    facet_wrap( ~ `Diversity index`, scales = "free") +
    #    theme_publication() +
    # theme(
    #     axis.text.x = element_blank(),
    #     axis.ticks.x = element_blank(),
    #     axis.title.x = element_blank()) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_boxplot(aes(x = donorType, y = value)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    # geom_signif(
    #     comparisons = list(
    #         c("healthy adults", "healthy children"),
    #         c("healthy adults", "transplant patients"),
    #         c("healthy children", "transplant patients"),
    #         c("healthy adults", "reference"),
    #         c("healthy children", "reference"),
    #         c("transplant patients", "reference")            
    #     ),
    #     step_increase = 0.15,
    #     map_signif_level = FALSE,
    #     size = 0.2,
    #     textsize = 2.25
    # ) +
    geom_signif(data = annotations %>% select(start, end, label, y, `Diversity index`), aes(
        xmin = start, 
        xmax = end, 
        annotations = label,
        y_position = y), manual = TRUE, textsize = 2.25, size = 0.2) +
    ylab("Diversity measure") + 
    xlab("Sample type")) %>%
    #ggsave(plot = ., filename = "/g/scb/zeller/karcher/PRISMA/plots/AAAKVT7HR/diversity_measures_naive_wilcox_test.pdf", width = 4, height = 3.5)
    ggsave(plot = ., filename = "/g/scb/zeller/karcher/PRISMA/plots/AAAKVT7HR/diversity_measures_lmm.pdf", width = 4, height = 3.5)

# (profile %>%
#     mutate(genus = str_split_fixed(taxon, "[|]", 7)[, 6]) %>%
#     group_by(stoolDonor, oxygen, genus, relAb) %>%
#     group_by(stoolDonor, oxygen, genus) %>%
#     summarize(relAb = sum(relAb)) %>%
#     group_by(stoolDonor, oxygen) %>%
#     summarize(richness = sum(relAb > pseudoCount), `Shannon` = diversity(relAb, 'shannon')) %>%
#     pivot_longer(-c(stoolDonor, oxygen)) %>%
#     rename(`Diversity index` = name) %>%
#     mutate(donorType = case_when(
#         str_detect(stoolDonor, "A-") ~ "healthy adults",
#         str_detect(stoolDonor, "T-") ~ "transplant patients",
#         str_detect(stoolDonor, "P-") ~ "healthy children",
#     .default = NA)) %>%
#     ggplot(aes(x = donorType, y = value)) +
#     facet_wrap(. ~ `Diversity index`, scales = "free") +
#     #    theme_publication() +
#     # theme(
#     #     axis.text.x = element_blank(),
#     #     axis.ticks.x = element_blank(),
#     #     axis.title.x = element_blank()) +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     geom_boxplot() +
#     scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15))) +
#     geom_signif(
#         comparisons = list(
#             c("healthy adults", "healthy children"),
#             c("healthy adults", "transplant patients"),
#             c("healthy children", "transplant patients")
#         ),
#         step_increase = 0.3,
#         map_signif_level = FALSE,
#         size = 0.2
#     ) +
#     ylab("Diversity Measure")) %>%
#     ggsave(plot = ., filename = "/g/scb/zeller/karcher/PRISMA/plots/AAAKVT7HR/diversity_measures.pdf", width = 4, height = 3.5)

drugMetabAA <- read_xlsx('/g/scb/zeller/karcher/PRISMA/data/WGS_metadata/Tables_STM.xlsx', sheet = 7, skip = 1)
drugMetabMA <- read_xlsx('/g/scb/zeller/karcher/PRISMA/data/WGS_metadata/Tables_STM.xlsx', sheet = 8, skip = 1)
drugMetab <- rbind(
    drugMetabAA %>%
        mutate(oxygen = "AA"),
    drugMetabMA %>%
        mutate(oxygen = "MA")) %>%
    rename(stoolDonor = Stool_donor)
testedBacteria <- read_xlsx('/g/scb/zeller/karcher/PRISMA/data/WGS_metadata/Tables_STM.xlsx', sheet = 4, skip = 2)
testedBacteria <- testedBacteria %>%
    #rename(Name = `Human gut bacteria tested for drug-metabolizing activity`) %>%
    filter(!is.na(`Growth condition`)) %>%
    select(Name, `Growth condition`, `NCBI tax ID`) %>%
    rename(
        bacterium = Name,
        growth_condition = `Growth condition`,
        NCBI_strain_tax_ID = `NCBI tax ID`,
    ) %>%
    mutate(MA = str_detect(str_to_lower(growth_condition), "micro")) %>%
    mutate(AA = str_detect(str_to_lower(growth_condition), "anaerobic", )) %>%
    # select(-growth_condition) %>%
    # This table was created using the script map_marals_taxids_to_species_taxids.py
    left_join(read_csv('/g/scb/zeller/karcher/PRISMA/data/WGS_metadata/straintaxid_speciestaxid.csv') %>%
        select(-`...1`), by = c('NCBI_strain_tax_ID' = "originaltaxID")) %>%
    mutate(speciesTaxID = as.character(speciesTaxID))

# Take from CAZy project...
##################################################################
# mOTUs3 depths
##################################################################

p <- ggplot(
    profile %>%
        group_by(sampleIDGenecore, stoolDonor, oxygen) %>%
        summarize(totalmOTUDepth = sum(count))
) + geom_boxplot(aes(x = oxygen, y = totalmOTUDepth)) +
    theme_classic()

ggsave(plot = p, filename = "/g/scb/zeller/karcher/PRISMA/plots/AAAKVT7HR/depth_boxplot.pdf", width = 4, height = 4)

##################################################################
# Are our samples representative of a healthy human gut microbiome?
##################################################################

pwdistances <- vegan::vegdist(profile %>%
    mutate(log10relAb = log10(relAb + pseudoCount)) %>%
    select(taxon, sampleIDGenecore, log10relAb) %>%
    pivot_wider(id_cols = taxon, names_from = sampleIDGenecore, values_from = log10relAb) %>%
    as.data.frame() %>%
    column_to_rownames('taxon') %>%
    t(), method = 'euclidean')

pcoa <- cmdscale(pwdistances, k = 2)
pcoa <- pcoa %>%
    as.data.frame()
colnames(pcoa) <- c("PCo 1", "PCo 2")
pcoa <- pcoa %>%
    rownames_to_column('sampleIDGenecore') %>%
    left_join(
        profile %>%
            ungroup() %>%
            select(sampleIDGenecore, oxygen, sampleID, stoolDonor) %>%
            distinct(),
        by = "sampleIDGenecore"
    ) %>%
    mutate(sampleType = case_when(
        str_detect(stoolDonor, "A-") ~ "healthy adults",
        str_detect(stoolDonor, "T-") ~ "transplant patients",
        str_detect(stoolDonor, "P-") ~ "healthy children"
    ))

meta <- data.frame(sampleID = rownames(as.matrix(pwdistances))) %>%
    left_join(pcoa %>% select(sampleIDGenecore, oxygen, sampleID, stoolDonor), by = c('sampleID' = 'sampleIDGenecore')) %>%
    mutate(sampleType = case_when(
        str_detect(stoolDonor, "A-") ~ "healthy adults",
        str_detect(stoolDonor, "T-") ~ "transplant patients",
        str_detect(stoolDonor, "P-") ~ "healthy children"
    )) %>%
    column_to_rownames('sampleID')

#meta %>% write_tsv('/g/scb2/zeller/karcher/tmp/meta.tsv')
#pwdistancesBig %>% as.matrix() %>% as.data.frame() %>% rownames_to_column('sampleID') %>% write_tsv('/g/scb2/zeller/karcher/tmp/b.tsv')

stopifnot(all(rownames(meta) == rownames(as.matrix(pwdistances))))
# This runs for a little while because we have large N
permanova <- adonis2(pwdistances ~ sampleType, data = meta)
permanova <- adonis2(pwdistances ~ oxygen, data = meta)


p <- ggplot(
    pcoa
) +
    geom_line(aes(x = `PCo 1`, y = `PCo 2`, group = sampleID)) +
    geom_point(aes(x = `PCo 1`, y = `PCo 2`, shape = oxygen, color = sampleType)) +
    geom_text_repel(data = pcoa %>% filter(oxygen == "MA"), aes(x = `PCo 1`, y = `PCo 2`, label = stoolDonor), size = 2) +
    theme_classic()

ggsave(plot = p, filename = '/g/scb/zeller/karcher/PRISMA/plots/AAAKVT7HR/ordination_only_AAAKVT7HR.pdf', width = 6, height = 4)



profilesTogether <- rbind(
    profile %>%
        ungroup() %>%
        select(stoolDonor, oxygen, mOTUs_ID, relAb) %>%
        rename(sampleID = stoolDonor) %>%
        mutate(profileType = "overnightCulture"),
    referenceProfiles %>%
        ungroup() %>%
        select(sampleID, mOTU_ID, relAb) %>%
        rename(mOTUs_ID = mOTU_ID) %>%
        mutate(profileType = 'WGSProfile', oxygen = "NA")
) %>%
    mutate(log10relAb = log10(relAb + pseudoCount)) %>%
    pivot_wider(id_cols = c(sampleID, oxygen, profileType), names_from = mOTUs_ID, values_from = log10relAb, values_fill = log10(pseudoCount))

# This runs for at least a few minutes, if not longer...
# CAREFUL! doing low-prev filtering for computational time savings...
# tmp2 <- profilesTogether %>% select(c(sampleID, oxygen, profileType))
# tmp3 <- profilesTogether %>% select(-c(sampleID, oxygen, profileType))
# profilesTogether <- cbind(tmp2, tmp3[, apply(tmp3, 2, function(x) mean(x > -5) > 0.1)])

# Runs for a minute or two.
if (reloadHealthyProfiles) {
    pwdistancesBig <- vegan::vegdist(profilesTogether %>%
        mutate(profilesTogether = str_c(sampleID, as.character(oxygen), profileType, sep = "___")) %>%
        select(-c(sampleID, oxygen, profileType)) %>%
        as.data.frame() %>%
        column_to_rownames('profilesTogether'), method = 'euclidean')
}

pcoaBig <- cmdscale(pwdistancesBig, k = 2)
pcoaBig <- pcoaBig %>%
    as.data.frame()
colnames(pcoaBig) <- c("PCo 1", "PCo 2")
pcoaBig <- pcoaBig %>%
    rownames_to_column('tmpSampleID') %>%
    left_join(
        profilesTogether %>%
            ungroup() %>%
            mutate(tmpSampleID = str_c(sampleID, as.character(oxygen), profileType, sep = "___")) %>%
            select(sampleID, oxygen, profileType, tmpSampleID) %>%
            distinct(),
        by = "tmpSampleID"
    )

pcoaBig <- pcoaBig %>%
    mutate(sampleType = case_when(
        str_detect(sampleID, "A-") ~ "healthy adults",
        str_detect(sampleID, "T-") ~ "transplant patients",
        str_detect(sampleID, "P-") ~ "healthy children"
    ))


meta <- data.frame(raw = rownames(as.matrix(pwdistancesBig))) %>%
    mutate(profileType = str_split_fixed(raw, "___", n = 3)[, 3]) %>%
    column_to_rownames('raw')

#meta %>% write_tsv('/g/scb2/zeller/karcher/tmp/meta.tsv')
#pwdistancesBig %>% as.matrix() %>% as.data.frame() %>% rownames_to_column('sampleID') %>% write_tsv('/g/scb2/zeller/karcher/tmp/b.tsv')

stopifnot(all(rownames(meta) == rownames(as.matrix(pwdistancesBig))))
# This runs for a little while because we have large N
permanova <- adonis2(pwdistancesBig ~ profileType, data = meta)

# Create a named vector with the colors
#ordination_color_vector <- c("healthy adults" = "#ACB78E", "healthy children" = "#F88379", "transplant patients" = "#B0C4DE")
# Update the color vector with the new blue
ordination_color_vector <- c("healthy adults" = "#19984A", "healthy children" = "#317EC2", "transplant patients" = "#E7872B")

p1Ordination <- ggplot() +
    geom_point(data = pcoaBig %>%
        filter(profileType == "WGSProfile"), aes(x = `PCo 1`, y = `PCo 2`), color = 'grey', alpha = 0.25) +
    geom_line(data = pcoaBig %>%
        filter(profileType != "WGSProfile"), aes(x = `PCo 1`, y = `PCo 2`, group = sampleID), color = 'black', alpha = 0.5) +
    geom_point(data = pcoaBig %>%
        filter(profileType != "WGSProfile"), aes(x = `PCo 1`, y = `PCo 2`, color = sampleType, shape = oxygen), alpha = 1) +
    theme_classic() +
    # scale_alpha(range  = c(0.25, 1)) +
    guides(alpha = 'none') +
    scale_color_manual(values = ordination_color_vector)


get_family_level_barplot <- function(pObj, dataB, taxLevel = 'Family', levelsToShow = NULL, pc = pseudoCount) {
    dataB <- dataB %>%
        mutate(relAb = (10^relAb) - pc) %>%
        mutate(taxa = .data[[taxLevel]]) %>%
        mutate(taxa = as.character(taxa)) %>%
        mutate(taxa = ifelse(taxa %in% levelsToShow, taxa, "other")) %>%
        # mutate(taxa = factor(taxa, levels = c(levelsToShow, "other", 'unclassified'))) %>%
        group_by(taxa, sampleID, oxygen) %>%
        summarize(relAb = sum(relAb))
    dataC <- dataB %>%
        ungroup() %>%
        group_by(sampleID, oxygen) %>%
        summarize(relAb = 1 - sum(relAb)) %>%
        mutate(taxa = "unclassified")
    dataB <- rbind(dataB, dataC) %>%
        mutate(taxa = factor(as.character(taxa), levels = c(levelsToShow, "other", 'unclassified')))

    dataB <- dataB %>%
        mutate(sampleID = str_c(sampleID, oxygen, sep = "__")) %>%
        mutate(oxygen = str_split_fixed(sampleID, "__", n = 2)[, 2]) %>%
        mutate(sampleID = str_split_fixed(sampleID, "__", n = 2)[, 1]) %>%
        mutate(oxygen = case_when(
            oxygen == "MA" ~ "Microaerobic",
            oxygen == "AA" ~ "Anaerobic",
        ))

    pObj <- pObj +
        geom_bar(data = dataB,
            aes(x = sampleID, y = relAb, fill = taxa), position = 'stack', stat = 'identity') +
        theme_classic() +
        facet_grid(. ~ oxygen)
    return(pObj)
}


profilesTogetherLong <- profilesTogether %>%
    pivot_longer(-c(sampleID, oxygen, profileType)) %>%
    rename(mOTUs_ID = name, relAb = value) %>%
    left_join(species_taxID %>%
        mutate(mOTUs_ID = str_replace(mOTUs_ID, "_v3_", "_v31_")) %>%
        mutate(Genus = str_split_fixed(Genus, " ", n = 7)[, 2]) %>%
        mutate(Family = str_split_fixed(Family, " ", n = 7)[, 3]) %>%
        mutate(Order = str_split_fixed(Order, " ", n = 7)[, 4]) %>%
        mutate(Class = str_split_fixed(Class, " ", n = 7)[, 5]) %>%
        mutate(Phylum = str_split_fixed(Phylum, " ", n = 7)[, 6]), by = 'mOTUs_ID') %>%
    select(-Species)

for (taxL in c("Phylum", "Class", "Order", "Family", "Genus")) {
    numTaxa <- 10
    levelsToShow <- profilesTogetherLong %>%
        mutate(relAb = 10^relAb - pseudoCount) %>%
        group_by(sampleID, .data[[taxL]]) %>%
        summarize(relAb = sum(relAb)) %>%
        group_by(.data[[taxL]]) %>%
        summarize(m = mean(relAb)) %>%
        filter(!str_detect(.data[[taxL]], "incertae")) %>%
        arrange(desc(m)) %>%
        head(numTaxa) %>%
        pull(.data[[taxL]])
    getPalette <- colorRampPalette(brewer.pal(9, "Paired"))
    set.seed(1231321)
    colors <- sample(getPalette(numTaxa))
    colors <- c(colors, "#808080", "#D3D3D3")

    p2 <- get_family_level_barplot(
        ggplot(),
        profilesTogetherLong %>%
            group_by(sampleID, oxygen, profileType) %>%
            nest() %>%
            filter(profileType == "overnightCulture") %>%
            ungroup() %>%
            # if you want a bit more resolution
            # slice_sample(n = 10) %>%
            unnest(),
        taxL,
        levelsToShow = levelsToShow) +
        scale_fill_manual(values = colors) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5)) +
        ggtitle("Overnight Cultures") +
        guides(fill = guide_legend(ncol = 2)) +
        theme_publication() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    # set.seed(12313312)
    # p3 <- get_family_level_barplot(
    #     ggplot(),
    #     profilesTogetherLong %>%
    #         group_by(sampleID, oxygen, profileType) %>%
    #         nest() %>%
    #         filter(profileType == "WGSProfile") %>%
    #         ungroup() %>%
    #         # Same number of bars as for overnightCulture..
    #         slice_sample(n = 38) %>%
    #         unnest(),
    #     'Family',
    #     levelsToShow = levelsToShow) +
    #                 scale_fill_manual(values = colors) +
    #                     theme(
    #                         axis.text.x = element_blank(),
    #                         axis.ticks.x = element_blank(),
    #                         axis.title.x = element_blank()
    #                     ) +
    #                     ggtitle("Stool WGS profiles") +
    #                     guides(fill = guide_legend(ncol = 1))

    ggsave(plot = p2,
        filename = str_c('/g/scb/zeller/karcher/PRISMA/plots/AAAKVT7HR/AAAKVT7HR_', taxL, '_level_barplot.pdf'), width = 6.5, height = 2)

}

#####################################################################################
# Are tested strains represntative of the tested unrelated communities?
#####################################################################################

testedSpeciesRealMetaG <- referenceProfiles %>%
    rename(mOTUs_ID = mOTU_ID) %>%
    #  We're losing a single 11 mOTUs here, presumably because they've been added in mOTUs3.
    # They are all basically always 0 though, so no reason to worry. Retrace by replacing inner_join with anti_join.
    inner_join(species_taxID %>% select(NCBI_tax_ID, mOTUs_ID), by = 'mOTUs_ID') %>%
    inner_join(testedBacteria %>% select(-species), by = c('NCBI_tax_ID' = 'speciesTaxID'))

print(str_c("Out of the ", dim(testedBacteria)[1], " in cohort metaGs, we can map ", length(unique(testedSpeciesRealMetaG$NCBI_tax_ID)), " species directly via NCBI species taxIDs."))
print("Let's investigate the 6 that are missing")
notDirectlyMappable <- anti_join(testedBacteria %>% select(bacterium, speciesTaxID), testedSpeciesRealMetaG %>% ungroup() %>% select(NCBI_tax_ID) %>% distinct(), by = c("speciesTaxID" = "NCBI_tax_ID"))
print("Some mOTUs can not be mapped because they are not having an NCBI tax id, in turn because they contain more than one species")
print("Let us add them here via their subspecies label")
for (rowIndex in 1:dim(notDirectlyMappable)[1]) {
    # print(notDirectlyMappable[rowIndex, ]$bacterium)
    print(str_c("For ", notDirectlyMappable[rowIndex, ]$bacterium, ", we can map ", str_c(rep(" ", 47 - str_length(str_c(str_c("For ", notDirectlyMappable[rowIndex, ]$bacterium, ", we can map ")))), collapse = ""),
        referenceProfiles %>%
            rename(mOTUs_ID = mOTU_ID) %>%
            #  We're losing a single 11 mOTUs here, presumably because they've been added in mOTUs3.
            # They are all basically always 0 though, so no reason to worry. Retrace by replacing inner_join with anti_join.
            inner_join(species_taxID, by = 'mOTUs_ID') %>%
            ungroup() %>%
            select(taxon) %>%
            distinct() %>%
            filter(str_detect(taxon, str_split(notDirectlyMappable[rowIndex, ]$bacterium, " ")[[1]][2])) %>%
            pull(taxon) %>%
            table() %>%
            length(), " motus via the subspecies label/identifier"))
    testedSpeciesRealMetaG <- rbind(testedSpeciesRealMetaG, referenceProfiles %>%
        rename(mOTUs_ID = mOTU_ID) %>%
        #  We're losing a single 11 mOTUs here, presumably because they've been added in mOTUs3.
        # They are all basically always 0 though, so no reason to worry. Retrace by replacing inner_join with anti_join.
        inner_join(species_taxID, by = 'mOTUs_ID') %>%
        filter(str_detect(taxon, str_split(notDirectlyMappable[rowIndex, ]$bacterium, " ")[[1]][2])) %>%
        mutate(bacterium = notDirectlyMappable[rowIndex, ]$bacterium))
}
testedSpeciesRealMetaG <- testedSpeciesRealMetaG %>% distinct(mOTUs_ID, .keep_all = T)
print(str_c("In the end, I have ", length(unique(testedSpeciesRealMetaG$NCBI_tax_ID)), " distinct species and ", length(unique(testedSpeciesRealMetaG$mOTUs_ID)), " distinct mOTUs in the tested communities"))
testedSpeciesRealMetaG %>% ungroup() %>% select(mOTUs_ID, `Original name`) %>% distinct() %>% mutate(Species  = map_chr(`Original name`, \(x) return(str_c(str_c(str_split(x, " ")[[1]][1:2], collapse =  " "))))) %>%write_tsv('/g/scb/zeller/karcher/PRISMA/data/WGS_metadata/species_mOTUs_link.tsv')


#############################################################
# Are tested strains represntative of the tested commmunities?
#############################################################

testedSpecies <- profile %>%
    #  We're losing a single 11 mOTUs here, presumably because they've been added in mOTUs3.
    # They are all basically always 0 though, so no reason to worry. Retrace by replacing inner_join with anti_join.
    inner_join(species_taxID %>% select(NCBI_tax_ID, mOTUs_ID), by = 'mOTUs_ID') %>%
    inner_join(testedBacteria %>% select(-species), by = c('NCBI_tax_ID' = 'speciesTaxID'))

print(str_c("Out of the ", dim(testedBacteria)[1], " Maral tested, we can map ", length(unique(testedSpecies$NCBI_tax_ID)), " species directly via NCBI species taxIDs."))
print("Let's investigate the 10 that are missing")
notDirectlyMappable <- anti_join(testedBacteria %>% select(bacterium, speciesTaxID), testedSpecies %>% ungroup() %>% select(NCBI_tax_ID) %>% distinct(), by = c("speciesTaxID" = "NCBI_tax_ID"))
print("Some mOTUs can not be mapped because they are not having an NCBI tax id, in turn because they contain more than one species")
print("Let us add them here via their subspecies label")
for (rowIndex in 1:dim(notDirectlyMappable)[1]) {
    # print(notDirectlyMappable[rowIndex, ]$bacterium)
    print(str_c("For ", notDirectlyMappable[rowIndex, ]$bacterium, ", we can map ", str_c(rep(" ", 47 - str_length(str_c(str_c("For ", notDirectlyMappable[rowIndex, ]$bacterium, ", we can map ")))), collapse = ""),
        profile %>%
            #  We're losing a single 11 mOTUs here, presumably because they've been added in mOTUs3.
            # They are all basically always 0 though, so no reason to worry. Retrace by replacing inner_join with anti_join.
            inner_join(species_taxID, by = 'mOTUs_ID') %>%
            filter(str_detect(taxon, str_split(notDirectlyMappable[rowIndex, ]$bacterium, " ")[[1]][2])) %>%
            pull(taxon) %>%
            table() %>%
            length(), " motus via the subspecies label/identifier"))
    testedSpecies <- rbind(testedSpecies, profile %>%
        #  We're losing a single 11 mOTUs here, presumably because they've been added in mOTUs3.
        # They are all basically always 0 though, so no reason to worry. Retrace by replacing inner_join with anti_join.
        inner_join(species_taxID, by = 'mOTUs_ID') %>%
        filter(str_detect(taxon, str_split(notDirectlyMappable[rowIndex, ]$bacterium, " ")[[1]][2])))
}
testedSpecies <- testedSpecies %>% distinct()
print(str_c("In the end, I have ", length(unique(testedSpecies$NCBI_tax_ID)), " distinct species and ", length(unique(testedSpecies$mOTUs_ID)), " distinct mOTUs in the tested communities"))


# p1 <- ggplot(
#     testedSpecies %>% group_by(sampleID, oxygen) %>% summarize(cumRelAb = sum(relAb))
#     ) +
#         geom_line(aes(x = oxygen, y = cumRelAb, group = sampleID), alpha = 0.35) +
#         geom_boxplot(aes(x = oxygen, y = cumRelAb)) +
#             theme_classic() +
#             ylab("Cumulative relative Abundance\n of tested strains in overnight cultures")

# p1.1 <- ggplot(
#     testedSpeciesRealMetaG %>% group_by(sampleID) %>% summarize(cumRelAb = sum(relAb))
#     ) +
#         # geom_line(aes(x = 1, y = cumRelAb, group = sampleID), alpha = 0.35) +
#         geom_boxplot(aes(x = 1, y = cumRelAb)) +
#             theme_classic() +
#             ylab("Cumulative relative Abundance\n of tested strains in gut microbial communities") +
#             theme(
#                 axis.text.x = element_blank(),
#                 axis.title.x = element_blank(),
#                 axis.ticks.x = element_blank()
#             )

p1 <- testedSpecies %>%
    group_by(sampleID, oxygen) %>%
    summarize(cumRelAb = sum(relAb)) %>%
    mutate(group = ifelse(oxygen == "MA", "ON cultures (MA)", "ON cultures (AA)")) %>%
    rbind(testedSpeciesRealMetaG %>%
        group_by(sampleID) %>%
        summarize(cumRelAb = sum(relAb)) %>%
        mutate(group = "gut microbial communities")) %>%
    ggplot(aes(x = group, y = cumRelAb)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Cumulative relative Abundance\n of tested strains in gut microbial communities") +
    geom_boxplot()

p2 <- profile %>%
    mutate(genus = str_replace(genus, "g__", "")) %>%
    inner_join(testedBacteria %>% mutate(genus = str_split_fixed(bacterium, " ", n = 2)[, 1]) %>% select(genus) %>% distinct()) %>%
    group_by(sampleID, oxygen) %>%
    summarize(cumRelAb = sum(relAb)) %>%
    mutate(group = ifelse(oxygen == "MA", "ON cultures (MA)", "ON cultures (AA)")) %>%
    rbind(referenceProfiles %>%
        mutate(genus = str_replace(genus, "g__", "")) %>%
        inner_join(testedBacteria %>% mutate(genus = str_split_fixed(bacterium, " ", n = 2)[, 1]) %>% select(genus) %>% distinct()) %>%
        group_by(sampleID) %>%
        summarize(cumRelAb = sum(relAb)) %>%
        mutate(group = "gut microbial communities")) %>%
    ggplot(aes(x = group, y = cumRelAb)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Cumulative relative Abundance\n of tested genera in gut microbial communities") +
    geom_boxplot()

ggsave(plot = p1 + p2, filename = "/g/scb/zeller/karcher/PRISMA/plots/AAAKVT7HR/cumulative_relAb.pdf", width = 6, height = 4.5)




# genusDiversity <- inner_join(
#     referenceProfiles %>%
#         group_by(sampleID, genus) %>%
#         summarize(relAb = sum(relAb)) %>%
#         filter(relAb > pseudoCount) %>%
#         group_by(sampleID) %>%
#         summarize(genusDiversityTotal = n()),
#     testedSpeciesRealMetaG %>%
#         group_by(sampleID, genus) %>%
#         summarize(relAb = sum(relAb)) %>%
#         filter(relAb > pseudoCount) %>%
#         group_by(sampleID) %>%
#         summarize(genusDiversityRepresented = n())
# ) %>%
#     mutate(genusDiversityRepresentedFractional = genusDiversityRepresented / genusDiversityTotal)
# p1.1.1 <- ggplot(
#     genusDiversity
#     ) +
#             geom_boxplot(aes(x = 1, y = genusDiversityRepresentedFractional)) +
#             theme_classic() +
#             ylab("tested genera/\nall genera in gut microbial communities") +
#         theme(
#             axis.text.x = element_blank(),
#             axis.title.x = element_blank(),
#             axis.ticks.x = element_blank()
#         )


# p2 <- ggplot(
#      testedSpecies %>% group_by(sampleID, oxygen) %>% mutate(detected = relAb > pseudoCount) %>% summarize(strainsDetected = sum(detected))
#     ) +
#         geom_line(aes(x = oxygen, y = strainsDetected, group = sampleID), alpha = 0.35) +
#         geom_boxplot(aes(x = oxygen, y = strainsDetected)) +
#             theme_classic() +
#             ylab("Number of mOTUs detected in tested communities\n(out of a total of 50; relAb > -1E-5)") +
#             ylim(c(0, 50))

# p2.1 <- ggplot(
#      testedSpeciesRealMetaG %>% group_by(sampleID) %>% mutate(detected = relAb > pseudoCount) %>% summarize(strainsDetected = sum(detected))
#     ) +
#         #geom_line(aes(x = oxygen, y = strainsDetected, group = sampleID), alpha = 0.35) +
#         geom_boxplot(aes(x = 1, y = strainsDetected)) +
#             theme_classic() +
#             ylab("Number of mOTUs detected in gut microbial communities\n(out of a total of 57; relAb > -1E-5)") +        theme(
#             axis.text.x = element_blank(),
#             axis.title.x = element_blank(),
#             axis.ticks.x = element_blank()
#         )

# genusDiversity <- inner_join(
#     profile %>%
#         mutate(genus = str_split_fixed(taxon, "[|]", n = 7)[, 6]) %>%
#         group_by(sampleID, genus, oxygen) %>%
#         summarize(relAb = sum(relAb)) %>%
#         filter(relAb > pseudoCount) %>%
#         group_by(sampleID, oxygen) %>%
#         summarize(genusDiversityTotal = n()),
#     testedSpecies %>%
#         mutate(genus = str_split_fixed(taxon, "[|]", n = 7)[, 6]) %>%
#         group_by(sampleID, genus) %>%
#         summarize(relAb = sum(relAb)) %>%
#         filter(relAb > pseudoCount) %>%
#         group_by(sampleID) %>%
#         summarize(genusDiversityRepresented = n())
# ) %>%
#     mutate(genusDiversityRepresentedFractional = genusDiversityRepresented / genusDiversityTotal)
# p2.1.1 <- ggplot(
#     genusDiversity
#     ) +
#             geom_boxplot(aes(x = oxygen, y = genusDiversityRepresentedFractional)) +
#             theme_classic() +
#             ylab("tested genera/\nall genera in tested communities")

tmp <- testedSpecies %>%
    ungroup() %>%
    select(taxon, sampleID, oxygen, relAb, growth_condition) %>%
    group_by(oxygen, taxon, growth_condition) %>%
    summarize(medianRelAb = median(relAb), prevalence = mean(relAb > pseudoCount)) %>%
    pivot_wider(id_cols = c(taxon, growth_condition), names_from = oxygen, values_from = c(medianRelAb, prevalence)) %>%
    mutate(MA_relAb_higher_AA = medianRelAb_MA > medianRelAb_AA) %>%
    mutate(prevalenceDiff = prevalence_MA - prevalence_AA) %>%
    mutate(absoluteRelAbDiff = medianRelAb_MA - medianRelAb_AA) %>%
    mutate(fc = medianRelAb_MA / medianRelAb_AA) %>%
    arrange(desc(fc)) %>%
    left_join(profile %>%
        group_by(taxon) %>%
        summarize(taxonMedianAbundance = median(relAb))) %>%
    mutate(log10TMA = log10(taxonMedianAbundance + pseudoCount)) %>%
    arrange(desc(prevalenceDiff)) %>%
    mutate(species = str_split_fixed(taxon, "[|]", n = 8)[, 7]) %>%
    mutate(species = map_chr(species, function(x) str_c(str_split(x, " ")[[1]][1:2], sep = " ", collapse = " "))) %>%
    mutate(species = str_replace(species, "\\[", "")) %>%
    mutate(species = str_replace(species, "\\]", ""))
tmp$taxon <- factor(tmp$taxon, levels = tmp$taxon)
p3 <- ggplot(
    tmp %>%
        mutate(`prevalence/abundance` = ifelse(prevalenceDiff > 0, "Higher in MA", "Higher in AA"))
) + geom_hline(yintercept = 0, alpha = 0.5) +
    geom_point(aes(x = taxon, y = prevalenceDiff, color = `prevalence/abundance`)) +
    theme_classic() +
    scale_x_discrete(labels = tmp$species) +
    coord_flip() +
    ylab("Prevalence diff.") +
    theme(legend.position = 'bottom', legend.direction = 'vertical') +
    guides(color = guide_legend(nrow = 2))

p4 <- ggplot(
    tmp %>%
        mutate(abundance = ifelse(absoluteRelAbDiff > 0, "Higher in MA", "Higher in AA"))
) + geom_hline(yintercept = 0, alpha = 0.5) +
    geom_point(aes(x = taxon, y = log2(fc), color = abundance)) +
    theme_classic() +
    scale_x_discrete(labels = tmp$species) +
    coord_flip() +
    theme(
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()
    ) +
    ylab("log2(fc)") +
    guides(color = 'none')

p5 <- ggplot() +
    geom_point(data = tmp %>%
        mutate(abundance = ifelse(absoluteRelAbDiff > 0, "Higher in MA", "Higher in AA")) %>%
        select(taxon, medianRelAb_MA, medianRelAb_AA) %>%
        pivot_longer(-taxon) %>%
        rename(relAb = value) %>%
        filter(name == "medianRelAb_AA"), aes(x = taxon, y = log10(relAb + 1E-5)), color = '#F8766D') +
    geom_point(data = tmp %>%
        mutate(abundance = ifelse(absoluteRelAbDiff > 0, "Higher in MA", "Higher in AA")) %>%
        select(taxon, medianRelAb_MA, medianRelAb_AA) %>%
        pivot_longer(-taxon) %>%
        rename(relAb = value) %>%
        filter(name == "medianRelAb_MA"), aes(x = taxon, y = log10(relAb + 1E-5)), color = '#00BFC4') +
    coord_flip() +
    theme_classic() +
    theme(
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()
    )


p6 <- ggplot(
    tmp %>%
        mutate(abundance = ifelse(absoluteRelAbDiff > 0, "Higher in MA", "Higher in AA"))
) + geom_tile(aes(x = taxon, y = 1, fill = growth_condition)) +
    coord_flip() + theme(
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()
    )



##################################################################
### Do differential abundance testing between MA and AA samples ###
##################################################################

calcGFC <- function(
    x,
    conditionCol = 'non_westernized',
    levs = c("Non_Westernized", "Westernized"),
    probs.fc = seq(.05, .95, .05)
    ) {
    x.pos <- x$relAb[x[[conditionCol]] == levs[1]]
    x.neg <- x$relAb[x[[conditionCol]] == levs[2]]
    q.p <- quantile(x.pos, probs = probs.fc)
    q.n <- quantile(x.neg, probs = probs.fc)
    return(sum(q.p - q.n) / length(q.p))
}


wilcox <- profile %>%
    # inner_join(profile %>% group_by(species) %>% summarize(prev = mean(relAb > pseudoCount)) %>% filter(prev > 0.1) %>% select(species)) %>%
    group_by(species) %>%
    mutate(relAb = log10(relAb + pseudoCount)) %>%
    filter(any(max(relAb) > -4) & mean(relAb > -5) > 0.1) %>%
    nest() %>%
    mutate(wilcox = map(data, \(x){
        wilcox.test(x = x$relAb[x$oxygen == "MA"], y = x$relAb[x$oxygen == "AA"])
    })) %>%
    mutate(wilcox.p = map_dbl(wilcox, \(x) return(x$p.value))) %>%
    arrange(wilcox.p) %>%
    mutate(gFC = map_dbl(data, \(x) {
        return(calcGFC(x, "oxygen", c("AA", "MA")))
    }))
wilcox$wilcox.p.adj <- p.adjust(wilcox$wilcox.p, method = 'BH')
wilcox <- wilcox %>%
    left_join(
        testedSpecies %>%
            ungroup() %>%
            select(species) %>%
            distinct() %>%
            mutate(tested = TRUE)
    ) %>%
    mutate(tested = ifelse(is.na(tested), FALSE, TRUE))

qu <- wilcox %>%
    arrange(desc(wilcox.p)) %>%
    filter(wilcox.p.adj < 0.1) %>%
    head(1) %>%
    pull(wilcox.p) # Nothing survives q-value correction...

library(ggrepel)
library(ggquantileplot)

# pVolcano <- ggplot(
# ) +
#     geom_point(data = wilcox %>%
#         arrange(tested), aes(x = gFC, y = -log10(wilcox.p), color = tested, alpha = wilcox.p < 0.05)) +
#     geom_text_repel(data = wilcox %>%
#         filter(tested) %>%
#         filter(wilcox.p < 0.05) %>%
#         mutate(species = map_chr(species, \(x) {
#             if (str_detect(species, "s__Coprococcus comes")) {
#                 return("Coprococcus comes")
#             } else if (str_detect(species, "s__Streptococcus sp.")) {
#                 return("Streptococcus mitis")
#             } else {
#                 return(x)
#             }
#         })),
#     aes(x = gFC, y = -log10(wilcox.p), label = species,
#         color = tested), show.legend = F) +
#     scale_color_manual(values = c("TRUE" = "#CD5C5C", "FALSE" = "darkgrey")) +
#     scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
#     theme_classic() +
#     geom_hline(yintercept = -log10(0.05)) +
#     annotate("text", x = -0.25, y = -log10(0.05 * 0.9), label = "unadjusted p = 0.05", size = 3) +
#     xlab("Enrichment in Anaerobic") +
#     ylab('-log10("p-value")') +    
#     NULL

pVolcano <- ggplot(
) +
    geom_point(data = wilcox %>%
        #arrange(tested), aes(x = gFC, y = -log10(wilcox.p), color = tested, alpha = wilcox.p < 0.05)) +
        arrange(tested), aes(x = gFC, y = -log10(wilcox.p), alpha = wilcox.p < 0.05), color = 'darkgrey') +
    # geom_text_repel(data = wilcox %>%
    #     filter(tested) %>%
    #     filter(wilcox.p < 0.05) %>%
    #     mutate(species = map_chr(species, \(x) {
    #         if (str_detect(species, "s__Coprococcus comes")) {
    #             return("Coprococcus comes")
    #         } else if (str_detect(species, "s__Streptococcus sp.")) {
    #             return("Streptococcus mitis")
    #         } else {
    #             return(x)
    #         }
    #     })),
    # aes(x = gFC, y = -log10(wilcox.p), label = species,
    #     color = tested), show.legend = F) +
    #scale_color_manual(values = c("TRUE" = "#CD5C5C", "FALSE" = "darkgrey")) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
    theme_classic() +
    geom_hline(yintercept = -log10(0.05)) +
    annotate("text", x = -0.25, y = -log10(0.05 * 0.9), label = "unadjusted p = 0.05", size = 3) +
    xlab("Enrichment in Anaerobic") +
    ylab('-log10("p-value")') +    
    NULL


 ggsave(plot = (p1Ordination + theme(legend.position = 'bottom') + guides(color = guide_legend(nrow = 3), shape = guide_legend(nrow = 2))) + (pVolcano + theme(legend.position = 'bottom')) + plot_layout(heights = c(1, 1)),
     filename = '/g/scb/zeller/karcher/PRISMA/plots/AAAKVT7HR/ordination_AAAKVT7HR_with_volcano.pdf', width = 6*0.8, height = 11*0.8)

 ggsave(plot = ((p2 + p2.1 + p1.1.1 + p1 + p1.1 + p2.1.1) | p3 | p4 | p5 | p6) + plot_layout(ntrow guides = "collect", widths = c(5.5, 1.5, 1.5, 1.5, 0.5)), filename = "/g/scb/zeller/karcher/PRISMA/plots/AAAKVT7HR/cumrelab_of_tested_strains.pdf", width = 15.5, height = 7.5)
# ggsave(plot = (p2 + ylim(c(0, 50))) + (p2.1 + ylim(c(0, 50))) + (p1 + ylim(c(0, 1))) + (p1.1 + ylim(c(0, 1))) + (p2.1.1 + ylim(c(0, 0.35))) + (p1.1.1 + ylim(c(0, 0.35))) + plot_layout(byrow = TRUE, nrow = 1, guides = "collect"), filename = "/g/scb/zeller/karcher/PRISMA/plots/AAAKVT7HR/cumrelab_of_tested_strains.pdf", width = 12, height = 5)
ggsave(plot = pVolcano, filename = "/g/scb/zeller/karcher/PRISMA/plots/AAAKVT7HR/volcano.pdf", width = 4, height = 4)


#########################################
### Generate code that writes heatmap ###
#########################################

library(tidyverse)
library(pheatmap)
library(data.table)

baseOutPath <- "/g/scb/zeller/karcher/PRISMA/plots/AAAKVT7HR/"
# plot_heatmap <- function(my_df, baseOutPath = NULL) {

#     to_plot_controls <- controls %>%
#         mutate(perc_metabolized_drug = 100 - round(perc_median_correctedAreaNormalized, 0)) %>%
#         filter(Time == 12 & grepl("Avg_Control", Strain)) %>%
#         group_by(Drug, Strain) %>%
#         select(Strain, Drug, perc_metabolized_drug)

#     to_plot_controls$Drug <- as.character(to_plot_controls$Drug)

#     column_factor = to_plot_controls$perc_metabolized_drug
#     names(column_factor) = to_plot_controls$Drug

#     mat_to_plot = pivot_wider(to_plot, names_from = Strain, values_from = perc_metabolized_drug)
#     mat_to_plot = as.data.frame(mat_to_plot)
#     rownames(mat_to_plot) = mat_to_plot$Drug
#     mat_to_plot = mat_to_plot[, -1]
#     mat_to_plot[mat_to_plot < 0] = 0

#     column_factor[column_factor < 0] = 0

#     temp = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(201)
#     temp = unique(temp[(sort(column_factor * 2) + 1)])
#     names(temp) = unique(sort(column_factor))
#     out <- pheatmap::pheatmap(mat_to_plot,
#         cutree_cols = 5,
#         cutree_rows = 8,
#         # cluster_rows = F,
#         # cluster_cols = F,
#         angle_col = 90,
#         fontsize = 13,
#         annotation_row = data.frame(ctrl = column_factor),
#         annotation_colors = list(ctrl = temp), # colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
#         annotation_legend = F,
#         # labels_col=paste0("Donor ", 1:10),
#     )
#     return(out)
# }

#### DF TO USE FOR PLOTTING HEATMAP
for (id in c("SS_AA", "SS_MA")) {
    uniq_min_df_to_plot <- read.csv(str_c("/g/scb/zeller/karcher/PRISMA/data/uniq_min_df_to_plot_", id, ".csv"))

    #### ADJUSTING DF FOT PLOTTING HEATMAP
    uniq_min_df_to_plot <- uniq_min_df_to_plot %>% rename(Sample = abb_strain_name)
    uniq_min_df_to_plot = subset(unique(uniq_min_df_to_plot))
    perc_degr_tmp <- uniq_min_df_to_plot
    uniq_min_df_to_plot$Sample[grepl('Control_', uniq_min_df_to_plot$Sample)] = 'Sterile control'
    uniq_min_df_to_plot <- uniq_min_df_to_plot %>% group_by(Drug, Sample) %>% summarize(FC = median(FC), .groups = "drop")

    uniq_min_df_to_plot$Drug <- factor(uniq_min_df_to_plot$Drug, levels = c("Azathioprine", "Mycophenolate Mofetil", "Sirolimus", "Tacrolimus",
        "Methotrexate", "Cyclosporine", "Everolimus", "Hydrocortisone", "Prednisolone", "Prednisone", "Cortisone",
        "Methylprednisolone", "Betamethasone", "Dexamethasone", "Budesonid", "Mycophenolic Acid",
        "Simvastatin", "Diltiazem", "Amlodipin", "Candesartan", "Ramipril", "Ranitidine", "Atenolol", "Cinacalcet", "Furosemide", "IS_WARFARIN", "IS_CAFFEINE", "IS_IPRIFLAVONE", "IS_LISINOPRIL", "IS_SULFAMETHOXAZOLE"))

    uniq_min_df_to_plot$Sample <- factor(uniq_min_df_to_plot$Sample, levels = c(unique(uniq_min_df_to_plot$Sample)[unique(uniq_min_df_to_plot$Sample) != "Sterile control"], "Sterile control"))

    uniq_min_df_to_plot <- as.data.table(uniq_min_df_to_plot)
    uniq_min_df_to_plot[FC > 1, FC := 1]
    uniq_min_df_to_plot[, Perc_degr := (100 - (FC * 100))]
    uniq_min_df_to_plot[, FC := NULL]

    perc_degr_tmp <- as.data.table(perc_degr_tmp)
    perc_degr_tmp[FC > 1, FC := 1]
    perc_degr_tmp[, Perc_degr := (100 - (FC * 100))]
    perc_degr_tmp[, FC := NULL]

    cast_heat_data = dcast(uniq_min_df_to_plot[, .(Sample, Drug, Perc_degr)], Sample ~ Drug, value.var = "Perc_degr")

    setDF(cast_heat_data)
    rownames(cast_heat_data) = cast_heat_data$Sample
    cast_heat_data = cast_heat_data[, -c(1)]

    # remove IS drugs
    interesting_drugs <- grep("IS_", colnames(cast_heat_data), fixed = T, invert = T, value = T)
    cast_heat_data_red <- cast_heat_data[, colnames(cast_heat_data) %in% interesting_drugs]
    rownames(cast_heat_data_red) <- str_replace_all(rownames(cast_heat_data_red), "  ", " ")
    # cast_heat_data_red <- left_join(
    if(str_detect(id, "MA")) {
        bla <- read_tsv('/g/scb/zeller/karcher/PRISMA/data/itol_tree_leaf_order_MA.txt', col_names = F)
    } else if (str_detect(id, "AA")) {
        bla <- read_tsv('/g/scb/zeller/karcher/PRISMA/data/itol_tree_leaf_order.txt', col_names = F)
    }
    cast_heat_data_red <- inner_join(
        bla %>%
            mutate(X1 = map_chr(X1, \(x){
                a <- str_split(x, " ")[[1]][1]
                a <- str_split(a, "")[[1]][1]
                a <- str_to_upper(a)
                b <- str_split(x, " ")[[1]][2]
                return(str_c(a, ". ", b, sep = ''))
            })) %>%
            rbind(., data.frame(X1 = "Sterile control")) %>%
            rename(species = X1) %>%
            mutate(species = ifelse(species == 'B. longum', "B. longum subsp. infantis", species)),
        cast_heat_data_red %>%
            rownames_to_column('species') %>%
            mutate(species = ifelse(species == 'A. naeslundi', "A. naeslundii", species)) %>%
            as_tibble()
    ) %>%
        column_to_rownames('species') %>%
        as.data.frame() %>%
        as.matrix()


    library(viridis)
    # pdf("/g/scb/zeller/karcher/PRISMA/plots/AAAKVT7HR/Heatmap_SS_MA.pdf", width = 19, height = 14)
    if (id == "SS_MA") {
        pdf(str_c("/g/scb/zeller/karcher/PRISMA/plots/AAAKVT7HR/Heatmap_", id, ".pdf"), width = 19, height = 7)
    } else if (id == "SS_AA") {
        pdf(str_c("/g/scb/zeller/karcher/PRISMA/plots/AAAKVT7HR/Heatmap_", id, ".pdf"), width = 19, height = 14)
    } else {
        asdadad
    }

    pheatmap(cast_heat_data_red, cluster_rows = F, cluster_cols = F, angle_col = 45, fontsize = 16, color = viridis(12))
    dev.off()

}

###########################################
### Generate genus-level tax-tree ###
###########################################

taxs <- referenceProfiles %>%
    ungroup() %>%
    select(kingdom, phylum, class, order, family, genus, species, mOTU_ID) %>%
    distinct()
colnames(taxs) <- c("kingdom", 'phylum', 'class', 'order', 'family', 'genus', 'species', 'mOTU_ID')
taxs <- taxs %>%
    as.data.frame() %>%
    as_tibble() %>%
    # distinct() %>%
    mutate(mOTU_ID = map_chr(mOTU_ID, function(x) str_replace(x, "_v31", "_v3")))

taxs <- taxs %>%
    select(kingdom, phylum, class, order, family, genus, species) %>%
    group_by(kingdom, phylum, class, order, family, genus, species) %>%
    tally() %>%
    arrange(desc(n)) %>%
    select(-n) %>%
    ungroup() %>%
    distinct(genus, .keep_all = T)

taxs <- taxs[apply(taxs, 1, function(x) !any(is.na(x))), ]
taxs <- taxs %>%
    # filter(!str_detect(family_GTDB, "annotated")) %>%
    filter(!str_detect(genus, "annotated")) %>%
    filter(!str_detect(genus, "sedis")) %>%
    filter(genus != '') %>%
    distinct()

### CHECK THIS OUT
# Here I write the tax file and manually clean it up and add the few tested species
# that we can't extract from the mOTUs3 taxonomy

########### IMPORANT #################
# Upong reading the tree into itol, you have to disable the "Leaf sorting" option in the "Advanced" Reiter
######################################

library(ape)
library(RColorBrewer)
library(tidyverse)
library(taxonomizr)

taxs <- read_tsv('/g/scb/zeller/karcher/PRISMA/data/AAAKVT7HR_tree_tax.tsv', col_names = F) %>%
    filter(!str_detect(X7, "graevenitzii")) %>%
    mutate(X7 = str_replace(X7, "s__", "")) %>%
    # Bring tree into desired order...
    inner_join(read_tsv('/g/scb/zeller/karcher/PRISMA/data/itol_tree_leaf_order.txt', col_names = F) %>%
        rename(X7 = X1), .) %>%
    mutate(X7 = str_replace_all(X7, " ", "_________")) %>%
    relocate(X1, X2, X3, X4, X5, X6, X7)

tree <- makeNewick(taxs %>%
    as.matrix())
fileConn <- file("/g/scb/zeller/karcher/PRISMA/results/AAAKVT7HR/tree.genus.ncbi.nwk")
writeLines(tree, fileConn)
close(fileConn)

tree.raw <- ape::read.tree('/g/scb/zeller/karcher/PRISMA/results/AAAKVT7HR/tree.genus.ncbi.nwk')
write.tree(phy = tree.raw, file = "/g/scb/zeller/karcher/PRISMA/results/AAAKVT7HR/AAAKVT7HR.tree.genus.ncbi.filtered.nwk")

tree.filtered <- read.tree("/g/scb/zeller/karcher/PRISMA/results/AAAKVT7HR/AAAKVT7HR.tree.genus.ncbi.filtered.nwk")
tree.filtered$tip.label <- str_replace_all(tree.filtered$tip.label, "_________", " ")
tree.filtered$tip.label <- str_replace_all(tree.filtered$tip.label, "_", " ")
write.tree(phy = tree.filtered, file = "/g/scb/zeller/karcher/PRISMA/results/AAAKVT7HR/AAAKVT7HR.tree.genus.ncbi.filtered.nwk")

# ... aaaand for MA I'm actually doing it manually
# It's here anyway

taxs <- read_tsv('/g/scb/zeller/karcher/PRISMA/data/AAAKVT7HR_tree_tax.tsv', col_names = F) %>%
    #filter(!str_detect(X7, "graevenitzii")) %>%
    mutate(X7 = str_replace(X7, "s__", "")) %>%
    # Bring tree into desired order...
    inner_join(read_tsv('/g/scb/zeller/karcher/PRISMA/data/itol_tree_leaf_order_MA.txt', col_names = F) %>%
        rename(X7 = X1), .) %>%
    mutate(X7 = str_replace_all(X7, " ", "_________")) %>%
    relocate(X1, X2, X3, X4, X5, X6, X7)

tree <- makeNewick(taxs %>%
    as.matrix())
fileConn <- file("/g/scb/zeller/karcher/PRISMA/results/AAAKVT7HR/tree.genus.ncbi.nwk")
writeLines(tree, fileConn)
close(fileConn)

tree.raw <- ape::read.tree('/g/scb/zeller/karcher/PRISMA/results/AAAKVT7HR/tree.genus.ncbi.nwk')
write.tree(phy = tree.raw, file = "/g/scb/zeller/karcher/PRISMA/results/AAAKVT7HR/AAAKVT7HR.tree.genus.ncbi.filtered_MA.nwk")


tree.filtered <- read.tree("/g/scb/zeller/karcher/PRISMA/results/AAAKVT7HR/AAAKVT7HR.tree.genus.ncbi.filtered_MA.nwk")
tree.filtered$tip.label <- str_replace_all(tree.filtered$tip.label, "_________", " ")
tree.filtered$tip.label <- str_replace_all(tree.filtered$tip.label, "_", " ")
write.tree(phy = tree.filtered, file = "/g/scb/zeller/karcher/PRISMA/results/AAAKVT7HR/AAAKVT7HR.tree.genus.ncbi.filtered_MA.nwk")


taxa <- tree.filtered$node.label
unique_phyla <- taxa[str_detect(taxa, "p__")]
colors <- colorRampPalette(brewer.pal(n = length(unique_phyla),
    "Pastel2"))(length(unique_phyla))
set.seed(6)
colors <- sample(colors, size = length(unique_phyla))

color_data <- data.frame(phyla = unique_phyla,
    colors = colors)
color_data$colors <- c("#b3e2cd",
    "#fdcdac",
    "#cbd5e8",
    "#f4cae4",
    "#e6f5c9")
# color_data$colors[color_data$phyla == "p__Actinobacteria"] <- "#8286ed"
# color_data$colors[color_data$phyla == "p__Tenericutes"] <- "#a65057"

# 1. Phylum-level labeling of inner nodes. Keep black and later on label in illustrator.
file.copy(from = '/g/scb2/zeller/karcher/CAZY_project_v2/scripts/generate_figure_1/template_files/dataset_symbols_template.txt',
    '/g/scb/zeller/karcher/PRISMA/results/AAAKVT7HR/itol_files_populated/phylum_symbols_variation.txt', overwrite = T)
filePath <- '/g/scb/zeller/karcher/PRISMA/results/AAAKVT7HR/itol_files_populated/phylum_symbols_variation.txt'
tx <- readLines(filePath)
tx2 <- gsub(pattern = "MAXIMUM_SIZE,50", replace = "MAXIMUM_SIZE,25", x = tx)
writeLines(tx, con = filePath)

for (i in 1:dim(color_data)[1]) {
    phylum <- color_data$phyla[i]
    write(str_c(phylum, "2", 10, "#000000", "1", "1", sep = ",", collapse = ","), filePath, append = T)
}

# 1: Phylum-level ranges indicating phylum membership
file.copy(from = '/g/scb2/zeller/karcher/CAZY_project_v2/scripts/generate_figure_1/template_files/dataset_ranges_template.txt',
    '/g/scb/zeller/karcher/PRISMA/results/AAAKVT7HR/itol_files_populated/phylum_ranges.txt', overwrite = T)
filePath <- '/g/scb/zeller/karcher/PRISMA/results/AAAKVT7HR/itol_files_populated/phylum_ranges.txt'
for (i in 1:dim(color_data)[1]) {
    phylum <- color_data$phyla[i]
    col <- color_data$colors[i]
    write(str_c(phylum, ",", phylum, ",", col, ",", col, ",#000000,dashed,0,,#0000ff,0,italic"), filePath, append = T)
}

# 2 Leaf nodes scaled to the genus mean relative abundance
file.copy(from = '/g/scb2/zeller/karcher/CAZY_project_v2/scripts/generate_figure_1/template_files/dataset_symbols_template.txt',
    '/g/scb/zeller/karcher/PRISMA/results/AAAKVT7HR/itol_files_populated/meanGenusAbundance_symbols_variation.txt', overwrite = T)
filePath <- '/g/scb/zeller/karcher/PRISMA/results/AAAKVT7HR/itol_files_populated/meanGenusAbundance_symbols_variation.txt'
tx <- readLines(filePath)
tx2 <- gsub(pattern = "MAXIMUM_SIZE,50", replace = "MAXIMUM_SIZE,25", x = tx)
writeLines(tx, con = filePath)

# tmp <- referenceProfiles %>%
#     # group_by(sampleID, genus) %>%
#     # summarize(relAb = sum(relAb)) %>%
#     group_by(species) %>%
#     summarize(meanRelAb = mean(relAb)) %>%
#     #mutate(genus = str_replace(genus, "g__", "")) %>%
#     mutate(species = str_replace(species, "s__", "")) %>%
#     filter(species %in% tree.filtered$tip.label)



tmp <- testedSpeciesRealMetaG %>%
    ungroup() %>%
    mutate(bacterium = ifelse(bacterium == "Eubacterium hallii\u00a0", "Eubacterium hallii", bacterium)) %>%
    mutate(bacterium = ifelse(bacterium == "Escherichia coli BW25113", "Escherichia coli", bacterium)) %>%
    mutate(bacterium = ifelse(bacterium == "Actinomyces naeslundi", "Actinomyces naeslundii", bacterium)) %>%
    mutate(bacterium = ifelse(bacterium == "Clostridium ramosum\u00a0", "Clostridium ramosum", bacterium)) %>%
    mutate(bacterium = ifelse(bacterium == "Streptococcus salivarius\u00a0", "Streptococcus salivarius", bacterium)) %>%
    select(kingdom, phylum, species, bacterium, sampleID, relAb) %>%
    group_by(phylum, bacterium) %>%
    summarize(relAb = mean(relAb)) %>%
    filter(bacterium %in% tree.filtered$tip.label) %>%
    mutate(relAb = 1 / log10(relAb + pseudoCount)) %>%
    mutate(bacterium = str_replace_all(bacterium, " ", "_"))

for (i in 1:dim(tmp)[1]) {
    family <- tmp$bacterium[i]
    val <- tmp$relAb[i[]]
    # write(str_c(family, "2", 1/-log10(val), "#000000", "1", "1",  sep = ",", collapse = ","), filePath, append = T)
    write(str_c(family, "2", val, "#000000", "1", "1", sep = ",", collapse = ","), filePath, append = T)
}

# Similarly as above, get the values from remote R session,
# then execute this code in a local R environment to get a legend for the heatmap

library(ggplot2)
ggplot(data = data.frame(x = 1:2, y = c(log10(0.00000748),
    log10(0.149)))) +
    geom_point(aes(x = x, y = y, size = y)) +
    # scale_radius(breaks = c(-0.01, -0.1)) +
    scale_radius(breaks = c(-4, -3, -2, -1, log10(0.149)), range = c(1, 7.5)) +
    theme_classic()
# scale_fill_gradientn(colours = c("#0000FF", "#FFFFFF", "#FF0000"), values = rescale(c(-1.447, 0, 5.5)), breaks = c(-1, 0, 1, 3, 5))


# 3. Bar-chart with the mOTUs/genus
file.copy(from = '/g/scb2/zeller/karcher/CAZY_project_v2/scripts/generate_figure_1/template_files/dataset_simplebar_template.txt',
    '/g/scb/zeller/karcher/PRISMA/results/AAAKVT7HR/itol_files_populated/numbergenomes_bar.txt', overwrite = T)
filePath <- '/g/scb/zeller/karcher/PRISMA/results/AAAKVT7HR/itol_files_populated/numbergenomes_bar.txt'
tx <- readLines(filePath)
# tx2  <- gsub(pattern = "abc", replace = "ccccccccccccccccccccc", x = tx)
# tx <- str_replace(tx, "#DATASET_SCALE,2000,10000,20000", "DATASET_SCALE,1,2,3,4,5")
tx <- str_replace(tx, "#DATASET_SCALE,2000,10000,20000", "DATASET_SCALE,0-0-#000000-2-1-3, DATASET_SCALE,5-5-#000000-2-1-3,10-10-#000000-2-1-3")
tx <- str_replace(tx, "DATASET_LABEL,label 1", "DATASET_LABEL, Number of tested species")
tx <- str_replace(tx, "COLOR,#ff0000", "COLOR,#000000")
writeLines(tx, con = filePath)

plotData <- read_csv('/g/scb/zeller/karcher/PRISMA/data/WGS_metadata/straintaxid_genustaxid.csv') %>%
    # We need some manual fixes due to differences in taxonomy (mOTUs3 calls Clostridium bolteae, but NCBI tax is actually Enterocloster bolteae)
    mutate(genusName = ifelse(species == "Enterocloster bolteae", 'Clostridium', genusName)) %>%
    # Anaerobutyricum hallii -> Eubacterium hallii
    mutate(genusName = ifelse(species == "Anaerobutyricum hallii", 'Eubacterium', genusName)) %>%
    # ...
    mutate(genusName = ifelse(species == "Lacticaseibacillus paracasei", 'Lactobacillus', genusName)) %>%
    mutate(genusName = ifelse(species == "Ligilactobacillus ruminis", 'Lactobacillus', genusName)) %>%
    mutate(genusName = ifelse(species == "Limosilactobacillus reuteri", 'Lactobacillus', genusName)) %>%
    mutate(genusName = ifelse(species == "Limosilactobacillus fermentum", 'Lactobacillus', genusName)) %>%
    mutate(genusName = ifelse(species == "Segatella copri", 'Prevotella', genusName)) %>%
    mutate(genusName = ifelse(species == "Thomasclavelia ramosa", 'Clostridium', genusName)) %>%
    # second round (might be fucked)
    mutate(genusName = ifelse(species == "Capnocytophaga ochracea", 'Bacteroides', genusName)) %>%
    mutate(genusName = ifelse(species == "Phocaeicola vulgatus", 'Bacteroides', genusName)) %>%
    select(genusName) %>%
    rename(genus = genusName) %>%
    group_by(genus) %>%
    tally()

for (i in 1:dim(plotData)[1]) {
    family <- plotData$genus[i]
    value <- plotData$n[i]
    write(str_c(family, ",", value), filePath, append = T)
}
