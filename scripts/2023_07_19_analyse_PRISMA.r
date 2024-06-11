library(tidyverse)
library(ggplot2)
library(patchwork)
library(vegan)
library(ggrepel)
library(RColorBrewer)
library(eeptools)
library(readxl)

setwd('/g/scb/zeller/karcher/PRISMA/scripts')
source('/home/karcher/utils/utils.r')

# Define some convenience functions

labelLink <- c(
    "1" = "pre-transplant (1)",
    "2" = "pre-transplant (2)",
    "3" = "post-immunosuppression\npre-transplant",
    "4" = "week 1 post-transplant",
    "5" = "week 4 post-transplant",
    "6" = "Month 3 post-transplant",
    "7" = "Month 6 post-transplant"
)

scale_x_discrete_prisma <- function(labelMap = labelLink, how = 'discrete', ...) {
    if (how == 'discrete') {
        scale_x_discrete(labels = labelMap, ...)
    } else if (how == "continuous") {
        scale_x_continuous(labels = labelMap, ...)
    }
}

rarefactionDepth <- 1E4
pseudoCount <- 1E-4

basedir <- '/g/scb/zeller/karcher/PRISMA/'
setwd(basedir)

profiles <- readRDS('/g/scb/zeller/karcher/PRISMA/profiles/16S/KLGPG_221206_PRISMA_MG/res_mapseq.rds')
meta <- read_tsv("/g/scb/zeller/karcher/PRISMA/data/16S_metadata/221227_PRISMA_16S_Overview.tsv") %>%
    rename(sampleID = Sample_ID, visit = Visit) %>%
    mutate(PSN = str_replace(PSN, "Mue", "M")) %>%
    mutate(PSN = str_replace(PSN, "NTXM", "NZMU")) %>%
    mutate(PSN = str_replace(PSN, "-0", "-"))
# Remove funny colnames with encoding
meta <- meta[, !str_detect(colnames(meta), "Conc")]
meta <- meta[, !str_detect(colnames(meta), "10 ng")]
fullTax <- read_tsv("/g/scb/zeller/karcher/PRISMA/data/MAPseq_AlessioCurated.tax", col_names = F)
fullTax <- fullTax %>%
    mutate(X6 = map_chr(X6, \(x) {
        if (str_detect(x, "\\[Eubacterium\\]")) {
            x <- "Eubacterium"
        } else {
            x <- x
        }
        return(x)
    }))
colnames(fullTax) <- c(
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
)
fullTax <- fullTax %>%
    select(-species) %>%
    distinct() %>%
    group_by(genus) %>%
    nest() %>%
    mutate(data = map2(data, genus, \(x, g) {
        if (g == "Eubacterium") {
            return(x %>% filter(family == "Lachnospiraceae"))
        } else {
            return(x)
        }
    })) %>%
    unnest()


profiles <- profiles[!rownames(profiles) == "Bacteria", ]
colnames(profiles) <- map_chr(colnames(profiles), function(x) str_replace(x, ".*lane1", ""))
profiles <- profiles[, !str_detect(colnames(profiles), "water")]
colnames(profiles) <- map_chr(colnames(profiles), function(x) str_replace(x, "MG", "MG_"))

rownames(profiles) <- map_chr(rownames(profiles), \(x) {
    if (str_detect(x, "\\[Eubacterium\\]")) {
        x <- "Eubacterium"
    } else {
        x <- x
    }
    return(x)
})

profiles <- profiles %>%
    as.data.frame() %>%
    mutate(genus = rownames(.)) %>%
    pivot_longer(-genus) %>%
    rename(sampleID = name, count = value) %>%
    group_by(genus, sampleID) %>%
    summarize(count = sum(count))

depths <- profiles %>%
    group_by(sampleID) %>%
    summarize(totalGenusReadCount = sum(count)) %>%
    left_join(meta, by = 'sampleID')

depth_histo <- ggplot(depths %>%
    mutate(visit = as.factor(visit))) +
    theme_classic() +
    geom_histogram(aes(x = totalGenusReadCount), alpha = 1, position = "identity") +
    geom_vline(xintercept = rarefactionDepth) +
    geom_text(data = depths %>%
        mutate(visit = as.factor(visit)) %>%
        group_by(visit) %>%
        filter(totalGenusReadCount < rarefactionDepth) %>%
        tally(), aes(x = rarefactionDepth * 0.9, y = 5, label = n)) +
    facet_grid(visit ~ .) +
    ylab("N samples") +
    scale_y_continuous(breaks = c(1, 3, 5, 7, 9, 11))

ggsave(plot = depth_histo, filename = "plots/KLGPG_221206/depth_histogram.pdf", width = 5, height = 7)


highDepthSamples <- depths %>%
    mutate(visit = as.factor(visit)) %>%
    group_by(sampleID) %>%
    filter(totalGenusReadCount >= rarefactionDepth) %>%
    select(sampleID)

profiles <- profiles %>%
    inner_join(highDepthSamples)

# rarefy
profiles <- profiles %>%
    pivot_wider(id_cols = sampleID, names_from = genus, values_from = count) %>%
    column_to_rownames("sampleID") %>%
    rrarefy(sample = rarefactionDepth) %>%
    as.data.frame() %>%
    rownames_to_column('sampleID') %>%
    pivot_longer(-sampleID) %>%
    rename(genus = name, count = value)

profilesCountsWithUnresolved <- profiles

profiles <- profiles %>%
    group_by(sampleID) %>%
    mutate(relAb = count / sum(count)) %>%
    select(-count) %>%
    filter(genus != "not_resolved")

cumRelAbToGenus <- profiles %>%
    left_join(meta, by = 'sampleID') %>%
    group_by(sampleID, visit) %>%
    summarize(sumRelAb = sum(relAb))

cumRelAbToGenus_box <- ggplot(cumRelAbToGenus %>%
    mutate(visit = as.factor(visit))) +
    theme_classic() +
    geom_boxplot(aes(x = visit, y = sumRelAb)) +
    ylab("Fraction reads\nannotated to Genus-level") +
    theme_classic()

ggsave(plot = cumRelAbToGenus_box, filename = "plots/KLGPG_221206/cumRelAb_boxplot.pdf", width = 4, height = 3.5)

profiles <- profiles %>%
    mutate(relAbOrig = relAb) %>%
    mutate(relAb = log10(relAb + pseudoCount)) %>%
    left_join(meta, by = "sampleID") %>%
    left_join(fullTax, by = 'genus')

profiles <- profiles %>%
    # mutate(visitMinusOne = visit - 1) %>%
    group_by(PSN, visit) %>%
    nest() %>%
    group_by(PSN) %>%
    nest() %>%
    mutate(data = map(data, function(x) {
        x <- x %>%
            mutate(visit = as.numeric(as.character(visit))) %>%
            arrange(visit) %>%
            mutate(visitMinusOne = c(NA, visit[1:(length(visit) - 1)]))
        return(x)
    })) %>%
    unnest() %>%
    unnest() %>%
    relocate(sampleID, genus, relAb, relAbOrig, PSN, visit, visitMinusOne)

pairwiseDistances <- pivot_wider(profiles, id_cols = genus, names_from = sampleID, values_from = relAb) %>%
    column_to_rownames("genus") %>%
    as.data.frame() %>%
    as.matrix() %>%
    t() %>%
    vegdist(method = "euclidean", k = 2)

pairwiseDistancesIdentityEuclidean <- pivot_wider(profiles, id_cols = genus, names_from = sampleID, values_from = relAbOrig) %>%
    column_to_rownames("genus") %>%
    as.data.frame() %>%
    as.matrix() %>%
    t() %>%
    vegdist(method = "euclidean", k = 2)

pcoa <- cmdscale(pairwiseDistances) %>%
    as.data.frame() %>%
    rownames_to_column("sampleID") %>%
    left_join(meta, by = "sampleID") %>%
    as_tibble() %>%
    mutate(visit = as.character(visit)) %>%
    mutate(visit = factor(visit, levels = 1:7))

# for (i in seq_along(names(labelLink))) {
#     pre <- names(labelLink)[i]
#     repl <- labelLink[i]
#     pcoa$visit <- ifelse(pcoa$visit == pre, repl, pcoa$visit)
# }

# pcoa <- pcoa %>%
#     mutate(visit = factor(visit, levels = labelLink))


pcoa_plot <- ggplot() +
    geom_point(data = pcoa, aes(x = V1, y = V2, color = visit)) +
    theme_classic() +
    xlab("PCo 1") +
    ylab("PCo 2") +
    scale_color_manual(values = c(
        "#5b92fc",
        "#5b92fc",
        "#3579fc",
        "#ff3636",
        "#f75e5e",
        "#facaca",
        "#fcebeb"))

ggsave(plot = pcoa_plot, filename = "plots/KLGPG_221206/pcoa_visit_v1.pdf", width = 5.5, height = 3.5)


pcoa_plot <- ggplot() +
    geom_point(data = pcoa %>% mutate(visit = factor(visit, levels = 1:7)), aes(x = V1, y = V2, color = visit)) +
    theme_classic() +
    xlab("PCo 1") +
    ylab("PCo 2") +
    scale_color_manual(values = c(
        "#C20000",
        "#FF0000",
        "#B1B1B4",
        "#00CDFF",
        "#00AFDA",
        "#0091B5",
        "#010A4E"))



add_line_with_ind <- function(pObj = NULL, dataB = NULL, indName = NULL) {
    dataB <- dataB %>%
        mutate(visitAr = as.numeric(as.character(visit))) %>%
        filter(PSN == indName) %>%
        arrange(visitAr)

    dataBSegment <- dataB %>%
        select(V1, V2, PSN, visit) %>%
        mutate(V1End = c(V1[2:length(V1)], V1[1])) %>%
        mutate(V2End = c(V2[2:length(V2)], V2[1]))

    pObj <- pObj +
        # geom_path(data = dataB, aes(x = V1, y = V2, group = PSN), color = 'black', alpha = 0.5, arrow = arrow(angle = 15, ends = "both", type = "closed"))  +
        # geom_path(data = dataB, aes(x = V1, y = V2, group = PSN, fill = 1), arrow = arrow(), color = 'black', alpha = 0.5)  +
        geom_segment(data = dataBSegment, aes(x = V1, xend = V1End, y = V2, yend = V2End), arrow = arrow(length = unit(0.01, "npc"))) +
        geom_label_repel(data = dataB, aes(x = V1, y = V2, label = as.numeric(visit)), alpha = 0.75)
    return(pObj)
}

for (indName in unique(pcoa$PSN)) {
    # for (indName in c('RiGr-NZHD-8')) {
    ggsave(
        plot = add_line_with_ind(ggplot(), pcoa, indName) +
            geom_point(data = pcoa, aes(x = V1, y = V2, color = visit)) +
            theme_classic() +
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
        filename = str_c("plots/KLGPG_221206/single_ind_pcoa_plots/pcoa_visit_v3_visit_3_removed_treatment_dichotomized_", indName, ".pdf", collapse = ""), width = 5.5, height = 3.5)
}

nice_looking_individuals <- c(
    "AbWa-KKHD-12",
    "RaCu-NZHD-17",
    "BaEr-NZHD-39",
    "DeZi-KKHD-2",
    "AnCz-KKHD-11",
    'MaBa-NZHD-10'
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
        theme_classic()
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
    for (indName in c("AbWa-KKHD-12")) {
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
            filename = str_c("plots/KLGPG_221206/single_ind_tax_barplots/taxa_barplots_visit_v3_visit_3_removed_treatment_dichotomized_", indName, "__", taxL, ".pdf", collapse = ""),
            width = 5.75,
            height = 3
        )
    }

}

library(tidyverse)

tmp <- profiles %>%
    filter(genus %in% c(
        "Bacteroides",
        "Eubacterium",
        "Ruinococcus",
        "Dorea",
        "Clostridium",
        "Prevotella"
    ))

genus_boxplot <- ggplot() +
    geom_boxplot(data = tmp %>%
        mutate(visit = factor(as.character(visit), levels = 1:7)), aes(x = genus, y = relAb, fill = visit)) +
    theme_classic()

ggsave(plot = genus_boxplot, filename = "plots/KLGPG_221206/genus_boxplots.pdf", width = 7.25, height = 4)

genus_boxplot <- ggplot() +
    geom_boxplot(data = tmp %>%
        filter(visit != 3) %>%
        mutate(visit = factor(as.character(visit), levels = 1:7)), aes(x = genus, y = relAb, fill = visit)) +
    theme_classic()

ggsave(plot = genus_boxplot, filename = "plots/KLGPG_221206/genus_boxplots_v2_visit_3_removed.pdf", width = 7.25, height = 4)

genus_boxplot <- ggplot() +
    geom_boxplot(data = tmp %>%
        filter(visit != 3) %>%
        mutate(treatment = case_when(visit <= 2 ~ "before\ntransplantation",
            visit > 2 ~ "after\ntransplanation")), aes(x = genus, y = relAb, fill = treatment)) +
    theme_classic()

ggsave(plot = genus_boxplot, filename = "plots/KLGPG_221206/genus_boxplots_v3_visit_3_removed_treatment_dichotomized.pdf", width = 6, height = 4)

#################################################################################
# 221024: Integrate clinical metadata and do first analysis of interims cohort
#################################################################################


# visits are uniquely identified by group_by(v65_pat_id, v62_visit_number)
outcomeInformation <- read_csv('/g/scb/zeller/karcher/PRISMA/data/16S_metadata/221024_PRISMA_clinical_metadata.csv')
# outcomeInformation <- read_csv('/g/scb/zeller/karcher/PRISMA/data/16S_metadata/221024_PRISMA_clinical_metadata_names_fixed_ACTUALLY_NEVERMIND_JUST_DO_IT_YOURSELF.csv')
clinMetCode <- read_tsv('/g/scb/zeller/karcher/PRISMA/data/16S_metadata/231024_PRISMA_clinical_metadata_codebook.tsv')

# this is the tibble containing
model_covariates <- read_excel('/g/scb/zeller/karcher/PRISMA/data/16S_metadata/covariable_columns.xlsx') %>%
    # ATTENTION: This will have to change at a later stage
    filter(C_D_Ratio_Relevance_prio)
tmp <- model_covariates$Column_newName
model_covariates <- model_covariates$Column
names(model_covariates) <- tmp

outcomeInformation <- outcomeInformation %>%
    rename(
        patientID = v65_pat_id,
        visitNumber = v62_visit_number,
        # birthday = v13_dob,
        ##################################
        # clinical outcome related stuff
        ##################################
        hospitalization = V30a_hospitalization,
        changeImmunosuppRegimen = v138_change_immunosupp,
        rejection = v97_rejection,
        ## Neither graft-loss nor death happened in initial 30-membre cohort
        graftLoss = v22a_graft_loss,
        death = v23a_death,
        ## stuff to calculate CD ratio
        ## Every patient got Tacrolimus but not everyone got the same preparation. So loop over groups and take non-NA values
        ### doses
        tacDosePrograf = v112b_Tacrolimus_Prograf_dose,
        tacDoseEnvarsus = v113b_Tacrolimus_Envarsus_dose,
        tacDoseAdvagraf = v114b_Tacrolimus_Advagraf_dose,
        tacDoseModigraf = v404b_Tacrolimus_Modigraf_dose,
        ### concentrations
        tacConcentration = v254_Tacrolimus,
        ### study center
        studyCenter = study_center,
        # cyp3a5star3 = v4_CYP3A5_3,
        # ATTENTION: cyp3a4star22 is (all == FALSE), in model building cohert
        # cyp3a4star22 = v64_CYP3A4_22,
    ) %>%
    rename(all_of(model_covariates)) %>%
    mutate(patientID = str_replace(patientID, "-0", '-')) %>%
    mutate(birthday = as.Date(birthday)) %>%
    mutate(age = age_calc(birthday, as.Date("2023-11-07"), "years")) %>%
    mutate(ageCategorical = ifelse(age > 18, "adult", 'non-adult'))
## [1]v4_cyp_genotype =  ";\"\";1;\"CYP3A5(*3) Positive\";2;\"CYP3A5(*3) Negative\""
## [2]v64_cyp_genotype_2 = ";\"\";1;\"CYP3A4(*22) Positive\";2;\"CYP3A4(*22) Negative\""
### Update: Information in codebook is wrong. They are already coded as True/False
# cyp3a5star3 = case_when(
#     cyp3a5star3 == 1 ~ "positive",
#     cyp3a5star3 == 2 ~ "negative",
#     .default = "unknown"
# ),
# cyp3a4star22 = case_when(
#     cyp3a4star22 == 1 ~ "positive",
#     cyp3a4star22 == 2 ~ "negative",
#     .default = "unknown"
# )
# abxInfo <- read.table('/g/scb/zeller/karcher/PRISMA/data/16S_metadata/221030_PRISMA_antiinfectives_metadata_v2.csv', sep = ";", header =TRUE)
abxInfo <- read.table('/g/scb/zeller/karcher/PRISMA/data/16S_metadata/antiinfectives_metadata.csv', comment.char = "#", sep = ",", header = TRUE) %>%
    ########################################################################
    # IMPORTANT: I completely disregard Cotrimoxazol, Cefriaxon, Nystatin
    ########################################################################
    filter(!drug_name %in% c("v160a_Cotrimoxazol", "v188a_Nystatin", "v155a_Ceftriaxon"))


tmp <- c()
for (abx in abxInfo$drug_name) {
    abxNewName <- str_split_fixed(abx, "_", n = 3)[2:length(str_split_fixed(abx, "_", n = 3))]
    abxNewName <- str_c(abxNewName, sep = ".", collapse = ".")
    abxNewName <- str_replace(abxNewName, "[.]$", "")
    if (sum(colnames(outcomeInformation) == abx) == 0) {
        print(str_c("Cannot find column ", abx, ' in metadata file.'))
    } else if (sum(colnames(outcomeInformation) == abx) > 1) {
        print("More than one column name is indetical??")
        exit()
    }
    colnames(outcomeInformation)[colnames(outcomeInformation) == abx] <- abxNewName
    tmp <- c(tmp, abxNewName)
}
abxInfo$allAbx <- tmp
abxInfo <- abxInfo %>%
    mutate(subclass = subclass_2)
abxInfo <- abxInfo
# mutate(subclass = ifelse(subclass %in% (abxInfo %>% group_by(subclass) %>% tally() %>% filter(n>=3) %>% pull(subclass)), subclass, "miscellaneous"))
abxInfo$X1 <- NULL
outcomeInformation <- outcomeInformation[, !str_detect(colnames(outcomeInformation), "_")]

# Merge dose columns to have only one meaningful one and then calc CD ratio
outcomeInformation <- outcomeInformation %>%
    mutate(finTacDose = pmap_dbl(list(tacDosePrograf, tacDoseEnvarsus, tacDoseAdvagraf, tacDoseModigraf), function(a, b, c, d) {
        tmp <- c(a, b, c, d)
        if (all(is.na(tmp))) {
            return(NA)
        }
        # stopifnot(sum(!is.na(tmp)) == 1)
        return(tmp[!is.na(tmp)])
    })) %>%
    select(-all_of(colnames(.)[str_detect(colnames(.), 'tacDose')])) %>%
    mutate(CD = tacConcentration / finTacDose)
stopifnot(all(abxInfo$allAbx == outcomeInformation %>% select(all_of(abxInfo$allAbx)) %>% colnames()))
# Same for ABx
####################################
# For ABx, I interpret NAs as FALSE
####################################
outcomeInformation$anyABx <- apply(
    outcomeInformation %>%
        select(all_of(abxInfo$allAbx)),
    1,
    \(x) return(any(ifelse(is.na(x), FALSE, x)))
)
# Generate ABx type link
outcomeInformation$ABxSubClass <- apply(
    outcomeInformation %>%
        select(all_of(abxInfo$allAbx)),
    1,
    # \(x) return(any(ifelse(is.na(x), FALSE, x)))
    \(x) {
        if (all(is.na(x))) {
            # If all are NA, that means probably visit1/2
            return("none")
        } else {
            if (all(!x)) {
                return("none")
            } else {
                tmp <- str_c(abxInfo$subclass[which(x)])
                # print('a')
                # If some entries of tmp are NA, there are ABx,s taken that do not correspond to any of the major subclasses we defined
                if (any(is.na(tmp))) {
                    types <- c('others')
                } else {
                    types <- c()
                }
                # print(types)
                tmp <- tmp[!is.na(tmp)]
                # print('c')
                # print(c(types, tmp))
                return(str_c(sort(unique(c(types, tmp))), sep = ',', collapse = ','))

            }
        }
    }
)

outcomeInformation <- outcomeInformation %>%
    # Finally, set 'others' and 'none' to NA
    mutate(ABxSubClass = ifelse(ABxSubClass %in% c("none"), NA, ABxSubClass))

outcomeInformation <- outcomeInformation %>%
    select(-all_of(abxInfo$allAbx))

outcomeInformation <- outcomeInformation %>%
    mutate(postTransplant = visitNumber >= 4)

outcomeInformation <- outcomeInformation %>%
    mutate(across(c(hospitalization, rejection, changeImmunosuppRegimen), \(x) ifelse(is.na(x), FALSE, x)))

outcomeInformation <- outcomeInformation %>% mutate(CDbinary = factor(ifelse(CD > 1, "high", "low"), levels = c('low', 'high')))

##################################
##################################
############ IMPORTANT ###########
##################################
##################################

# clinicalMetadata summarizes, well, clinical metadata. i.e. clinical model
# it will NOT contain outcomes (those will be)
clinicalMetadata <- outcomeInformation %>%
    group_by(patientID) %>%
    nest() %>%
    mutate(data = map(data, function(x) {
        x <- x %>%
            mutate(anyComplication = pmap_lgl(list(hospitalization, rejection, changeImmunosuppRegimen), function(a, b, c) {
                return(any(c(a, b, c)))
            }))
        return(x)
    })) %>%
    mutate(firstVisitWithComplication = map_int(data, function(x) {
        x <- x %>%
            mutate(anyComplication = pmap_lgl(list(hospitalization, rejection, changeImmunosuppRegimen), function(a, b, c) {
                return(any(c(a, b, c)))
            }))
        if (any(x$anyComplication)) {
            return(
                x %>%
                    filter(anyComplication) %>%
                    arrange(visitNumber) %>%
                    head(1) %>%
                    pull(visitNumber)
            )
        } else {
            return(NA)
        }
    })) %>%
    mutate(firstAlbuminMeasurement = map_dbl(data, function(x) {
        x <- x %>%
            filter(!is.na(albumin)) %>%
            arrange(visitNumber) %>%
            head(1) %>%
            pull(albumin)
        if (length(x) == 0) {
            return(NA)
        } else {
            return(x)
        }
    })) %>%
    mutate(firstHematocritMeasurement = map_dbl(data, function(x) {
        x <- x %>%
            filter(!is.na(hematocrit)) %>%
            arrange(visitNumber) %>%
            head(1) %>%
            pull(hematocrit)
        if (length(x) == 0) {
            return(NA)
        } else {
            return(x)
        }
    })) %>%
    mutate(firstVisitWithRejection = map_int(data, function(x) {
        x <- x %>%
            mutate(anyRejection = pmap_lgl(list(rejection), function(r) {
                return(any(c(r)))
            }))
        if (any(x$anyRejection)) {
            return(
                x %>%
                    filter(rejection) %>%
                    arrange(visitNumber) %>%
                    head(1) %>%
                    pull(visitNumber)
            )
        } else {
            return(NA)
        }
    })) %>%
    mutate(allVisitsAfterFirstComplication = map(data, function(x) {
        x <- x %>%
            mutate(anyComplication = pmap_lgl(list(hospitalization, rejection, changeImmunosuppRegimen), function(a, b, c) {
                return(any(c(a, b, c)))
            }))
        if (any(x$anyComplication)) {
            # return(
            #     x %>%
            #         filter(anyComplication) %>%
            #         arrange(visitNumber) %>%
            #         pull(visitNumber)
            # )
            firstCompVisit <- x %>%
                filter(anyComplication) %>%
                arrange(visitNumber) %>%
                head(1) %>%
                pull(visitNumber)
            return(
                x %>%
                    filter(visitNumber >= firstCompVisit) %>%
                    pull(visitNumber)
            )
        } else {
            return(NA)
        }
    })) %>%
    mutate(nrVisitsWithComplication = map_int(data, function(x) {
        x <- x %>%
            mutate(anyComplication = pmap_lgl(list(hospitalization, rejection, changeImmunosuppRegimen), function(a, b, c) {
                return(any(c(a, b, c)))
            }))
        if (any(x$anyComplication)) {
            return(
                x %>%
                    filter(anyComplication) %>%
                    arrange(visitNumber) %>%
                    pull(visitNumber) %>%
                    length()
            )
        } else {
            return(0)
        }
    })) %>%
    mutate(typeFirstComplication = map_chr(data, function(x) {
        x <- x %>%
            mutate(anyComplication = pmap_lgl(list(hospitalization, rejection, changeImmunosuppRegimen), function(a, b, c) {
                return(any(c(a, b, c)))
            }))
        if (any(x$anyComplication)) {
            vecSel <- c("hospitalization", "rejection", "changeImmunosuppRegimen")
            # return(
            x <- x %>%
                filter(anyComplication) %>%
                arrange(visitNumber) %>%
                head(1) %>%
                select(all_of(vecSel)) %>%
                .[, vecSel]
            # apply(., 1, function(x) str_c(vecSel[!is.na(x)], sep = ",", collapse = ","))
            # x <- colnames(x)[x[1,]]
            x <- colnames(x)[unlist(as.data.frame(x)[1, , drop = T])]
            if (length(x) == 1) {
                return(x)
            } else {
                return("multComplications")
            }
        } else {
            return("None")
        }
    })) %>%
    mutate(anyComplicationEver = !is.na(firstVisitWithComplication)) %>%
    select(patientID,
        data,
        anyComplicationEver,
        typeFirstComplication,
        allVisitsAfterFirstComplication,
        nrVisitsWithComplication,
        firstVisitWithRejection,
        firstAlbuminMeasurement,
        firstHematocritMeasurement) %>%
    left_join(outcomeInformation %>% select(patientID, visitNumber) %>% rename(visit = visitNumber), by = c('patientID')) %>%
    full_join(pcoa %>% select(PSN, visit, V1, V2) %>% mutate(visit = as.numeric(as.character(visit))), by = c('patientID' = "PSN", "visit")) %>%
    mutate(anyComplicationEver2 = pmap_chr(list(anyComplicationEver, allVisitsAfterFirstComplication, visit), function(anyComp, allVisAfterFirstComp, vis) {
        if (!anyComp) {
            return("None")
        } else {
            if (vis %in% allVisAfterFirstComp) {
                return("AfterFirstComp")
            } else {
                return("BeforeFirstComp")
            }
        }
    })) %>%
    left_join(outcomeInformation %>% select(
        patientID,
        visitNumber,
        age,
        weight,
        ageCategorical,
        anyABx,
        ABxSubClass,
        postTransplant,
        studyCenter,
        all_of(names(model_covariates))) %>%
        mutate(visit = factor(visitNumber, levels = 1:7)) %>% select(-visitNumber) %>% mutate(visit = as.numeric(as.character(visit))), by = c('patientID', 'visit')) %>%
    group_by(patientID, visit) %>%
    nest() %>%
    group_by(patientID) %>%
    nest() %>%
    mutate(data = map(data, function(x) {
        x <- x %>%
            mutate(visit = as.numeric(as.character(visit))) %>%
            arrange(visit) %>%
            mutate(visitMinusOne = c(NA, visit[1:(length(visit) - 1)]))
        return(x)
    })) %>%
    unnest() %>%
    unnest() %>%
    # mutate(visitMinusOne = visit - 1) %>%
    mutate(visit = factor(visit, levels = 1:7)) %>%
    mutate(visitMinusOne = factor(visitMinusOne, levels = 1:6)) %>%
    # I'm removing data nested data as I'm fiddling with the visit column and I don't want to overwrite things by mistake...
    select(-data)

# For reason that will become clear later on (essentially, in certain situations, I want to predict outcome at T from microbiome + metadata at T-1),
# I hear distinguish between outcomeInformation and clinicalMetadata

# So to avoid confusion, outcomeInforamtion really only contains information relating to outcomes.
# In this object I will never mess with the time
outcomeInformation <- outcomeInformation %>%
    ungroup() %>%
    mutate(anyComplication = pmap_lgl(list(hospitalization, rejection, changeImmunosuppRegimen), function(a, b, c) {
        return(any(c(a, b, c)))
    })) %>%
    select(patientID,
        visitNumber,
        anyComplication,
        hospitalization,
        rejection,
        changeImmunosuppRegimen,
        CD,
        CDbinary)
# cyp3a4star22,
# cyp3a5star3)

# Transform patientID into factor where ordering corresponds to a meaningful ordering
orderDFPatientID <- clinicalMetadata %>%
    rename(PSN = patientID) %>%
    group_by(PSN) %>%
    left_join(outcomeInformation %>% select(patientID, visitNumber, CD, rejection, changeImmunosuppRegimen, hospitalization) %>% mutate(visitNumber = as.factor(visitNumber)), by = c("PSN" = 'patientID', "visit" = 'visitNumber')) %>%
    summarize(v = case_when(
        any(rejection) ~ "patientHadRejection",
        any(changeImmunosuppRegimen) ~ "patientHadChangeImmunoSuppRegime",
        any(hospitalization) ~ "patientWasHospitalized",
        .default = "NoComplication")) %>%
    mutate(v = factor(v, levels = rev(c(
        'patientHadRejection',
        "patientHadChangeImmunoSuppRegime",
        "patientWasHospitalized",
        "NoComplication"
    )), ordered = TRUE)) %>%
    arrange(v)

(ggplot(data = clinicalMetadata, aes(x = V1, y = V2, color = anyComplicationEver)) +
    geom_point() +
    theme_classic()) %>%
    ggsave(filename = "/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/AnycomplicationEverOrdination.pdf", width = 6.5, height = 5)

v <- RColorBrewer::brewer.pal((length(unique(clinicalMetadata %>% pull(anyComplicationEver2))) - 1), "Set2")
clinicalMetadata$anyComplicationEver2 <- factor(clinicalMetadata$anyComplicationEver2, levels = c("BeforeFirstComp", "AfterFirstComp", "None"))
(ggplot(data = clinicalMetadata, aes(x = V1, y = V2, color = anyComplicationEver2)) +
    geom_point() +
    theme_classic() +
    scale_color_manual(values = c(v[c(1, 2)], "#808080"))) %>%
    ggsave(filename = "/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/AnycomplicationEverOrdination2.pdf", width = 6.5, height = 5)

v <- RColorBrewer::brewer.pal((length(unique(clinicalMetadata %>% pull(typeFirstComplication))) - 2), "Set2")
clinicalMetadata$typeFirstComplication <- factor(clinicalMetadata$typeFirstComplication, levels = c(unique(clinicalMetadata$typeFirstComplication)[!unique(clinicalMetadata$typeFirstComplication) %in% c("None", "multComplications")], "multComplications", "None"))
(ggplot(data = clinicalMetadata, aes(x = V1, y = V2, color = typeFirstComplication)) +
    geom_point() +
    theme_classic() +
    scale_color_manual(values = c(v, 'red', "#808080"))) %>%
    ggsave(filename = "/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/TypeFirstComplicationOrdination.pdf", width = 7.5, height = 5)

v <- RColorBrewer::brewer.pal(length(unique(orderDFPatientID$v)), "Set1")
# clinicalMetadata$typeFirstComplication <- factor(clinicalMetadata$typeFirstComplication, levels = c(unique(clinicalMetadata$typeFirstComplication)[!unique(clinicalMetadata$typeFirstComplication) %in% c("None", "multComplications")], "multComplications", "None"))
(ggplot(data = clinicalMetadata %>%
    left_join(orderDFPatientID %>%
        rename(complicationsOrdered = v, patientID = PSN) %>%
        mutate(complicationsOrdered = factor(complicationsOrdered, levels = rev(levels(complicationsOrdered)), ordered = TRUE)), by = 'patientID'), aes(x = V1, y = V2, color = complicationsOrdered)) +
    geom_point() +
    theme_classic() +
    scale_color_manual(values = c("red", v[2:length(v)]))) %>%
    ggsave(filename = "/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/OrderComplicationOrdination.pdf", width = 7.5, height = 5)



p1 <- (ggplot(data = clinicalMetadata %>%
    left_join(profiles %>% group_by(PSN, visit) %>% summarize(richness = sum(relAb > log10(pseudoCount))) %>%
        # For ordination, limit richness to 150 max (there is a single outlier sample...)
        mutate(richness = ifelse(richness > 150, 150, richness)) %>% mutate(visit = factor(visit, levels = 1:7)), by = c('patientID' = "PSN", 'visit'))
, aes(x = V1, y = V2, color = richness)) +
    geom_point() +
    theme_classic() +
    ggtitle("richness capped at 150") +
    NULL +
    xlab("Principle Coordinate 1") +
    ylab("Principle Coordinate 2"))

# scale_color_manual(values = c("red", v[2:length(v)])))  %>%
ggsave(plot = p1, filename = "/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/RichnessOrdination.pdf", width = 7.5, height = 5)

tmp <- clinicalMetadata %>%
    left_join(profiles %>% group_by(PSN, visit) %>% summarize(richness = sum(relAb > log10(pseudoCount))) %>%
        # For ordination, limit richness to 150 max (there is a single outlier sample...)%>%
        mutate(visit = factor(visit, levels = 1:7)), by = c('patientID' = "PSN", 'visit')) %>%
    select(patientID, visit, V1, V2, richness) %>%
    left_join(orderDFPatientID %>%
        rename(complicationsOrdered = v, patientID = PSN) %>%
        mutate(complicationsOrdered = factor(complicationsOrdered, levels = rev(levels(complicationsOrdered)), ordered = TRUE)), by = 'patientID') %>%
    mutate(complicationsOrdered = as.character(complicationsOrdered)) %>%
    mutate(complicationsOrdered = ifelse(complicationsOrdered == "patientHadRejection", "patientHadRejection", "Other Patients")) %>%
    mutate(patientI = ifelse(complicationsOrdered != "Other Patients", patientID, NA)) %>%
    mutate(patientID = ifelse(complicationsOrdered != "Other Patients", patientID, "Other Patients")) %>%
    filter(!is.na(richness))
tmp$patientID <- factor(tmp$patientID, levels = c(unique(tmp$patientID)[unique(tmp$patientID) != "Other Patients"], "Other Patients"))
v <- RColorBrewer::brewer.pal(length(unique(tmp$patientID)) - 1, "Set2")
tmp <- tmp %>%
    mutate(v = as.numeric(as.character(visit))) %>%
    arrange(visit)
p2 <- (ggplot(data = tmp,
    aes(x = V1, y = V2, color = patientID, alpha = complicationsOrdered)) +
    geom_point() +
    geom_path(data = tmp %>% filter(!is.na(patientI)), aes(group = patientI)) +
    geom_point(data = tmp %>%
        filter(!is.na(patientI)) %>% inner_join(clinicalMetadata %>% ungroup() %>% filter(!is.na(firstVisitWithRejection)) %>% select(patientID, visit, firstVisitWithRejection) %>% distinct() %>% filter(visit == (firstVisitWithRejection - 1)), by = c('patientID', "visit")), color = 'red', size = 5, shape = 1) +
    # geom_point(data = tmp %>% filter(complicationsOrdered != "Other Patients"), inherit.aes = FALSE, aes(x = visit, y = richness, color = patientID)) +
    theme_classic() +
    geom_text(data = tmp %>% filter(!is.na(patientI)), aes(label = visit), size = 3, nudge_y = -0.25) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_alpha_manual(values = c('patientHadRejection' = 1, "Other Patients" = 0.1)) +
    NULL +
    scale_color_manual(values = c(v, "darkgrey")) +
    xlab("Principle Coordinate 1") +
    ylab("Principle Coordinate 2")
)
ggsave(plot = p1 / p2, filename = "/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/RejectionOrdination.pdf", width = 7.5, height = 8)


tmp <- clinicalMetadata %>%
    left_join(profiles %>% group_by(PSN, visit) %>% summarize(richness = sum(relAb > log10(pseudoCount))) %>%
        # For ordination, limit richness to 150 max (there is a single outlier sample...)%>%
        mutate(visit = factor(visit, levels = 1:7)), by = c('patientID' = "PSN", 'visit')) %>%
    select(patientID, visit, richness) %>%
    left_join(orderDFPatientID %>%
        rename(complicationsOrdered = v, patientID = PSN) %>%
        mutate(complicationsOrdered = factor(complicationsOrdered, levels = rev(levels(complicationsOrdered)), ordered = TRUE)), by = 'patientID') %>%
    mutate(complicationsOrdered = as.character(complicationsOrdered)) %>%
    mutate(complicationsOrdered = ifelse(complicationsOrdered == "patientHadRejection", "patientHadRejection", "Other Patients")) %>%
    mutate(patientI = patientID) %>%
    mutate(patientID = ifelse(complicationsOrdered != "Other Patients", patientID, "Other Patients")) %>%
    filter(!is.na(richness))
tmp$patientID <- factor(tmp$patientID, levels = c(unique(tmp$patientID)[unique(tmp$patientID) != "Other Patients"], "Other Patients"))
v <- RColorBrewer::brewer.pal(length(unique(tmp$patientID)) - 1, "Set2")
(ggplot(data = tmp,
    aes(x = visit, y = richness, color = patientID, group = patientI, alpha = complicationsOrdered)) +
    # geom_point() +
    # geom_boxplot() +
    geom_line() +
    geom_point(data = tmp %>% filter(complicationsOrdered != "Other Patients"), inherit.aes = FALSE, aes(x = visit, y = richness, color = patientID)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_alpha_manual(values = c('patientHadRejection' = 1, "Other Patients" = 0.4)) +
    NULL +
    scale_color_manual(values = c(v, "darkgrey"))) %>%
    ggsave(filename = "/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/RichnessOverTimeVsRejection_v3.pdf", width = 6.0, height = 3.5)

tmp <- clinicalMetadata %>%
    left_join(profiles %>% group_by(PSN, visit) %>% summarize(richness = sum(relAb > log10(pseudoCount))) %>%
        # For ordination, limit richness to 150 max (there is a single outlier sample...)%>%
        mutate(visit = factor(visit, levels = 1:7)), by = c('patientID' = "PSN", 'visit')) %>%
    select(patientID, visit, richness) %>%
    left_join(orderDFPatientID %>%
        rename(complicationsOrdered = v, patientID = PSN) %>%
        mutate(complicationsOrdered = factor(complicationsOrdered, levels = rev(levels(complicationsOrdered)), ordered = TRUE)), by = 'patientID') %>%
    mutate(complicationsOrdered = as.character(complicationsOrdered)) %>%
    mutate(complicationsOrdered = ifelse(complicationsOrdered == "patientWasHospitalized", "patientWasHospitalized", "Other Patients")) %>%
    mutate(patientI = patientID) %>%
    mutate(patientID = ifelse(complicationsOrdered != "Other Patients", patientID, "Other Patients")) %>%
    filter(!is.na(richness))
tmp$patientID <- factor(tmp$patientID, levels = c(unique(tmp$patientID)[unique(tmp$patientID) != "Other Patients"], "Other Patients"))
v <- RColorBrewer::brewer.pal(length(unique(tmp$patientID)) - 1, "Set2")
(ggplot(data = tmp,
    aes(x = visit, y = richness, color = patientID, group = patientI, alpha = complicationsOrdered)) +
    # geom_point() +
    # geom_boxplot() +
    geom_line() +
    geom_point(data = tmp %>% filter(complicationsOrdered != "Other Patients"), inherit.aes = FALSE, aes(x = visit, y = richness, color = patientID)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_alpha_manual(values = c('patientWasHospitalized' = 1, "Other Patients" = 0.4)) +
    NULL +
    scale_color_manual(values = c(v, "darkgrey"))) %>%
    ggsave(filename = "/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/RichnessOverTimeVsHospitalization_v3.pdf", width = 6.0, height = 3.5)

tmp <- clinicalMetadata %>%
    left_join(profiles %>% group_by(PSN, visit) %>% summarize(richness = sum(relAb > log10(pseudoCount))) %>%
        # For ordination, limit richness to 150 max (there is a single outlier sample...)%>%
        mutate(visit = factor(visit, levels = 1:7)), by = c('patientID' = "PSN", 'visit')) %>%
    select(patientID, visit, richness) %>%
    left_join(orderDFPatientID %>%
        rename(complicationsOrdered = v, patientID = PSN) %>%
        mutate(complicationsOrdered = factor(complicationsOrdered, levels = rev(levels(complicationsOrdered)), ordered = TRUE)), by = 'patientID') %>%
    mutate(complicationsOrdered = as.character(complicationsOrdered)) %>%
    mutate(complicationsOrdered = ifelse(complicationsOrdered == "patientHadChangeImmunoSuppRegime", "patientHadChangeImmunoSuppRegime", "Other Patients")) %>%
    mutate(patientI = patientID) %>%
    mutate(patientID = ifelse(complicationsOrdered != "Other Patients", patientID, "Other Patients")) %>%
    filter(!is.na(richness))
tmp$patientID <- factor(tmp$patientID, levels = c(unique(tmp$patientID)[unique(tmp$patientID) != "Other Patients"], "Other Patients"))
# v <- RColorBrewer::brewer.pal(length(unique(tmp$patientID)) - 1, "Set2")
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}
v <- gg_color_hue(length(unique(tmp$patientID)) - 1)
(ggplot(data = tmp,
    aes(x = visit, y = richness, color = patientID, group = patientI, alpha = complicationsOrdered)) +
    # geom_point() +
    # geom_boxplot() +
    geom_line() +
    geom_point(data = tmp %>% filter(complicationsOrdered != "Other Patients"), inherit.aes = FALSE, aes(x = visit, y = richness, color = patientID)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_alpha_manual(values = c('patientHadChangeImmunoSuppRegime' = 1, "Other Patients" = 0.4)) +
    NULL +
    scale_color_manual(values = c(v, "darkgrey"))) %>%
    ggsave(filename = "/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/RichnessOverTimeVsChangeImmunoSupp_v3.pdf", width = 6.0, height = 3.5)


tmp <- clinicalMetadata %>%
    left_join(profiles %>% group_by(PSN, visit) %>% summarize(richness = sum(relAb > log10(pseudoCount))) %>%
        # For ordination, limit richness to 150 max (there is a single outlier sample...)%>%
        mutate(visit = factor(visit, levels = 1:7)), by = c('patientID' = "PSN", 'visit')) %>%
    select(patientID, visit, richness) %>%
    # left_join(orderDFPatientID %>%
    #     rename(complicationsOrdered = v, patientID = PSN) %>%
    #     mutate(complicationsOrdered = factor(complicationsOrdered, levels = rev(levels(complicationsOrdered)), ordered = TRUE)), by = 'patientID') %>%
    # mutate(complicationsOrdered = as.character(complicationsOrdered)) %>%
    #     mutate(complicationsOrdered = ifelse(complicationsOrdered == "patientHadRejection", "patientHadRejection", "Other Patients")) %>%
    #     mutate(patientI = patientID) %>%
    #     mutate(patientID = ifelse(complicationsOrdered != "Other Patients", patientID, "Other Patients")) %>%
    filter(!is.na(richness))
# tmp$patientID <- factor(tmp$patientID, levels = c(unique(tmp$patientID)[unique(tmp$patientID)!="Other Patients"], "Other Patients"))
# v <- RColorBrewer::brewer.pal(length(unique(tmp$patientID)) - 1, "Set2")
(ggplot(data = tmp,
    aes(x = visit, y = richness)) +
    geom_jitter(height = 0, alpha = 0.3, width = 0.2) +
    geom_line(data = tmp %>% group_by(visit) %>% summarize(richness = mean(richness)), aes(group = 1)) +
    # geom_boxplot() +
    # geom_line() +
    # geom_point(data = tmp %>% filter(complicationsOrdered != "Other Patients"), inherit.aes = FALSE, aes(x = visit, y = richness, color = patientID)) +
    theme_classic() +
    NULL
# theme(axis.text.x = element_text(angle =45, hjust = 1)) +
# scale_alpha_manual(values = c('patientHadRejection' = 1, "Other Patients" = 0.4)) +
# NULL +
# scale_color_manual(values = c(v, "darkgrey"))
) %T>%
    ggsave(filename = "/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/RichnessOverTimePlusAverageRichness_v1.pdf", width = 6.0, height = 3.5) %>%
    ggsave(filename = "/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/RichnessOverTimePlusAverageRichness_v1.png", width = 6.0, height = 3.5)

pcoa_plot <- ggplot() +
    geom_point(data = pcoa %>%
        left_join(clinicalMetadata %>%
            select(patientID, studyCenter) %>% distinct(), by = c("PSN" = 'patientID')) %>%
        mutate(studyCenter = case_when(
            studyCenter == "k" ~ "Kinderklinik",
            studyCenter == "m" ~ "Muenster",
            studyCenter == "n" ~ "Nierenzentrum",
        )), aes(x = V1, y = V2, color = studyCenter)) +
    theme_classic() +
    xlab("PCo 1") +
    ylab("PCo 2") +
    NULL

ggsave(plot = pcoa_plot, filename = "plots/KLGPG_221206/pcoa_visit_study_center.pdf", width = 6.25, height = 5)
ggsave(plot = pcoa_plot, filename = "plots/KLGPG_221206/pcoa_visit_study_center.png", width = 6.25, height = 5)


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
    theme_classic() +
    xlab("PCo 1") +
    ylab("PCo 2") +
    NULL

ggsave(plot = pcoa_plot, filename = "plots/KLGPG_221206/pcoa_visit_age.pdf", width = 6.25, height = 5)
ggsave(plot = pcoa_plot, filename = "plots/KLGPG_221206/pcoa_visit_age.png", width = 6.25, height = 5)


tmp <- pairwiseDistances %>%
    as.matrix()
tmp[upper.tri(tmp)] <- NA
tmp <- tmp %>%
    as.data.frame() %>%
    rownames_to_column('sampleID') %>%
    pivot_longer(-sampleID) %>%
    rename(sampleID_1 = sampleID, sampleID_2 = name, distance = value) %>%
    filter(!is.na(distance)) %>%
    left_join(meta %>%
        select(sampleID, PSN, visit), by = c("sampleID_1" = "sampleID")) %>%
    left_join(meta %>%
        select(sampleID, PSN, visit), by = c("sampleID_2" = "sampleID"), suffix = c("_1", "_2")) %>%
    mutate(visit_1 = as.factor(visit_1)) %>%
    mutate(visit_2 = as.factor(visit_2)) %>%
    filter(PSN_1 == PSN_2) %>%
    rename(PSN = PSN_1) %>%
    select(-PSN_2) %>%
    filter(visit_1 != visit_2) %>%
    as_tibble()

tmp %>%
    group_by(PSN) %>%
    arrange(PSN, visit_2) %>%
    group_by(PSN) %>%
    filter(visit_2 == 1) %>%
    select(-visit_2)

get_family_level_barplot_for_all_samples_WIDE <- function(pObj, dataB, taxLevel = 'family', levelsToShow = NULL) {
    dataB <- dataB %>%
        # filter(PSN == indName) %>%
        mutate(relAb = (10^relAb) - pseudoCount) %>%
        mutate(visit = as.factor(visit)) %>%
        mutate(taxa = .data[[taxLevel]]) %>%
        mutate(taxa = as.character(taxa)) %>%
        mutate(taxa = ifelse(taxa %in% levelsToShow, taxa, "other")) %>%
        group_by(taxa, PSN, visit) %>%
        summarize(relAb = sum(relAb))
    dataC <- dataB %>%
        ungroup() %>%
        group_by(PSN, visit) %>%
        summarize(relAb = 1 - sum(relAb)) %>%
        mutate(taxa = "unclassified")
    dataB <- rbind(dataB, dataC) %>%
        mutate(taxa = factor(as.character(taxa), levels = c(levelsToShow, "other", 'unclassified')))

    dataB <- dataB %>%
        # left_join(
        # clinicalMetadata %>% select(patientID, visit, anyComplicationEver2, anyABx, ABxSubClass, hospitalization, rejection, changeImmunosuppRegimen),
        left_join(
            outcomeInformation %>%
                select(patientID, visitNumber, rejection, changeImmunosuppRegimen, hospitalization) %>%
                mutate(visitNumber = as.factor(visitNumber)),
            by = c('PSN' = 'patientID', "visit" = 'visitNumber'))

    dataB$PSN <- factor(dataB$PSN, levels = orderDFPatientID$PSN)

    dataB$visit <- factor(dataB$visit, levels = 1:7)

    dataB$visit <- factor(dataB$visit, levels = 1:7)

    pObj <- pObj +
        geom_bar(data = dataB,
            aes(x = PSN, y = relAb, fill = taxa), position = 'stack', stat = 'identity') +
        geom_text(data = outcomeInformation %>%
            rename(PSN = patientID, visit = visitNumber) %>%
            # select(PSN, visit, rejection) %>%
            mutate(rejection = ifelse(rejection, "R", "")) %>%
            mutate(PSN = factor(PSN, levels = orderDFPatientID$PSN)),
        aes(x = PSN, y = 1.525, label = rejection), size = 2.5) +
        geom_text(data = outcomeInformation %>%
            rename(PSN = patientID, visit = visitNumber) %>%
            # select(PSN, visit, changeImmunosuppRegimen) %>%
            mutate(changeImmunosuppRegimen = ifelse(changeImmunosuppRegimen, "C", "")) %>%
            mutate(PSN = factor(PSN, levels = orderDFPatientID$PSN)),
        aes(x = PSN, y = 1.375, label = changeImmunosuppRegimen), size = 2.5) +
        geom_text(data = outcomeInformation %>%
            rename(PSN = patientID, visit = visitNumber) %>%
            # select(PSN, visit, changeImmunosuppRegimen) %>%
            mutate(hospitalization = ifelse(hospitalization, "H", "")) %>%
            mutate(PSN = factor(PSN, levels = orderDFPatientID$PSN)),
        aes(x = PSN, y = 1.225, label = hospitalization), size = 2.5) +
        # geom_point(data = dataB %>% ungroup() %>% filter(anyABx) %>% select(PSN, ABxSubClass) %>% distinct()  ,
        # aes(x = PSN, y = 1.1, color = ABxSubClass)) +
        theme_classic() +
        # facet_grid(anyComplicationEver2 ~ visit)
        facet_grid(visit ~ ., scales = "fixed") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
        scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1.725)) +
        # ylim(c(0, 1.725)) +
        ylab("relative bacterial abundance")

    return(pObj)
}



# LONG PLOT IS OUTDATED
# for (taxLevel in c(
#     #"phylum",
#     #"class",
#     #"order",
#     "family"
#     #"genus"
#     )
#     ) {
#         set.seed(2)
#         #taxL <- "family"
#         numTaxa <- 10
#         levelsToShow <- profiles %>%
#             mutate(relAb = 10^relAb - pseudoCount) %>%
#             group_by(sampleID, .data[[taxL]]) %>%
#             summarize(relAb = sum(relAb)) %>%
#             group_by(.data[[taxL]]) %>%
#             summarize(m = mean(relAb)) %>%
#             arrange(desc(m)) %>%
#             head(numTaxa) %>%
#             pull(.data[[taxL]])
#         getPalette <- colorRampPalette(brewer.pal(9, "Paired"))
#         colors <- sample(getPalette(numTaxa))
#         colors <- c(colors, "#808080", "#D3D3D3")

#         ggsave(
#             plot = get_family_level_barplot_for_all_samples_LONG(ggplot(),
#             profiles %>%
#                 #filter(visit != 3) %>%
#                 mutate(visit = as.numeric(as.character(visit))),
#             taxL,
#             levelsToShow = levelsToShow) +
#             scale_fill_manual(values = colors),

#             filename = str_c("plots/KLGPG_221206/taxa_barplots_visit_v3_visit_3_removed_treatment_dichotomized_all__LONG_VERSION_", taxL, ".pdf", collapse = ""),
#                 width = 6,
#                 height = 28
#                 )

#     }

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
        filename = str_c("plots/KLGPG_221206/taxa_barplots_visit_v3_visit_3_removed_treatment_dichotomized_all__WIDE_VERSION_", taxL, ".pdf", collapse = ""),
        width = 8,
        height = 10
    )

}

#####################################
# LMMs predicting complication status
#####################################

library(dplyr)
library(lmerTest)
library(lme4)

lmms <- list()
lmmsAdjusted <- list()

# TODO: thse dont have weight
formulaStrings <- c(
    # "lmmsAdjusted" = function(outcomeString) return(paste(outcomeString, "~ log10(relAb + pseudoCount) + anyABx + ageCategorical + albumin + hematocrit + sex + (1 | PSN)")),
    "lmmsAdjusted" = function(outcomeString) return(paste(outcomeString, "~ log10(relAb + pseudoCount) + anyABx + ageCategorical + (1 | PSN)")),
    "lmms" = function(outcomeString) return(paste(outcomeString, "~ log10(relAb + pseudoCount) + (1 | PSN)"))
)

# formulaStringAdjusted <- paste(outcomeString, "~ log10(relAb + pseudoCount) + anyABx + ageCategorical + (1 | PSN)")
# formulaString <- paste(outcomeString, "~ log10(relAb + pseudoCount) + (1 | PSN)")

for (taxLevel in c(
    # "phylum",
    # "class",
    # "order",
    # "family"
    "genus"
)
) {
    # for (visitType in c("visit", "visitMinusOne")) {
    for (visitType in c("visit", "visitMinusOne")) {
        importantTaxa <- profiles %>%
            mutate(relAb = (10^relAb) - pseudoCount) %>%
            mutate(taxa = .data[[taxLevel]]) %>%
            mutate(taxa = as.character(taxa)) %>%
            group_by(taxa, PSN, visit) %>%
            summarize(relAb = sum(relAb)) %>%
            group_by(taxa) %>%
            summarize(m = mean(relAb > pseudoCount) > 0.2 && any(relAb > 0.01)) %>%
            filter(m) %>%
            select(taxa) %>%
            filter(taxa %in% c("Roseburia", "Coprococcus", "Anaerostipes", "Enterococcus"))

        genusProfiles <- profiles %>%
            mutate(relAb = (10^relAb) - pseudoCount) %>%
            mutate(taxa = .data[[taxLevel]]) %>%
            mutate(taxa = as.character(taxa)) %>%
            ## CRUCIAL FOR HOW MODEL IS FIT.
            # DISTINGUISHES BETWEEN STOOL(T) <-> META(T) and STOOL(T-1) <-> META(T)
            mutate(visit = .data[[visitType]]) %>%
            group_by(phylum, taxa, PSN, visit) %>%
            summarize(relAb = sum(relAb))
        tmp <- genusProfiles %>%
            inner_join(importantTaxa)
        # Also add richness
        tmp <- rbind(
            tmp,
            genusProfiles %>% group_by(PSN, visit) %>% summarize(relAb = sum(relAb >= 1E-4)) %>% mutate(taxa = 'richness'),
            genusProfiles %>% group_by(PSN, visit) %>% filter(phylum == "Proteobacteria") %>% summarize(relAb = sum(relAb >= 1E-4)) %>% mutate(taxa = 'richnessProteobacteria'),
            genusProfiles %>% group_by(PSN, visit) %>% filter(phylum == "Firmicutes") %>% summarize(relAb = sum(relAb >= 1E-4)) %>% mutate(taxa = 'richnessFirmicutes'),
            genusProfiles %>% group_by(PSN, visit) %>% filter(phylum == "Bacteroidetes") %>% summarize(relAb = sum(relAb >= 1E-4)) %>% mutate(taxa = 'richnessBacteroidetes'),
            genusProfiles %>% group_by(PSN, visit) %>% filter(phylum == "Actinobacteria") %>% summarize(relAb = sum(relAb >= 1E-4)) %>% mutate(taxa = 'richnessActinobacteria')) %>%
            inner_join(
                # Join the CLINICAL METADATA from clinicalMetadata via the appropriate visit type (T or T-1)
                clinicalMetadata %>%
                    ungroup() %>%
                    select(
                        patientID,
                        .data[[visitType]],
                        age,
                        ageCategorical,
                        anyABx,
                        albumin,
                        hematocrit,
                        sex
                    ) %>%
                    # Next I overwrite visit (afterwards used to join complication metadata) with T or T-1
                    mutate(visit = .data[[visitType]]) %>%
                    distinct() %>%
                    mutate(visit = as.numeric(as.character(visit))),
                by = c('PSN' = 'patientID', "visit")) %>%
            # join complication metadata (will always remain at T)
            left_join(outcomeInformation, by = c("PSN" = 'patientID', "visit" = 'visitNumber')) %>%
            mutate(visit = factor(visit, levels = 1:7))
        # for (outcomeString in c("anyComplication", "rejection", "changeImmunosuppRegimen", "hospitalization", "CD")) {
        for (outcomeString in c("anyComplication", "rejection", "hospitalization", "CD", "CDbinary")) {
            # for (outcomeString in c("anyComplication", "rejection", "hospitalization", "CD")) {
            print(outcomeString)
            for (taxon in unique(tmp$taxa)) {

                tmp2 <- tmp %>%
                    filter(taxa == taxon)

                if (outcomeString != "CDbinary") {
                    lmmAdjusted <- lmer(formula = as.formula(formulaStrings[['lmmsAdjusted']](outcomeString)), data = tmp2)
                } else {
                    lmmAdjusted <- tryCatch(
                        # Sometimes the glmer function fails o converge, in which case nonsensically small p-values are returned.
                        # In order to give the model a fair chance to converge, I'm increasing the parameters here.
                        glmer(formula = as.formula(formulaStrings[['lmmsAdjusted']](outcomeString)), data = tmp2, family = 'binomial', control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000))),
                        error = function(e) e,
                        # Sometimes the glmer function fails o converge, in which case nonsensically small p-values are returned.
                        # I'm catching this here
                        warning = function(w) {
                            print(w$message)
                            if (str_detect(w$message, "failed to converge")) {
                                return(NA)
                            } else if (str_detect(w$message, "unable to evaluate scaled gradient")) {
                                return(NA)
                            } else if (str_detect(w$message, "Model is nearly unidentifiable: large eigenvalue rati")) {
                                return(NA)
                            } else {
                                print("This should never be reached.. but it's fine for now!")
                                return(NA)
                            }
                        }
                    )
                }

                lmmsAdjusted[[length(lmmsAdjusted) + 1]] <- lmmAdjusted
                names(lmmsAdjusted)[length(lmmsAdjusted)] <- str_c(outcomeString, taxon, visitType, sep = "___", collapse = "___")

                if (outcomeString != "CDbinary") {
                    lmm <- lmer(formula = as.formula(formulaStrings[['lmms']](outcomeString)), data = tmp2)
                } else {
                    lmm <- tryCatch(
                        glmer(formula = as.formula(formulaStrings[['lmms']](outcomeString)), data = tmp2, family = 'binomial', control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000))),
                        error = function(e) e,
                        warning = function(w) {
                            print(w$message)
                            if (str_detect(w$message, "failed to converge")) {
                                return(NA)
                            } else if (str_detect(w$message, "unable to evaluate scaled gradient")) {
                                return(NA)
                            } else if (str_detect(w$message, "Model is nearly unidentifiable: large eigenvalue rati")) {
                                return(NA)
                            } else {
                                print("This should never be reached.. but it's fine for now!")
                                return(NA)
                            }
                        })
                }
                lmms[[length(lmms) + 1]] <- lmm
                names(lmms)[length(lmms)] <- str_c(outcomeString, taxon, visitType, sep = "___", collapse = "___")
            }
        }
    }
}

for (lmmType in c("lmms", "lmmsAdjusted")) {
    for (compType in c('primary', "secondary")) {
        lmmsTibble <- tibble(outcomeTaxa = names(lmms), lmms = lmms, lmmsAdjusted = lmmsAdjusted) %>%
            mutate(visitType = str_split_fixed(outcomeTaxa, "___", n = 3)[, 3]) %>%
            mutate(taxa = str_split_fixed(outcomeTaxa, "___", n = 3)[, 2]) %>%
            mutate(outcome = str_split_fixed(outcomeTaxa, "___", n = 3)[, 1]) %>%
            # CHANGE THIS DEPENDING ON WHAT YOU WANT TO EVALUATE
            mutate(coefs = map(.data[[lmmType]], \(x) {
                if (is.na(x)) {
                    return(NA)
                } else {
                    coefficients(summary(x))
                }
            })) %>%
            mutate(pValTaxon = map2_dbl(coefs, outcome, \(x, y) {
                if (is.na(x)) {
                    return(NA)
                }
                if (y != 'CDbinary') {
                    return(x[str_detect(rownames(x), 'relAb'), 5])
                } else {
                    return(x[str_detect(rownames(x), 'relAb'), 4])
                }

            })) %>%
            mutate(effectSizeTaxon = map_dbl(coefs, \(x) {
                if (is.na(x)) {
                    return(NA)
                }
                return(x[str_detect(rownames(x), 'relAb'), 1])
            })) %>%
            group_by(outcome) %>%
            nest() %>%
            mutate(data = map(data, function(x) {
                x$pValTaxonAdjusted <- p.adjust(x$pValTaxon, method = "BH")
                return(x)
            })) %>%
            unnest() %>%
            inner_join(data.frame(
                outcome = c(
                    "anyComplication",
                    "hospitalization",
                    "rejection",
                    "changeImmunosuppRegimen",
                    "CD",
                    "CDbinary"
                ),
                complicationType = c(
                    "secondary",
                    "secondary",
                    "secondary",
                    "secondary",
                    "primary",
                    "primary")) %>%
                filter(complicationType == compType))

        lmmsTibble <- lmmsTibble %>%
            # filter(pValTaxonAdjusted < 0.2) %>%
            arrange(pValTaxonAdjusted)

        lmmsTibble <- lmmsTibble %>%
            select(taxa, outcome, pValTaxon, pValTaxonAdjusted, effectSizeTaxon, visitType) %>%
            filter(!str_detect(taxa, '\\[')) %>%
            filter(pValTaxon < 0.05) %>%
            # group_by(taxa) %>%
            # nest() %>%
            # filter(map_lgl(data, \(x) {
            #     any(x$pValTaxon < 0.05)
            # })) %>%
            # unnest() %>%

            mutate(pValTaxonSigned = ifelse(effectSizeTaxon > 0, pValTaxon, -pValTaxon)) %>%
            mutate(direction = ifelse(
                effectSizeTaxon > 0,
                "Bacterium more abundant in complication/when CD high",
                "Bacterium more abundant in complication/when CD low")) %>%
            mutate(featureGroup = ifelse(str_detect(taxa, 'richness'), "diversity", "taxon")) %>%
            mutate(fillColor = case_when(
                (pValTaxon < 0.05 & pValTaxonAdjusted < 0.1) ~ "p << 0.05",
                (pValTaxon < 0.05 & pValTaxonAdjusted >= 0.1) ~ "p < 0.05",
                (pValTaxon >= 0.05 & pValTaxonAdjusted >= 0.1) ~ "p >= 0.05",
                .default = NA
            )
            ) %>%
            mutate(fillColor = factor(fillColor, levels = c("p << 0.05", "p < 0.05", "p >= 0.05"))) %>%
            mutate(outcomeGroup = ifelse(outcome == "CD", "primary", "secondary")) %>%
            mutate(outcome = factor(outcome, levels = c("CD", unique(outcome)[unique(outcome) != "CD"]))) %>%
            left_join(
                fullTax %>%
                    ungroup() %>%
                    add_row(phylum = NA, genus = "richness"), by = c('taxa' = 'genus')) %>%
            arrange(phylum)
        lmmsTibble$taxa <- factor(lmmsTibble$taxa, levels = unique(lmmsTibble$taxa))

        lmmPhylumBar <- ggplot() +
            geom_tile(data = lmmsTibble, aes(x = 1, y = taxa, fill = phylum)) +
            theme_classic() +
            NULL +
            # facet_grid(featureGroup ~ ., scales = 'free', space = 'free') +
            theme(
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.x = element_blank()
            )

        numTaxa <- length(unique(lmmsTibble$taxa))
        stopifnot(class(lmmsTibble$taxa) == "factor")
        gridDf <- data.frame(taxa = levels(lmmsTibble$taxa), horizontalBars = rep(c("white", "gray"), 1000)[1:numTaxa])
        gridDf$taxa <- factor(gridDf$taxa, levels = levels(lmmsTibble$taxa))
        gridDf <- rbind(
            gridDf %>% mutate(visitType = 'visit'),
            gridDf %>% mutate(visitType = 'visitMinusOne'))
        #     mutate(featureGroup = ifelse(taxa %in% c("richness"), "diversity", "taxon")) %>%
        #     filter(taxa !='richness')

        lmmHeatmap <- ggplot() +
            geom_tile(data = lmmsTibble,
                aes(x = outcome, y = taxa, fill = fillColor, color = direction), width = 0.8, height = 0.8, linewidth = 1.25) +
            geom_rect(data = gridDf %>% filter(horizontalBars == "white"),
                aes(xmin = -Inf, xmax = +Inf, ymin = as.numeric(taxa) - 0.5, ymax = as.numeric(taxa) + 0.5), fill = "grey", alpha = 0.3) +
            # geom_tile(data = lmmsTibble,
            # aes(x = outcome, y = taxa, fill = fillColor, color = direction), width = 0.8, height = 0.8, linewidth = 1.25) +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            facet_grid(. ~ visitType, scales = 'free', space = 'free') +
            # scale_fill_gradient(low = "red", high = 'blue') +
            scale_fill_manual(values = c("black", "grey", "white")) +
            scale_color_manual(values = c("darkgreen", "darkred")) +
            ggtitle(formulaStrings[[lmmType]]('outcome')) +
            theme(
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title.y = element_blank()
            )
        if (compType == "primary") {
            ggsave(plot = lmmPhylumBar + lmmHeatmap + plot_layout(widths = c(1, 10), guides = 'collect'), filename = str_c("/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/LmmHeatmap_", lmmType, "__", compType, ".pdf"), width = 7.5, height = 7)
        } else if (compType == "secondary") {
            ggsave(plot = lmmPhylumBar + lmmHeatmap + plot_layout(widths = c(1, 10), guides = 'collect'), filename = str_c("/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/LmmHeatmap_", lmmType, "__", compType, ".pdf"), width = 8, height = 7)
        }
    }
}

visualize_taxon_vs_complication_lineplot <- function(
    data = NULL,
    taxon = NULL,
    complication = NULL,
    visitType = NULL,
    taxLevel = "genus"
    ) {
    genusProfiles <- profiles %>%
        mutate(relAb = (10^relAb) - pseudoCount) %>%
        mutate(taxa = .data[[taxLevel]]) %>%
        mutate(taxa = as.character(taxa)) %>%
        ## CRUCIAL FOR HOW MODEL IS FIT.
        # DISTINGUISHES BETWEEN STOOL(T) <-> META(T) and STOOL(T-1) <-> META(T)
        mutate(visit = .data[[visitType]]) %>%
        group_by(phylum, taxa, PSN, visit) %>%
        summarize(relAb = sum(relAb))
    tmp <- genusProfiles %>%
        inner_join(importantTaxa)
    # Also add richness
    tmp <- rbind(
        tmp,
        genusProfiles %>% group_by(PSN, visit) %>% summarize(relAb = sum(relAb >= 1E-4)) %>% mutate(taxa = 'richness'),
        genusProfiles %>% group_by(PSN, visit) %>% filter(phylum == "Proteobacteria") %>% summarize(relAb = sum(relAb >= 1E-4)) %>% mutate(taxa = 'richnessProteobacteria'),
        genusProfiles %>% group_by(PSN, visit) %>% filter(phylum == "Firmicutes") %>% summarize(relAb = sum(relAb >= 1E-4)) %>% mutate(taxa = 'richnessFirmicutes'),
        genusProfiles %>% group_by(PSN, visit) %>% filter(phylum == "Bacteroidetes") %>% summarize(relAb = sum(relAb >= 1E-4)) %>% mutate(taxa = 'richnessBacteroidetes'),
        genusProfiles %>% group_by(PSN, visit) %>% filter(phylum == "Actinobacteria") %>% summarize(relAb = sum(relAb >= 1E-4)) %>% mutate(taxa = 'richnessActinobacteria')) %>%
        inner_join(
            # Join the CLINICAL METADATA from clinicalMetadata via the appropriate visit type (T or T-1)
            clinicalMetadata %>%
                ungroup() %>%
                select(
                    patientID,
                    .data[[visitType]],
                    age,
                    ageCategorical,
                    anyABx
                ) %>%
                # Next I overwrite visit (afterwards used to join complication metadata) with T or T-1
                mutate(visit = .data[[visitType]]) %>%
                distinct() %>%
                mutate(visit = as.numeric(as.character(visit))),
            by = c('PSN' = 'patientID', "visit")) %>%
        # join complication metadata (will always remain at T)
        left_join(outcomeInformation, by = c("PSN" = 'patientID', "visit" = 'visitNumber')) %>%
        mutate(visit = factor(visit, levels = 1:7)) %>%
        filter(taxa == taxon) %>%
        mutate(`log10(relAb)` = log10(relAb + pseudoCount)) %>%
        group_by(PSN) %>%
        nest() %>%
        mutate(data = map(data, \(x) {
            x['complicationEver'] <- any(x[complication])
            return(x)
        })) %>%
        unnest()
    # return(tmp)
    tmp2 <- tmp %>%
        # filter(.data[[complication]]) %>%
        group_by(PSN) %>%
        nest() %>%
        mutate(data = map(data, \(x) {
            m <- min(as.numeric(as.character(x %>% filter(.data[[complication]]) %>% pull(visit))))
            x <- x %>%
                arrange(as.numeric(as.character(visit)))
            x <- rbind(x %>% filter(.data[[complication]]), x %>% filter(visit == m - 1))
            return(x)
        })) %>%
        unnest()
    print(tmp2)
    p <- ggplot() +
        geom_line(data = tmp, aes_string(x = "visit", y = "log10(relAb)", group = "PSN", color = "complicationEver"), show.legend = FALSE) +
        # geom_point(data = tmp %>% filter(rejection), aes_string(x = "visit", y = "log10(relAb)", color = complication)) +
        geom_line(data = tmp2, aes_string(x = "visit", y = "log10(relAb)", group = "PSN"), color = 'green') +
        theme_classic() +
        scale_x_discrete_prisma() +
        scale_color_manual(values = c(
            "FALSE" = "grey",
            "TRUE" = "red")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    return(p)
}

ggsave(
    plot = visualize_taxon_vs_complication_lineplot(taxon = "Roseburia", complication = 'CD', visitType = "visit", taxLevel = "genus"),
    filename = "/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/taxon_vs_complication_test.pdf", width = 5, height = 4
)

test

# tmp <- profiles %>%
#             mutate(relAb = (10^relAb) - pseudoCount) %>%
#             mutate(taxa = .data[[taxLevel]]) %>%
#             mutate(taxa = as.character(taxa)) %>%
#             #mutate(taxa = ifelse(taxa %in% levelsToShow, taxa, "other")) %>%
#             # mutate(taxa = factor(taxa, levels = c(levelsToShow, "other", 'unclassified'))) %>%
#             group_by(taxa, PSN, visit) %>%
#             summarize(relAb = sum(relAb)) %>%
#             inner_join(importantTaxa) %>%
#             left_join(
#                 # clinicalMetadata %>% select(patientID, visit, anyComplicationEver2, anyABx),
#                 clinicalMetadata %>% ungroup() %>% select(patientID, data) %>%
#                     unnest() %>%
#                     ungroup() %>%
#                     select(patientID, visitNumber, anyComplication, anyABx) %>%
#                     rename(visit = visitNumber) %>%
#                     distinct(),
#             by = c('PSN' = 'patientID', "visit")) %>%
#                 mutate(visit = factor(visit, levels = 1:7))

# compareTaxAssocsBoxplots <- function(taxon = NULL) {
#     tmp2 <- tmp %>%
#         filter(taxa == taxon) %>%
#         mutate(anyABx = ifelse(anyABx, "ABx", "no ABx"))
#     p <- ggplot() +
#         geom_boxplot(data = tmp2, aes(x = anyComplication, y = log10(relAb + pseudoCount)), outlier.colour = NA) +
#         geom_text(data =
#         tmp2 %>%
#         group_by(anyComplication, anyABx) %>%
#         tally(), aes(x = anyComplication, y = -(abs(max(log10(tmp2$relAb + pseudoCount))) * 0.95), label = n)) +
#         geom_jitter(data = tmp2, aes(x = anyComplication, y = log10(relAb + pseudoCount)), alpha = 0.5, width = 0.2) +
#         theme_classic() +
#         facet_grid(. ~ anyABx)
#     p <- p +
#         ggtitle(str_c("LMM BH-adj.p-val for ", taxon, ': ', round(lmmsTibble %>% filter(taxa == taxon) %>% pull(pValTaxonAdjusted), 4)))
#     return(p)
# }


# get_quantile_data <- function(i, cont = "relAb", ...) {
#     args <- c(...)
#     quantileData <- list()
#     #SSS <- inputData
#     # inputData <- inputData %>%
#     #  filter(type == "reference (CV)")
#     # medians <- i %>%
#     # group_by(...) %>%
#     # summarize(median = median(relAb + pc))
#     for (quantile in quantiles) {
#         tmp <- i %>%
#         #select(cyl, wt) %>%
#         group_by(across(all_of(args))) %>%
#         #mutate(cyl) %>%
#         summarize(value_max = quantile(.data[[cont]]+1E-5, probs = quantile),
#                 value_min = quantile(.data[[cont]]+1E-5, probs = 1-quantile)) %>%
#         mutate(quantile = quantile)
#         quantileData[[length(quantileData) + 1]] <- tmp
#     }
#     var <- args[1]
#     quantileData <- do.call('rbind', quantileData) %>%
#     #left_join(medians, by = c("genus", "label")) %>%
#     mutate({{var}} := as.factor(.data[[var]])) %>%
#     arrange(desc(quantile)) %>%
#     mutate(Quantiles = map2_chr(quantile, .data[[args[length(args)]]], function(x, y) return(str_c((1-x) * 100, '% - ', x * 100, "% - ", y))))
#     return(quantileData)
# }

get_quantile_plot <- function(inputData, axisColumn, labelColumn, valueColumn, plotGroup = "postTransplant", expectedNumLevels = 3, xlab = 'Genera', ylab = "Relative Abundances (log10)") {

    colnames(inputData)[colnames(inputData) == axisColumn] <- "genus"
    colnames(inputData)[colnames(inputData) == labelColumn] <- "label"
    colnames(inputData)[colnames(inputData) == valueColumn] <- "relAb"

    if (is.na(plotGroup)) {
        inputData$plotGroup <- 1
    } else {
        inputData$plotGroup <- inputData[[plotGroup]]
    }

    # print(head(inputData))
    # print(dim(inputData))
    # groupColors <- c("red", "blue")
    # quantiles <- c(0.5, 0.7, 0.9, 0.95)
    quantileData <- list()
    # SSS <- inputData
    # inputData <- inputData %>%
    #  filter(type == "reference (CV)")
    # medians <- inputData %>%
    # group_by(genus, label, plotGroup) %>%
    # summarize(median = median(relAb))
    for (quantile in quantiles) {
        tmp <- inputData %>%
            # select(cyl, wt) %>%
            group_by(genus, label, plotGroup) %>%
            # mutate(cyl) %>%
            summarize(value_max = quantile(relAb, probs = quantile),
                value_min = quantile(relAb, probs = 1 - quantile)) %>%
            mutate(quantile = quantile)
        quantileData[[length(quantileData) + 1]] <- tmp
    }

    quantileData <- do.call('rbind', quantileData) %>%
        # left_join(medians %>% select(-plotGroup), by = c("genus", "label")) %>%
        mutate(genus = as.factor(genus)) %>%
        arrange(desc(quantile)) %>%
        # mutate(Quantiles = map2_chr(quantile, label, function(x, y) return(str_c((1-x) * 100, '% - ', x * 100, "% - ", y))))
        mutate(Quantiles = map2_chr(quantile, label, function(x, y) {
            return(str_c((x * 100), "% - ", y))
        }))

    # print(head(quantileData))


    print(head(quantileData))
    # return(quantileData)
    l <- quantileData %>%
        ungroup() %>%
        select(Quantiles, quantile, label) %>%
        distinct() %>%
        group_by(label) %>%
        nest() %>%
        mutate(data = map(data, function(x) return(x %>% arrange(quantile)))) %>%
        unnest() %>%
        # arrange(quantile) %>%
        ungroup() %>%
        pull(Quantiles)
    l <- c(l[!str_detect(l, "no ")], l[str_detect(l, "no ")])
    print(l)
    names(colorVec) <- l

    quantileData <- quantileData %>%
        mutate(Quantiles = factor(Quantiles, levels = (l))) %>%
        arrange(desc(quantile))


    labelMap <- levels(quantileData$genus)
    names(labelMap) <- 1:length(labelMap)

    # print(levels(quantileData$Quantiles))
    # print(length(unique(quantileData$label)))
    # print(expectedNumLevels)
    if (length(unique(quantileData$label)) != expectedNumLevels) {
        print(unique(quantileData$label))
        asdaddads
    }

    groupLevels <- unique(quantileData$label)
    # groupLevels <- labelLevels

    colorVec <- c(colorRampPalette(c(groupColors[1], "white"))(length(quantiles)),
        colorRampPalette(c(groupColors[2], "white"))(length(quantiles)))
    names(colorVec) <- l
    print(head(quantileData))
    # print(dim(quantileData %>% ungroup() %>% group_by()))
    print(inputData %>%
        # select(cyl, wt) %>%
        group_by(genus, label, plotGroup) %>%
        tally())

    p <- ggplot(data = quantileData %>%
        mutate(Quantiles = factor(Quantiles, levels = l))) +
        geom_rect(data = quantileData %>%
            filter(label == groupLevels[1]) %>%
            mutate(Quantiles = factor(Quantiles, levels = l)),
        aes(xmin = as.integer(genus) - 0.1 - 0.15, xmax = as.integer(genus) + 0.1 - 0.15, ymin = value_max, ymax = value_min, fill = Quantiles), color = 'black') +
        geom_text(data = inputData %>%
            # select(cyl, wt) %>%
            group_by(genus, label, plotGroup) %>%
            filter(label == groupLevels[1]) %>%
            mutate(genus = factor(genus, levels = levels(quantileData$genus))) %>%
            tally(),
        aes(x = as.integer(genus) - 0.15, y = 0, label = n), color = 'black') +
        geom_rect(data = quantileData %>%
            filter(label == groupLevels[2]) %>%
            mutate(Quantiles = factor(Quantiles, levels = l)),
        aes(xmin = as.integer(genus) - 0.1 + 0.15, xmax = as.integer(genus) + 0.1 + 0.15, ymin = value_max, ymax = value_min, fill = Quantiles), color = 'black') +
        geom_text(data = inputData %>%
            # select(cyl, wt) %>%
            group_by(genus, label, plotGroup) %>%
            filter(label == groupLevels[2]) %>%
            mutate(genus = factor(genus, levels = levels(quantileData$genus))) %>%
            tally(),
        aes(x = as.integer(genus) + 0.15, y = 0, label = n), color = 'black') +
        geom_point(data = quantileData %>%
            filter(label == groupLevels[1]) %>%
            filter(str_detect(Quantiles, "50")) %>%
            mutate(Quantiles = factor(Quantiles, levels = l)),
        aes(x = as.integer(genus) - 0.15, y = value_max), fill = 'darkgreen', size = 3, pch = 23, color = 'white') +
        geom_point(data = quantileData %>%
            filter(label == groupLevels[2]) %>%
            filter(str_detect(Quantiles, "50")) %>%
            mutate(Quantiles = factor(Quantiles, levels = l)),
        aes(x = as.integer(genus) + 0.15, y = value_max), fill = 'darkgreen', size = 3, pch = 23, color = 'white') +
        scale_x_continuous(breaks = as.integer(names(labelMap)), labels = labelMap) +
        # scale_fill_manual(breaks = c("Dark","DarkLight","Medium","LightDark","Light"),
        #                values=c("red", "orange","yellow","cadetblue2","dodgerblue"))
        scale_fill_manual(values = colorVec[!str_detect(names(colorVec), '50')], breaks = names(colorVec[!str_detect(names(colorVec), '50')])) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        # coord_flip() +
        # facet_wrap(~plotGroup, nrow = 2, scales = "free") +
        xlab(xlab) +
        ylab(ylab) +
        NULL
    return(p)
}


quantiles <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
groupColors <- c("#852020", '#21520e', "#090938")
colorVec <- c(colorRampPalette(c(groupColors[1], "white"))(length(quantiles)),
    colorRampPalette(c(groupColors[3], "white"))(length(quantiles)),
    colorRampPalette(c(groupColors[2], "white"))(length(quantiles)))
colorVecGreen <- colorRampPalette(c(groupColors[2], "white"))(length(quantiles))

############################################################################################
## IMPORTANT: Choose visitType to make sure you're correlating/looking at the right thing
############################################################################################

visitType <- "visit"
# visitType <- "visitMinusOne"
genusProfiles <- profiles %>%
    mutate(relAb = (10^relAb) - pseudoCount) %>%
    mutate(taxa = .data[[taxLevel]]) %>%
    mutate(taxa = as.character(taxa)) %>%
    ## CRUCIAL FOR HOW MODEL IS FIT.
    # DISTINGUISHES BETWEEN STOOL(T) <-> META(T) and STOOL(T-1) <-> META(T)
    mutate(visit = .data[[visitType]]) %>%
    group_by(taxa, PSN, visit) %>%
    summarize(relAb = sum(relAb))
genusProfiles <- rbind(genusProfiles, genusProfiles %>% group_by(PSN, visit) %>% summarize(relAb = sum(relAb >= 1E-4)) %>% mutate(taxa = 'richness')) %>%
    inner_join(
        # Join the CLINICAL METADATA from clinicalMetadata via the appropriate visit type (T or T-1)
        clinicalMetadata %>%
            ungroup() %>%
            select(
                patientID,
                .data[[visitType]],
                # visit,
                age,
                anyABx
            ) %>%
            # Next I overwrite visit (afterwards used to join complication metadata) with T or T-1
            mutate(visit = .data[[visitType]]) %>%
            distinct() %>%
            mutate(visit = as.numeric(as.character(visit))),
        by = c('PSN' = 'patientID', "visit")) %>%
    # join complication metadata (will always remain at T)
    left_join(outcomeInformation, by = c("PSN" = 'patientID', "visit" = 'visitNumber')) %>%
    mutate(visit = factor(visit, levels = 1:7))

compareTaxAssocsQuantilePlots <- function(taxon = NULL, complication = "anyComplication", simpleLegend = TRUE) {
    tmp2 <- genusProfiles %>%
        filter(taxa == taxon) %>%
        mutate(anyABx = ifelse(anyABx, "ABx", "no ABx")) %>%
        mutate(anyABx = factor(anyABx, levels = c("no ABx", "ABx"))) %>%
        mutate(relAb = log10(relAb + pseudoCount))
    # mutate(postTransplant = ifelse(postTransplant, "post-Tx", "pre-Tx")) %>%
    # mutate(postTransplant = factor(postTransplant, levels = c("pre-Tx", "post-Tx")))
    if (is.logical(tmp2[[complication]])) {
        if (simpleLegend) {
            tmp2[[complication]] <- ifelse(tmp2[[complication]], "comp.", "no comp.")
            tmp2[[complication]] <- factor(tmp2[[complication]], levels = c("comp.", "no comp."))
        } else {
            tmp2[[complication]] <- ifelse(tmp2[[complication]], complication, str_c("no ", complication, sep = "", collapse = ""))
            tmp2[[complication]] <- factor(tmp2[[complication]], levels = c(complication, str_c("no ", complication, sep = "", collapse = "")))
        }
    }
    # Remove NAs
    tmp2 <- tmp2 %>%
        filter(!if_any(all_of(c('anyABx', complication, "relAb")), is.na))
    p <- get_quantile_plot(tmp2, "anyABx", complication, "relAb", plotGroup = NA, expectedNumLevels = 2)
    p <- p +
        ggtitle(str_c("LMM raw p-val for\n",
            taxon, '\n',
            complication, ': ',
            round(lmmsTibble %>% filter(taxa == taxon) %>% rename(v = visitType) %>% filter(v == visitType) %>% filter(outcome == complication) %>% pull(pValTaxon), 4),
            "\n",
            "Comparison: ",
            visitType
        ))
    return(p)
}

compareTaxAssocsScatter <- function(taxon = NULL, outcomeMeasure = "CD") {
    tmp2 <- genusProfiles %>%
        filter(taxa == taxon) %>%
        mutate(anyABx = ifelse(anyABx, "ABx", "no ABx")) %>%
        mutate(anyABx = factor(anyABx, levels = c("no ABx", "ABx"))) %>%
        mutate(relAb = log10(relAb + pseudoCount))
    # mutate(postTransplant = ifelse(postTransplant, "post-Tx", "pre-Tx")) %>%
    # mutate(postTransplant = factor(postTransplant, levels = c("pre-Tx", "post-Tx")))
    tmp2$outcome <- tmp2[[outcomeMeasure]]
    tmp2 <- tmp2[!is.na(tmp2$outcome), ]
    # print(tmp2 %>% select(relAb, outcome, postTransplant, anyABx))
    # return(tmp2)
    p <- ggplot() +
        geom_point(data = tmp2, aes(x = relAb + pseudoCount, y = outcome)) +
        theme_classic() +
        facet_grid(. ~ anyABx) +
        NULL
    p <- p +
        ggtitle(str_c("LMM raw p-val for\n",
            taxon, '\n',
            complication, ': ',
            round(lmmsTibble %>% filter(taxa == taxon) %>% rename(v = visitType) %>% filter(v == visitType) %>% filter(outcome == complication) %>% pull(pValTaxon), 4),
            "\n",
            "Comparison: ",
            visitType
        ))
    return(p)
}



# for (taxon in unique(lmmsTibble$taxa)) {
# for (taxon in lmmsTibble %>% filter(pValTaxon < 0.05) %>% pull(taxa) %>% unique()) {
#     for (complication in unique(lmmsTibble$outcome)) {
#         dir.create(str_c("/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/complication_LMMs_results/", complication), showWarnings = FALSE)
#         yLabel <- ifelse(taxon != "richness", "Relative Abundance (log10)", "richness (log10)")
#         if (complication == "CD") {
#             next
#         } else{
#             ggsave(plot = compareTaxAssocsQuantilePlots(taxon = taxon, complication = complication) + ylab(yLabel), filename = str_c('/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/complication_LMMs_results/',complication, "/", taxon, "_LMM_quantilePlots.pdf"), width = 5, height = 5)
#         }
#     }
# }

# for (taxon in unique(lmmsTibble$taxa)) {

### CAREFUL: THis code is non-functional in that you cannot distinguish between visit and visitMinusOne here...
plots <- list()
for (g in list(
    c("Roseburia", "CDbinary"),
    c("Eubacterium", "CDbinary"),
    c("Enterococcus", "CD"),
    c("Escherichia", "rejection"),
    c("Clostridium", "rejection"),
    c("Alistipes", "rejection")
)) {
    print(g)
    taxon <- g[1]
    complication <- g[2]
    dir.create(str_c("/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/complication_LMMs_results/", complication), showWarnings = FALSE)
    yLabel <- ifelse(taxon != "richness", "Relative Abundance (log10)", "richness (log10)")
    if (complication == "CD") {
        plots[[length(plots) + 1]] <- compareTaxAssocsScatter(taxon = taxon, outcomeMeasure = complication) + ylab("CD ratio")
    } else {
        plots[[length(plots) + 1]] <- compareTaxAssocsQuantilePlots(taxon = taxon, complication = complication) + ylab(yLabel) #            ggsave(plot = compareTaxAssocsQuantilePlots(taxon = taxon, complication = complication) + ylab(yLabel), filename = str_c('/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/complication_LMMs_results/',complication, "/", taxon, "_LMM_quantilePlots.pdf"), width = 5, height = 5)
    }
}

library(patchwork)
ggsave(plot = wrap_plots(plots, ncol = 2, nrow = 3, guides = 'collect'), filename = "/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/complicationAssocsSeleceted.pdf", width = 9, height = 10)

##############################
#### CD RATIO VARIABILITY ####
##############################

metabCDThreshold <- 1
allowDifference <- 1
tpFilter <- 3

cdMetabColors <- c(
    "high" = "#fc5852",
    "low" = "#6bc64a",
    'mixed' = "#6c6c6c"
)

p <- ggplot(outcomeInformation %>%
    mutate(visitNumber = factor(visitNumber, levels = 1:7, ordered = TRUE)) %>%
    group_by(patientID) %>%
    filter(!is.na(CD)) %>%
    arrange(visitNumber) %>%
    filter(visitNumber > tpFilter) %>%
    mutate(visitNumber = factor(visitNumber, levels = levels(visitNumber)[(tpFilter + 1):length(levels(visitNumber))])) %>%
    nest() %>%
    # mutate(varianceCD = map_dbl(data, \(x) {
    #     return(var(x$CD))
    # })) %>%
    mutate(`cdMetabolism` = map_chr(data, \(x) {
        samples <- dim(x)[1]
        case_when(
            sum(x$CD < metabCDThreshold) >= (samples - allowDifference) ~ "high",
            sum(x$CD > metabCDThreshold) >= (samples - allowDifference) ~ "low",
            .default = "mixed"
        )
    }
    )) %T>%
    write_tsv("/g/scb/zeller/karcher/PRISMA/results/CD_metabolism_map.tsv") %>%
    unnest() %>%
    identity()
, aes(x = visitNumber, y = CD)) +
    geom_hline(yintercept = 1, linetype = 'dotted') +
    geom_boxplot(outlier.color = NA) +
    geom_jitter(width = 0.05, height = 0, alpha = 0.3) +
    # geom_path(aes(group = patientID, color = sqrt(varianceCD), alpha = sqrt(varianceCD))) +
    geom_path(aes(group = patientID, color = `cdMetabolism`), alpha = 0.5) +
    theme_classic() +
    # scale_colour_gradient(low = "grey", high = "red") +
    scale_color_manual(values = cdMetabColors
    ) +
    scale_alpha(range = c(0.2, 1)) +
    #        annotate(geom = "text", x = 0.5, y = 1.5, label = "high\nmetab", hjust = 0) +
    #        annotate(geom = "text", x = 0.5, y = 0.5, label = "low\nmetab", hjust = 0) +
    ggtitle(str_c("Regarding cdMetabolism classification:\nAllowing, ",
        allowDifference,
        "sample(s) to disagree with label\nRemoving samples TP <=", tpFilter)) +
    # scale_x_discrete(drop = FALSE) +
    scale_x_discrete_prisma(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    NULL

ggsave(plot = p, filename = str_c("/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/CDOverTime_allowDifference_", allowDifference, ".pdf"), width = 5, height = 5.5)
# ggsave(plot = p, filename = "/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/CDOverTime.png", width = 5, height = 5.5)

(clinicalMetadata %>%
    select(patientID, cyp3a5star3, cyp3a4star22) %>%
    distinct() %>%
    select(-patientID) %>%
    table() %>%
    as.data.frame() %>%
    mutate(across(all_of(c('cyp3a4star22', 'cyp3a5star3')), \(x) factor(x, levels = c(FALSE, TRUE)))) %>%
    ggplot(aes(x = cyp3a4star22, y = Freq, fill = cyp3a5star3)) +
    geom_bar(stat = "identity") +
    scale_x_discrete(drop = FALSE) +
    theme_classic() +
    ylab("Number of patients")) %>%
    ggsave(filename = "/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/cyp_genotypes.pdf", width = 4, height = 2.5)



(read_tsv("/g/scb/zeller/karcher/PRISMA/results/CD_metabolism_map.tsv") %>%
    select(-data) %>%
    left_join(
        clinicalMetadata %>%
            select(patientID, cyp3a5star3, cyp3a4star22) %>%
            distinct()
    ) %>%
    pivot_longer(-c(patientID, `cdMetabolism`)) %>%
    rename(cyp_genotype = name) %>%
    identity() %>%
    select(-patientID) %>%
    relocate(`cdMetabolism`, value) %>%
    table() %>%
    as.data.frame() %>%
    ggplot() + geom_bar(aes(x = value, y = Freq, fill = cdMetabolism), stat = 'identity') +
    theme_classic() +
    facet_grid(. ~ cyp_genotype) +
    scale_fill_manual(values = cdMetabColors) +
    ylab("Number of patients")) %>%
    ggsave(filename = "/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/cyp_genotype_CD_metabolism.pdf")


cdModelDataSmall <- read_tsv("/g/scb/zeller/karcher/PRISMA/results/CD_metabolism_map.tsv") %>%
    select(-data) %>%
    left_join(
        clinicalMetadata %>%
            select(patientID, visit, cyp3a5star3, cyp3a4star22, firstAlbuminMeasurement, ageCategorical, sex, weight, firstHematocritMeasurement) %>%
            # for weight
            filter(visit == 1) %>%
            select(-visit) %>%
            distinct()
    ) %>%
    filter(cdMetabolism != 'mixed') %>%
    mutate(cdMetabolism = factor(cdMetabolism, levels = c('low', 'high'))) %>%
    mutate(`CD-ratio` = factor(ifelse(cdMetabolism == 'low', "high", 'low'), levels = c('low', "high"))) %>%
    mutate(sex = as.factor(sex))

# ATTENTION: I have to artificially include some noise cause otherwise the logistic model won't fit...
cdModelDataSmall$cyp3a5star3[c(1)] <- TRUE
# ATTENTION: I impute albumin/hematocrit with the mean
cdModelDataSmall$firstAlbuminMeasurement[is.na(cdModelDataSmall$firstAlbuminMeasurement)] <- mean(cdModelDataSmall$firstAlbuminMeasurement[!is.na(cdModelDataSmall$firstAlbuminMeasurement)])
# cdModelSmall <- glm(data = cdModelDataSmall ,
#     formula= cdMetabolism ~ cyp3a5star3 + firstAlbuminMeasurement + ageCategorical + firstHematocritMeasurement, family = 'binomial', maxit = 100000)
# cdModelSmall <- randomForest(cdMetabolism ~ cyp3a5star3 + firstAlbuminMeasurement + ageCategorical + firstHematocritMeasurement, data=cdModelDataSmall, proximity=TRUE)
library(randomForest)
summary(cdModelSmall)

ps <- list()
set.seed(2)
for (patientID in cdModelDataSmall$patientID) {
    test <- cdModelDataSmall[cdModelDataSmall$patientID == patientID, ]
    train <- cdModelDataSmall[cdModelDataSmall$patientID != patientID, ]
    cdModelSmall <- randomForest(ntree = 50, mtry = 1, cdMetabolism ~ cyp3a5star3 + firstAlbuminMeasurement + ageCategorical + firstHematocritMeasurement + sex + weight, data = train, proximity = TRUE)
    p <- predict(cdModelSmall, test, type = 'prob')[, 1]
    ps[[length(ps) + 1]] <- p
}
rocObjectModelSmall <- roc(predictor = unlist(ps), response = as.numeric(cdModelDataSmall$cdMetabolism))
rocObjectModelSmall

res <- list()
preTransplantProfiles <- profiles %>%
    inner_join(data.frame(visit = c(1))) %>%
    inner_join(importantTaxa %>% rename(genus = taxa))
# preTransplantProfiles <- profiles %>%
#     inner_join(data.frame(visit = c(1))) %>%
#     group_by(sampleID, family, relAb, PSN, visit) %>%
#     mutate(relAb = 10^(relAb) - pseudoCount) %>%
#     summarize(relAb = sum(relAb)) %>%
#     group_by(family) %>%
#     filter(mean(relAb > 0.01) > 0.1)
# preTransplantProfiles <- profiles %>%
#     inner_join(data.frame(visit = c(1))) %>%
#     group_by(sampleID, phylum, relAb, PSN, visit) %>%
#     mutate(relAb = 10^(relAb) - pseudoCount) %>%
#     summarize(relAb = sum(relAb)) %>%
#     group_by(phylum) %>%
#     filter(mean(relAb > 0.005) > 0.1)


# This doesn't really make much sense anymore, and is also overfitting (see later)
modelDataAll <- list()
for (g in unique(preTransplantProfiles$genus)) {
    # for (g in unique(preTransplantProfiles$family)) {
    # for (g in unique(preTransplantProfiles$phylum)) {
    cdModelData <- cdModelDataSmall %>%
        left_join(preTransplantProfiles %>%
            filter(genus == g) %>%
            select(genus, relAb, PSN) %>%
            rename(patientID = PSN))
    cdModel <- glm(data = cdModelData,
        # relAb here is already log10-scaled...
        formula = cdMetabolism ~ cyp3a5star3 + firstAlbuminMeasurement + ageCategorical + firstHematocritMeasurement + relAb, family = 'binomial')

    res[[length(res) + 1]] <- cdModel
    names(res)[length(res)] <- g

    modelDataAll[[length(modelDataAll) + 1]] <- cdModelData
    names(modelDataAll)[length(modelDataAll)] <- g
}
resTibble <- tibble(genus = names(res), models = res) %>%
    mutate(summary = map(models, summary)) %>%
    mutate(cyp3a5star3_pvalue = map_dbl(summary, \(x) {
        x$coefficients[rownames(x$coefficients) == "cyp3a5star3TRUE", 4]
    })) %>%
    mutate(taxon_pvalue = map_dbl(summary, \(x) {
        x$coefficients[rownames(x$coefficients) == "relAb", 4]
    })) %>%
    mutate(taxon_estimate = map_dbl(summary, \(x) {
        x$coefficients[rownames(x$coefficients) == "relAb", 1]
    })) %>%
    left_join(profiles %>% ungroup() %>% select(genus, phylum) %>% distinct(), by = c('genus' = 'genus')) %>%
    relocate(genus, phylum) %>%
    arrange(taxon_pvalue)

library(gridBase)
library(gridExtra)


(resTibble %>%
    arrange(taxon_pvalue) %>%
    head(20) %>%
    mutate(genus = factor(genus, levels = genus)) %>%
    ggplot(aes(x = genus, y = taxon_pvalue, fill = phylum)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    geom_bar(stat = 'identity')) %>%
    ggsave(filename = "/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/glm_cd_cyp_tax_profiles.pdf", width = 8, height = 3.5)

illustrate_taxon_hit <- function(modelData = NULL, taxon = NULL) {
    modelData <- modelData %>%
        filter(genus == taxon) %>%
        rename(`Tacrolimus\nmetabolism` = cdMetabolism) %>%
        mutate(relAb = (10^relAb) * 100)
    p <- ggplot() +
        geom_boxplot(
            data = modelData,
            aes(x = `CD-ratio`, y = relAb, fill = `CD-ratio`), outlier.color = NA) +
        geom_jitter(
            data = modelData,
            aes(x = `CD-ratio`, y = relAb, fill = `CD-ratio`), position = position_jitter()) +
        theme_classic() +
        ylab("Bacterial\nrelative abundance [%]") +
        scale_fill_manual(values = c('low' = "#4a5dca", "high" = "#d43e3e")) +
        scale_y_continuous(trans = 'log10', limits = c(0.005, 5))
    return(p)
}


plots <- list()
for (g in c("Erysipelatoclostridium", "Enterococcus", "Roseburia", "Coprococcus")) {
    # ggsave(filename = str_c("/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/", g, ".pdf"), width = 5, height = 4)
    plots[[length(plots) + 1]] <- illustrate_taxon_hit(do.call('rbind', modelDataAll), g) + ggtitle(g) + theme(plot.title = element_text(size = 8, face = "bold"))
}

ggsave(plot = wrap_plots(plots, guides = 'collect', nrow = 1),
    filename = "/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/cd_metabolism_hits.pdf", width = 9, height = 3.5)

candidateGenera <- c("Erysipelatoc.", "Enterococcus", "Roseburia", "Coprococcus")

cdModelDataBig <- cdModelDataSmall %>%
    left_join(preTransplantProfiles %>%
        filter(genus %in% candidateGenera | genus %in% "Erysipelatoclostridium") %>%
        select(genus, relAb, PSN) %>%
        rename(patientID = PSN) %>%
        pivot_wider(id_cols = patientID, names_from = genus, values_from = relAb)) %>%
    rename(`Erysipelatoc` = Erysipelatoclostridium)


# cdModels <- list()
# for (genus in candidateGenera) {
# cdModelBig <- glm(data = cdModelDataBig,
#     formula = as.formula(
#         str_c("cdMetabolism ~ cyp3a5star3 + firstAlbuminMeasurement + ageCategorical + firstHematocritMeasurement + ", str_c(genus, sep = " + ", collapse = " + "))), family = 'binomial', maxit = 100000)
ps <- list()
set.seed(2)
for (patientID in cdModelDataBig$patientID) {
    test <- cdModelDataBig[cdModelDataBig$patientID == patientID, ]
    train <- cdModelDataBig[cdModelDataBig$patientID != patientID, ]
    cdModelBig <- randomForest(ntree = 50, mtry = 1, cdMetabolism ~ cyp3a5star3 + firstAlbuminMeasurement + ageCategorical + firstHematocritMeasurement + sex + weight + Roseburia + Coprococcus + Enterococcus + Erysipelatoc, data = train, proximity = FALSE)
    p <- predict(cdModelBig, test, type = 'prob')[, 1]
    ps[[length(ps) + 1]] <- p
}
rocObjectModelBig <- roc(predictor = unlist(ps), response = as.numeric(cdModelDataBig$cdMetabolism))
print(rocObjectModelBig$auc)

cdModels <- tibble(
    model = list(cdModelSmall, cdModelBig),
    group = list("Clinical model", "CM + microbiome"),
    roc = list(rocObjectModelSmall, rocObjectModelBig))
cdModels <- cdModels %>%
    mutate(specs = map(roc, \(x) {
        return(data.frame(TPR = x$specificities, FPR = 1 - x$sensitivities))
    })) %>%
    mutate(auc = map_dbl(roc, \(x) {
        return(x$auc)
    })) %>%
    mutate(group = factor(group, levels = rev(c('Clinical model', "CM + microbiome")), ordered = TRUE)) %>%
    arrange(group) %>%
    mutate(y = seq(0.125, 0.05, length.out = length(levels(group)))) %>%
    mutate(label = str_c(group, ": ", round(auc, 3))) %>%
    rename(Features = group)

# colVec <- # Define colors
grey_color <- "#888888" # Grey
blue_color <- "#3498db" # Blue
green_color <- "#2ecc71" # Green
purple_color <- "#9b59b6" # Purple
red_color <- "#e74c3c" # Red

# Display the colors
colors <- c(blue_color, red_color, green_color, red_color, grey_color)
names(colors) <- levels(cdModels$Features)

pClinical <- ggplot() +
    geom_line(data = cdModels %>%
        select(Features, specs) %>%
        unnest() %>%
        filter(Features == 'Clinical model'), aes(x = FPR, y = TPR, group = Features, color = Features)) +
    theme_classic() +
    scale_color_manual(values = colors) +
    xlab("False Positive Rate") +
    ylab("True Positive Rate") +
    geom_text(data = cdModels %>%
        filter(Features == 'Clinical model'), aes(x = 0.1, y = y, label = label), inherit.aes = FALSE, hjust = 0)

ggsave(
    plot = pClinical,
    filename = "/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/cdMetabolismPredictionOnlyClinical.pdf", width = 5, height = 3.25)

pC <- ggplot() +
    geom_line(data = cdModels %>%
        select(Features, specs) %>%
        unnest() %>%
        filter(Features == 'CM + microbiome'), aes(x = FPR, y = TPR, group = Features, color = Features)) +
    theme_classic() +
    scale_color_manual(values = colors) +
    xlab("False Positive Rate") +
    ylab("True Positive Rate") +
    geom_text(data = cdModels %>%
        filter(Features == 'CM + microbiome'), aes(x = 0.1, y = y, label = label), inherit.aes = FALSE, hjust = 0)

ggsave(
    plot = pC,
    filename = "/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/cdMetabolismPredictionClinicalOnlyMicrobiome.pdf", width = 5, height = 3.25)

pAll <- ggplot() +
    geom_line(data = cdModels %>%
        select(Features, specs) %>%
        unnest() %>%
        # mutate(lt = ifelse(Features == 'clinical model', "2", "1")) %>%
        arrange(desc(Features)), aes(x = FPR, y = TPR, group = Features, color = Features)) +
    theme_classic() +
    scale_color_manual(values = colors) +
    xlab("False Positive Rate") +
    ylab("True Positive Rate") +
    # annotate('text', x = 0.1, y = 0.4, label = str_c("DeLong's test p-value:\n", round(roc.test(cdModels$roc[[1]], cdModels$roc[[2]])$p.value, 3)), hjust = 0) +
    geom_text(data = cdModels, aes(x = 0.1, y = y, label = label), inherit.aes = FALSE, hjust = 0)

ggsave(
    plot = pAll,
    filename = "/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/cdMetabolismPrediction.pdf", width = 5, height = 3.25)


library(pROC)


# cdModelDataBig %>%
#     select(patientID, Roseburia, Enterococcus, `Erysipelatoc.`, Coprococcus) %>%
#     left_join(clinicalMetadata)
data <- rbind(profiles %>%
    select(PSN, visit, genus, relAb) %>%
    rename(patientID = PSN) %>%
    mutate(relAb = (10^relAb) * 100),
outcomeInformation %>%
    select(patientID, visitNumber, CD) %>%
    filter(!is.na(CD)) %>%
    rename(visit = visitNumber) %>%
    mutate(genus = "CD") %>%
    rename(relAb = CD)) %>%
    arrange(patientID, visit) %>%
    # filter(genus == "CD") %>%
    inner_join(data.frame(genus = c("CD", candidateGenera, "Erysipelatoclostridium"))) %>%
    group_by(genus) %>%
    # mutate(relAb = scale(relAb)[, 1]) %>%
    identity()
# colnames(data)[colnames(data) == "relAb[,1]"] <- 'relAb'

get_line_plot <- function(PSN = NULL) {

    tmp <- c(
        "4" = "week 1",
        "5" = "week 4",
        "6" = "Month 3",
        "7" = "Month 6"
    )

    data$measurement <- data$genus

    data <- data %>%
        filter(patientID == PSN) %>%
        filter(visit > 3)

    data <- data %>%
        filter(measurement %in% c("CD", "Roseburia", "Coprococcus", "Erysipelatoclostridium"))

    scaleFactor <- max(data %>% filter(measurement != "CD") %>% pull(relAb)) / max(data %>% filter(measurement == "CD") %>% pull(relAb))

    p <- ggplot() +
        geom_line(data = data %>% filter(measurement != "CD"), aes(x = visit, y = relAb, color = measurement, group = measurement)) +
        geom_line(data = data %>% filter(measurement == "CD"), aes(x = visit, y = relAb * scaleFactor, color = measurement, group = measurement), color = 'grey', linetype = 2) +
        theme_classic() +
        ylab("Bacterial\nrelative abundances [%]") +
        # ylab("bacterial abundance/\nCD ratio") +
        scale_y_continuous(name = "Bacterial\nrelative abundances [%]", sec.axis = sec_axis(~ . / scaleFactor, name = "CD-ratio")) +
        ggtitle(str_c(PSN)) +
        scale_x_continuous(labels = tmp, breaks = as.numeric(names(tmp))) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1))

    return(p)
}

plots <- list()
for (PSN in unique(data$patientID)) {
    p <- get_line_plot(PSN = PSN)
    plots[[length(plots) + 1]] <- p
    llAppend(plots, p, PSN)
}
names(plots) <- unique(data$patientID)

map2(plots, names(plots), \(p, n){
    ggsave(plot = p, filename = str_c('/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/cd_lineplots_per_patient/', n, ".pdf"), width = 5.5, height = 3.5)
})

# library(ggplot2)

# scaleFactor <- max(mtcars$cyl) / max(mtcars$hp)

# ggplot(mtcars, aes(x=disp)) +
#   geom_smooth(aes(y=cyl), method="loess", col="blue") +
#   geom_smooth(aes(y=hp * scaleFactor), method="loess", col="red") +
#   scale_y_continuous(name="cyl", sec.axis=sec_axis(~./scaleFactor, name="hp")) +
#   theme(
#     axis.title.y.left=element_text(color="blue"),
#     axis.text.y.left=element_text(color="blue"),
#     axis.title.y.right=element_text(color="red"),
#     axis.text.y.right=element_text(color="red")
#   )

# a <- rownames(as.matrix(pairwiseDistances))
# b <- colnames(as.matrix(pairwiseDistances))
tmp <- pairwiseDistances %>%
    # tmp <- pairwiseDistancesIdentityEuclidean %>%
    as.matrix()
tmp[lower.tri(tmp)] <- NA
tmp <- tmp %>%
    as.data.frame() %>%
    rownames_to_column('sampleID') %>%
    pivot_longer(-sampleID) %>%
    rename(sampleID1 = sampleID, sampleID2 = name, distance = value) %>%
    left_join(pcoa %>% select(sampleID, PSN, visit) %>% mutate(visit = as.numeric(as.character(visit))), by = c("sampleID1" = 'sampleID')) %>%
    left_join(pcoa %>% select(sampleID, PSN, visit) %>% mutate(visit = as.numeric(as.character(visit))), by = c("sampleID2" = 'sampleID'), suffix = c("1", "2")) %>%
    filter(PSN1 == PSN2) %>%
    filter(!visit1 %in% c(1, 2, 3)) %>%
    filter(!visit2 %in% c(1, 2, 3)) %>%
    filter(!is.na(distance)) %>%
    filter(visit1 != visit2)

# cdShift <- (outcomeInformation %>%
CDDiffVsMicrobiomeDiff <- outcomeInformation %>%
    # mutate(visitNumber = factor(visitNumber, levels = 1:7, ordered = TRUE)) %>%
    group_by(patientID) %>%
    filter(!is.na(CD)) %>%
    arrange(visitNumber) %>%
    filter(!visitNumber %in% c(1, 2, 3)) %>%
    nest() %>%
    mutate(d = map(data, \(x) data.frame(refVisitNumber = x$visitNumber, visitNumber = c(NA, x$visitNumber[1:(length(x$visitNumber) - 1)])))) %>%
    left_join(tmp %>% group_by(PSN1) %>% nest() %>% rename(patientID = PSN1) %>% rename(pwd = data)) %>%
    mutate(pwd = map2(pwd, d, \(p, dd) {
        inner_join(dd %>%
            mutate(refVisitNumber = as.numeric(as.character(refVisitNumber))), p, by = c("refVisitNumber" = "visit2", "visitNumber" = "visit1"))
    })) %>%
    mutate(data = map(data, \(x) {
        x <- x %>%
            select(visitNumber, CD) %>%
            mutate(CDpreviousVisit = c(NA, CD[1:(length(CD) - 1)])) %>%
            mutate(CDShiftWRTPreviousVisit = CD - CDpreviousVisit)
    })) %>%
    mutate(data = map2(data, pwd, function(d, p) {
        inner_join(d, p, by = c("visitNumber" = "refVisitNumber"))
    })) %>%
    select(patientID, data) %>%
    unnest() %>%
    # select(patientID, datavisitNumber, CDShiftWRTPreviousVisit) %>%
    select(patientID, visitNumber, CDShiftWRTPreviousVisit, distance) %>%
    rename(microbiomeDistanceWRTToPreviousVisit = distance) %>%
    mutate(visitNumber = factor(visitNumber, levels = 1:7))


(CDDiffVsMicrobiomeDiff %>%
    ggplot() +
    geom_point(aes(x = microbiomeDistanceWRTToPreviousVisit, y = abs(CDShiftWRTPreviousVisit), color = visitNumber)) +
    geom_text_repel(aes(x = microbiomeDistanceWRTToPreviousVisit, y = abs(CDShiftWRTPreviousVisit), label = patientID)) +
    theme_classic()) %>%
    # ggsave(filename = '/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/CDChangeAgainstMicrobiomeChange__identity__Euclidean.pdf', width = 8, height = 8)
    ggsave(filename = '/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/CDChangeAgainstMicrobiomeChange__log10__Euclidean_with_patientID.pdf', width = 8, height = 8)
# ggsave(filename = '/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/CDChangeAgainstMicrobiomeChange__identity__Euclidean.png', width = 8, height = 8)

(CDDiffVsMicrobiomeDiff %>%
    ggplot() +
    geom_point(aes(x = microbiomeDistanceWRTToPreviousVisit, y = abs(CDShiftWRTPreviousVisit), color = visitNumber)) +
    # geom_text_repel(aes(x = microbiomeDistanceWRTToPreviousVisit, y = abs(CDShiftWRTPreviousVisit), label = patientID)) +
    theme_classic()) %>%
    # ggsave(filename = '/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/CDChangeAgainstMicrobiomeChange__identity__Euclidean.pdf', width = 8, height = 8)
    ggsave(filename = '/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/CDChangeAgainstMicrobiomeChange__log10__Euclidean.pdf', width = 8, height = 8)


(CDDiffVsMicrobiomeDiff %>%
    mutate(CDShiftGroup = ifelse(CDShiftWRTPreviousVisit > 1, "large (>1)", "small (<1)")) %>%
    ggplot() +
    geom_boxplot(aes(x = CDShiftGroup, y = microbiomeDistanceWRTToPreviousVisit)) +
    theme_classic() +
    xlab("CD difference WRT to previous time point")) %>%
    ggsave(filename = '/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/CDChangeAgainstMicrobiomeChange__log10__Euclidean__discrete_boxplot.pdf', width = 4, height = 4)

# (CDDiffVsMicrobiomeDiff %>%
#             ggplot() +
#             geom_point(aes(x = microbiomeDistanceWRTToPreviousVisit, y = abs(CDShiftWRTPreviousVisit), color = patientID)) +
#             geom_text_repel(aes(x = microbiomeDistanceWRTToPreviousVisit, y = abs(CDShiftWRTPreviousVisit), label = patientID)) +
#     theme_classic()) %>%
#     ggsave(filename = '/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/CDChangeAgainstMicrobiomeChange__log10__Euclidean_col_by_patientID.pdf', width = 8, height = 8)
# ggsave(filename = '/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/CDChangeAgainstMicrobiomeChange__identity__Euclidean_col_by_patientID.png', width = 8, height = 8)


(outcomeInformation %>%
    mutate(visit = as.factor(visitNumber)) %>%
    inner_join(profiles %>% group_by(PSN, visit) %>% summarize(richness = sum(relAb > log10(pseudoCount))) %>%
        # For ordination, limit richness to 150 max (there is a single outlier sample...)
        mutate(richness = ifelse(richness > 150, 150, richness)) %>%
        mutate(visit = as.factor(visit)), by = c("patientID" = "PSN", "visit")) %>%
    select(patientID, visit, richness, CD) %>%
    ggplot(aes(x = richness, y = CD)) +
    geom_point() +
    theme_classic()) %>%
    ggsave(filename = "/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/CDagainstRichness.pdf", width = 4.5, height = 4.5)

library(ggsignif)

tmp <- outcomeInformation %>%
    mutate(visit = as.factor(visitNumber)) %>%
    inner_join(profiles %>% group_by(PSN, visit) %>% summarize(richness = sum(relAb > log10(pseudoCount))) %>%
        # For ordination, limit richness to 150 max (there is a single outlier sample...)
        mutate(richness = ifelse(richness > 150, 150, richness)) %>%
        mutate(visit = as.factor(visit)), by = c("patientID" = "PSN", "visit")) %>%
    filter(!is.na(CD))
q <- quantile(tmp$richness, probs = c(0.33, 0.66))
tmp <- tmp %>%
    mutate(richnessGroup = case_when(
        richness < q[1] ~ "lower tercile",
        richness > q[1] & richness < q[2] ~ "middle tercile",
        richness > q[2] ~ "upper tercile",
    )) %>%
    mutate(richnessGroup = as.factor(richnessGroup))

(ggplot(data = tmp, aes(x = richnessGroup, y = CD)) +
    geom_boxplot() +
    geom_signif(comparisons = list(
        # c("lower tercile", 'middle tercile'),
        # c("lower tercile", 'upper tercile'),
        c("middle tercile", 'upper tercile'))) +
    theme_classic() +
    xlab("CD ratio (discrete)")
) %>%
    ggsave(filename = "/g/scb/zeller/karcher/PRISMA/plots/KLGPG_221206/CDagainstRichnessDiscrete.pdf", width = 4.5, height = 4.5)
