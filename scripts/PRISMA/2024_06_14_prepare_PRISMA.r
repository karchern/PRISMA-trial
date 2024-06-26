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

obj_path <- here('objects/PRISMA_idtaxa.rdata')

print("Loading profiles...")
# profiles_interim <- readRDS(here('profiles/16S/240718_PRISMA_INTERIMS_KLGPG_WITHOUT_QC/Results/collated/res_mapseq.rds'))
profiles_interim <- readRDS(here('profiles/16S/240718_PRISMA_INTERIMS_KLGPG_WITHOUT_QC/Results/collated/res_IDTaxa.rds'))
# profiles_tmp_batch_A <- readRDS(here('profiles/16S/240718_PRISMA_MODELLING_BATCHA_WITHOUT_QC/Results/collated/res_mapseq.rds'))
profiles_tmp_batch_A <- readRDS(here('profiles/16S/240718_PRISMA_MODELLING_BATCHA_WITHOUT_QC/Results/collated/res_IDTaxa.rds'))
# profiles_tmp_batch_B <- readRDS(here('profiles/16S/240718_PRISMA_MODELLING_BATCHB_WITHOUT_QC/Results/collated/res_mapseq.rds'))
profiles_tmp_batch_B <- readRDS(here('profiles/16S/240718_PRISMA_MODELLING_BATCHB_WITHOUT_QC/Results/collated/res_IDTaxa.rds'))
profiles_tmp_batch_A_genus <- .f_resolve_taxonomy(profiles_tmp_batch_A, 'genus')
profiles_tmp_batch_B_genus <- .f_resolve_taxonomy(profiles_tmp_batch_B, 'genus')
profiles_interim <- .f_resolve_taxonomy(profiles_interim, 'genus')
stopifnot(!any(colnames(profiles_tmp_batch_A_genus) %in% colnames(profiles_tmp_batch_B_genus)))
profiles <- map2(list('modellingBatchA', 'modellingBatchB', 'interim'), list(profiles_tmp_batch_A_genus, profiles_tmp_batch_B_genus, profiles_interim), \(batch_name, x) x %>%
    as.data.frame() %>%
    rownames_to_column('taxon') %>%
    pivot_longer(-taxon) %>%
    rename(sampleID = name, count = value) %>%
    mutate(batch = batch_name)) %>%
    do.call('rbind', .) %>%
    pivot_wider(id_cols = c(taxon), names_from = c(sampleID, batch), values_from = count, values_fill = 0, names_sep = "___") %>%
    as.data.frame() %>% column_to_rownames('taxon') %>% as.matrix()

# meta <- read_tsv(here("data/16S_metadata/221227_PRISMA_16S_Overview.tsv")) %>%
meta <- read_tsv(here("data/16S_metadata/221227_PRISMA_16S_modelling_cohort_Batch_A_and_Batch_B_overview.tsv"), show_col_types = FALSE) %>%
    select(PSN, Visit, ID, Sample_ID, `Sequencing Batch for 16S miSeq`) %>%
    mutate(batch = case_when(
        `Sequencing Batch for 16S miSeq` == "A" ~ "modellingBatchA",
        `Sequencing Batch for 16S miSeq` == "B" ~ "modellingBatchB",
        .default = NA
    )) %>%
    select(-`Sequencing Batch for 16S miSeq`) %>%
    rbind(read_tsv(here("data/16S_metadata/221227_PRISMA_16S_Overview.tsv"), show_col_types = FALSE) %>%
        select(PSN, Visit, ID, Sample_ID) %>%
        mutate(batch = "interim")) %>%
    rename(sampleID = Sample_ID, visit = Visit) %>%
    mutate(PSN = str_replace(PSN, "Mue", "M")) %>%
    mutate(PSN = str_replace(PSN, "MÜ", "M")) %>%
    mutate(PSN = str_replace(PSN, "ä", "ae")) %>%
    mutate(PSN = str_replace(PSN, "ö", "oe")) %>%
    mutate(PSN = str_replace(PSN, "ü", "ue")) %>%
    mutate(PSN = str_replace(PSN, "Ä", "AE")) %>%
    mutate(PSN = str_replace(PSN, "Ö", "OE")) %>%
    mutate(PSN = str_replace(PSN, "Ü", "UE")) %>%
    mutate(PSN = str_replace(PSN, "NTXM", "NZMU")) %>%
    mutate(PSN = str_replace(PSN, "-0", "-"))
# Remove funny colnames with encoding
meta <- meta[, !str_detect(colnames(meta), "Conc")]
meta <- meta[, !str_detect(colnames(meta), "10 ng")]
# fullTax <- read_tsv(here("data/MAPseq_AlessioCurated.tax"), col_names = F, show_col_types = FALSE)
fullTax <- read_tsv(here('data/gtdbk_tax_r207.tab'), col_names = F, show_col_types = FALSE)
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


if ("Bacteria" %in% rownames(profiles)) {
    profiles <- profiles[!rownames(profiles) == "Bacteria", ]
}
colnames(profiles) <- map_chr(colnames(profiles), function(x) str_replace(x, ".*lane1", ""))
profiles <- profiles[, !str_detect(colnames(profiles), "water")]
colnames(profiles) <- map_chr(colnames(profiles), function(x) str_replace(x, "MG", "MG_"))
colnames(profiles) <- map_chr(colnames(profiles), function(x) str_replace(x, "[.]singles$", ""))

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
    # rename(sampleID = name, count = value) %>%
    mutate(sampleID = str_split_fixed(name, "___", n = 2)[, 1]) %>%
    mutate(batch = str_split_fixed(name, "___", n = 2)[, 2]) %>%
    group_by(genus, sampleID, batch) %>%
    summarize(count = sum(value)) %>%
    inner_join(meta, by = c('sampleID', 'batch'))

print("Computing depths...")
depths <- profiles %>%
    group_by(sampleID, batch, visit) %>%
    summarize(totalGenusReadCount = sum(count))

print("Plotting depth histogram...")
depth_histo <- ggplot(depths %>%
    mutate(visit = as.factor(visit))) +
    theme_classic() +
    geom_histogram(aes(x = totalGenusReadCount), alpha = 1, position = "identity") +
    geom_vline(xintercept = rarefactionDepth) +
    geom_text(data = depths %>%
        mutate(visit = as.factor(visit)) %>%
        group_by(visit, batch) %>%
        filter(totalGenusReadCount < rarefactionDepth) %>%
        tally(), aes(x = rarefactionDepth * 0.5, y = 5, label = n)) +
    facet_grid(visit ~ batch) +
    ylab("Number of\nsamples") +
    xlab("Sequencing depth") +
    theme_presentation() +
    scale_y_continuous(breaks = c(1, 3, 5, 7, 9, 11))

ggsave(plot = depth_histo, filename = here("plots/KLGPG_221206/depth_histogram.pdf"), width = 8, height = 6)

highDepthSamples <- depths %>%
    mutate(visit = as.factor(visit)) %>%
    group_by(sampleID, batch) %>%
    filter(totalGenusReadCount >= rarefactionDepth) %>%
    select(sampleID)

print("Keeping samples with reasonably depth...")
print(str_c("Rarefaction depth: ", rarefactionDepth))
for (B in c("modellingBatchA", "modellingBatchB", "interim")) {
    n <- profiles %>% ungroup() %>% filter(batch == B) %>% anti_join(highDepthSamples, by = c("sampleID", "batch")) %>% select(sampleID) %>% distinct() %>% nrow()
    N <- profiles %>% ungroup() %>% filter(batch == B) %>% select(sampleID) %>% distinct() %>% nrow()
    f <- n / N
    print(str_c("Removing samples with low depth in batch ", B, ": ", n, ' corresponding to ', round(f, 3), ' of all samples'))
}

profiles <- profiles %>%
    inner_join(highDepthSamples, by = c("sampleID", "batch"))

print("Rarefying...")
set.seed(112312)
profiles <- profiles %>%
    pivot_wider(id_cols = c(sampleID, batch), names_from = genus, values_from = count) %>%
    mutate(tmp = str_c(sampleID, batch, sep = "___")) %>%
    relocate(tmp) %>%
    column_to_rownames("tmp") %>%
    select(-sampleID, -batch) %>%
    rrarefy(sample = rarefactionDepth) %>%
    as.data.frame() %>%
    rownames_to_column('tmp') %>%
    mutate(sampleID = str_split_fixed(tmp, "___", n = 2)[, 1]) %>%
    mutate(batch = str_split_fixed(tmp, "___", n = 2)[, 2]) %>%
    select(-tmp) %>%
    pivot_longer(-c(sampleID, batch)) %>%
    rename(genus = name, count = value)

profilesCountsWithUnresolved <- profiles

profiles <- profiles %>%
    group_by(sampleID, batch) %>%
    mutate(relAb = count / sum(count)) %>%
    select(-count) %>%
    filter(genus != "not_resolved") %>%
    filter(genus != "-1") # Should never happen - but just in case :)

print("Plotting boxplots of relative abundances annotated to genus level")
cumRelAbToGenus <- profiles %>%
    left_join(meta %>% select(visit, sampleID, batch), by = c("sampleID", "batch")) %>%
    group_by(sampleID, batch, visit) %>%
    summarize(sumRelAb = sum(relAb))

cumRelAbToGenus_box <- ggplot(cumRelAbToGenus %>%
    mutate(visit = factor(visit, levels = 1:7))) +
    theme_classic() +
    geom_boxplot(aes(x = visit, y = sumRelAb, fill = batch)) +
    ylab("Fraction reads\nannotated to Genus-level") +
    theme_classic()

ggsave(plot = cumRelAbToGenus_box, filename = here("plots/KLGPG_221206/cumRelAbToGenus_boxplot.pdf"), width = 6, height = 3.5)

print("Computing pairwise distances and PCOA object...")

profiles <- profiles %>%
    mutate(relAbOrig = relAb) %>%
    mutate(relAb = log10(relAb + pseudoCount)) %>%
    inner_join(meta, by = c("sampleID", "batch")) %>%
    left_join(fullTax %>% mutate(genus = str_replace(genus, "^g__", "")), by = 'genus')

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

pairwiseDistances <- pivot_wider(profiles, id_cols = genus, names_from = c(sampleID, batch), values_from = relAb, names_sep = "___") %>%
    column_to_rownames("genus") %>%
    as.data.frame() %>%
    as.matrix() %>%
    t() %>%
    vegdist(method = "euclidean", k = 2)

pairwiseDistancesIdentityEuclidean <- pivot_wider(profiles, id_cols = genus, names_from = c(sampleID, batch), values_from = relAbOrig, names_sep = "___") %>%
    column_to_rownames("genus") %>%
    as.data.frame() %>%
    as.matrix() %>%
    t() %>%
    vegdist(method = "euclidean", k = 2)

pcoa <- cmdscale(pairwiseDistances) %>%
    as.data.frame() %>%
    rownames_to_column("tmp") %>%
    mutate(sampleID = str_split_fixed(tmp, "___", n = 2)[, 1]) %>%
    mutate(batch = str_split_fixed(tmp, "___", n = 2)[, 2]) %>%
    select(-tmp) %>%
    left_join(meta, by = c("sampleID", 'batch')) %>%
    as_tibble() %>%
    mutate(visit = as.character(visit)) %>%
    mutate(visit = factor(visit, levels = 1:7))


#################################################################################
# 221024: Integrate clinical metadata and do first analysis of interims cohort
#################################################################################

print("Preparing clinical metadata...")
set.seed(2)
outcomeInformationInterim <- clean_patient_clinical_metadata(read_csv('/g/scb/zeller/karcher/PRISMA/data/16S_metadata/221024_PRISMA_clinical_metadata_names_fixed_ACTUALLY_NEVERMIND_JUST_DO_IT_YOURSELF.csv')) %>%
    filter(!is.na(v62_visit_number) & !is.na(v61_visit_date))
outcomeInformation <- clean_patient_clinical_metadata(read_csv('/g/scb/zeller/karcher/PRISMA/data/16S_metadata/221024_PRISMA_clinical_metadata_BatchA_BatchB.csv') %>%
    mutate(v65_pat_id = ifelse(v65_pat_id == "jobe-nzmu-08", "jobl-nzmu-08", v65_pat_id)) %>%
    # Some patients might have weird, mostly emptry entries. According to maral this can go.
    filter(!is.na(v62_visit_number) & !is.na(v61_visit_date)) %>%
    # This shit is only for RoVo, who has 2 metadata entries for visit 6 and I cannot be fucked this shit anymore
    group_by(v65_pat_id, v62_visit_number) %>%
    sample_n(1) %>%
    ungroup(), how = 'from_maral') %>%
    # the metadata file contains the metadata for the entire modelling cohort, not just the new batches. so this
    anti_join(outcomeInformationInterim %>% select(v65_pat_id) %>% distinct()) %>%
    mutate(v13_dob = as.Date(v13_dob))
outcomeInformation <- rbind(outcomeInformation, outcomeInformationInterim)
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
    mutate(birthday = as.Date(birthday)) %>%
    mutate(age = age_calc(birthday, as.Date("2023-11-07"), "years")) %>%
    mutate(ageCategorical = ifelse(age > 18, "adult", 'non-adult'))

stopifnot(all(outcomeInformation$patientID %in% profiles$PSN))
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
# abxInfo <- read.table(here("data/16S_metadata/221030_PRISMA_antiinfectives_metadata_v2.csv', sep = ");", header =TRUE)
abxInfo <- read.table(here('data/16S_metadata/antiinfectives_metadata.csv'), comment.char = "#", sep = ",", header = TRUE) %>%
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
# it will NOT contain outcomes (those will be stored in outcomeInformation)
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
        CDbinary) %>%
    # I'm kapping CD ratio at 5
    mutate(CD = ifelse(CD > 5, 5, CD))
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

importantTaxaGenus <- profiles %>%
    mutate(relAb = (10^relAb) - pseudoCount) %>%
    mutate(taxa = genus) %>%
    mutate(taxa = as.character(taxa)) %>%
    group_by(taxa, PSN, visit) %>%
    summarize(relAb = sum(relAb)) %>%
    group_by(taxa) %>%
    summarize(m = mean(relAb > pseudoCount) > 0.2, mm = any(relAb > 0.01)) %>%
    filter(m & mm) %>%
    select(taxa) %>%
    # filter(taxa %in% c("Roseburia", "Coprococcus", "Anaerostipes", "Enterococcus"))
    identity()

taxa_failed <- profiles %>%
    filter(is.na(class)) %>%
    ungroup() %>%
    select(genus) %>%
    distinct()

profiles <- profiles %>%
    anti_join(taxa_failed)

print(str_c("Saving all objects for downstream analysis to object", obj_path))

dataList <- c('meta', 'profiles', "pcoa", 'pairwiseDistances', 'pairwiseDistancesIdentityEuclidean', "outcomeInformation", "clinicalMetadata", "orderDFPatientID", "fullTax", "importantTaxaGenus")
save(list = dataList, file = obj_path)
