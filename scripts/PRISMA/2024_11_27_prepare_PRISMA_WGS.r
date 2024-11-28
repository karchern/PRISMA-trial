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

taxonomy_annot <- "ncbi_motus"
# taxonomy_annot <- "ncbi_mapseq"
# taxonomy_annot <- "gtdb_idtaxa"

obj_path <- here(str_c('objects/PRISMA_', taxonomy_annot, '.rdata'))

print("Loading profiles...")

PRISMA_interim_Batch1 <- readRDS(here("profiles/WGS/230111_PRISMA_Batch1_P3/Results/collated/res_mOTUs.rds"))
PRISMA_interim_Batch2 <- readRDS(here("profiles/WGS/230727_PRISMA_Batch2_P3/Results/collated/res_mOTUs.rds"))
PRISMA_interim_Batch3 <- readRDS(here("profiles/WGS/230727_PRISMA_Batch3_P3/Results/collated/res_mOTUs.rds"))
PRISMA_interim_Batch4 <- readRDS(here("profiles/WGS/230804_PRISMA_Batch4_P3/Results/collated/res_mOTUs.rds"))
PRISMA_interim_Batch5 <- readRDS(here("profiles/WGS/230804_PRISMA_Batch5_P3/Results/collated/res_mOTUs.rds"))
PRISMA_modelling <- readRDS(here("profiles/WGS/240813_MG011_PRISMA_NovaSeq_final/Results/collated/res_mOTUs.rds"))

# put everything into a named list
profiles <- list(
    interim_Batch1 = PRISMA_interim_Batch1,
    interim_Batch2 = PRISMA_interim_Batch2,
    interim_Batch3 = PRISMA_interim_Batch3,
    interim_Batch4 = PRISMA_interim_Batch4,
    interim_Batch5 = PRISMA_interim_Batch5,
    modelling = PRISMA_modelling
)

profiles <- map2(profiles, names(profiles), \(x, na) {
    cn <- colnames(x)
    if (str_detect(na, "Batch4") | str_detect(na, "Batch5")) {
        cn <- map_chr(cn, \(x) str_replace(x, ".*lane1", "MG_"))
        #cn <- str_replace(cn, "MG", "MG_")
    } else {
        cn <- map_chr(cn, \(x) str_replace(x, ".*lane1", ""))
        cn <- str_replace(cn, "MG", "MG_")
    }

    colnames(x) <- cn
    return(x)
    }
)

map(profiles, \(x) head(colnames(x))) 

metadata_files <- map(names(profiles), \(x) {
    tmp <- read_tsv(here('profiles/WGS/', str_c(x, '_meta.tsv')))
})

names(metadata_files) <- names(profiles)
metadata_files <- map2(metadata_files, names(metadata_files), \(x, me) {
    x <- x %>%
        mutate(PSN = str_replace(PSN, "Mue", "M")) %>%
        mutate(PSN = str_replace(PSN, "MÜ", "M")) %>%
        mutate(PSN = str_replace(PSN, "ä", "ae")) %>%
        mutate(PSN = str_replace(PSN, "ö", "oe")) %>%
        mutate(PSN = str_replace(PSN, "ü", "ue")) %>%
        mutate(PSN = str_replace(PSN, "Ä", "AE")) %>%
        mutate(PSN = str_replace(PSN, "Ö", "OE")) %>%
        mutate(PSN = str_replace(PSN, "Ü", "UE")) %>%
        mutate(PSN = str_replace(PSN, "NTXM", "NZMU")) %>%
        mutate(PSN = str_replace(PSN, "-0", "-")) %>%
        # Hard code this below - I'm having checks and balanced in place later on...
        mutate(PSN = ifelse(PSN == "JoeFr-NZMUue-5", "JoeFr-NZMU-5", PSN)) %>% 
        mutate(PSN = ifelse(PSN == "SaOs-NZMUue-22", "SaOs-NZMU-22", PSN)) %>%
        mutate(PSN = ifelse(PSN == "DoJE-NZHD-57", "DoJe-NZHD-57", PSN)) %>%
        mutate(batch = me)
    return(x)
})

for (batch in c(
    "interim_Batch1",
    "interim_Batch2",
    "interim_Batch3",
    "interim_Batch4",
    "interim_Batch5")) {
        metadata_files[[batch]] <- metadata_files[[batch]] %>%
            mutate(Visit = as.numeric(str_replace(SampleName, ".*_Visit_", ""))) %>%
            select(-SampleName, -ID) %>%
            relocate(PSN, Visit, Sample_ID)
    }

map(metadata_files, \(x) head(colnames(x)))

## Some sanity checks
pmap(list(profiles, metadata_files, names(metadata_files)), \(x, y, na) {
    print(na)
    print(colnames(x) %in% y$Sample_ID)
    print((y$Sample_ID %in% colnames(x)))
})

profiles <- map2(names(profiles), profiles, \(batch_name, x) x %>%
    as.data.frame() %>%
    rownames_to_column('taxon') %>%
    pivot_longer(-taxon) %>%
    rename(sampleID = name, count = value) %>%
    mutate(batch = batch_name)) %>%
    do.call('rbind', .) %>%
    pivot_wider(id_cols = c(taxon), names_from = c(sampleID, batch), values_from = count, values_fill = 0, names_sep = "___") %>%
    as.data.frame() %>% column_to_rownames('taxon') %>% as.matrix()

meta <- do.call('rbind', metadata_files) %>%
    rename(sampleID = Sample_ID) %>%
    rename(visit = Visit)

# god forgive me...

profiles <- cbind(profiles, rownames(profiles))
colnames(profiles)[dim(profiles)[2]] <- "motu_raw"
profiles <- as.data.frame(profiles)
# profiles$mOTU_ID <- str_split_fixed(rownames(profiles), "_", n = 8)[, 8]
# profiles$species <- str_split_fixed(rownames(profiles), "_", n = 8)[, 7]
# profiles$genus <- str_split_fixed(rownames(profiles), "_", n = 8)[, 6]
# profiles$family <- str_split_fixed(rownames(profiles), "_", n = 8)[, 5]

profiles <- profiles %>%
    as.data.frame() %>%
    pivot_longer(
        -c(
            motu_raw
            )) %>%
    mutate(sampleID = str_split_fixed(name, "___", n = 2)[, 1]) %>%
    mutate(batch = str_split_fixed(name, "___", n = 2)[, 2]) %>%
    group_by(
        motu_raw, 
        sampleID, 
        batch) %>%
    mutate(value = as.numeric(value)) %>%
    summarize(count = sum(value)) %>%
    inner_join(meta, by = c('sampleID', 'batch'))

print("Computing depths...")
depths <- profiles %>%
    group_by(sampleID, batch, visit) %>%
    summarize(totalReadCount = sum(count))

print("Plotting depth histogram...")
depth_histo <- ggplot(depths %>%
    #mutate(visit = as.factor(visit))) +
    mutate(visit = factor(visit, levels = 1:7))) +
    theme_classic() +
    geom_histogram(aes(x = totalReadCount), alpha = 1, position = "identity") +
    geom_vline(xintercept = rarefactionDepthWGS) +
    geom_text(data = depths %>%
        mutate(visit = as.factor(visit)) %>%
        group_by(visit, batch) %>%
        mutate(visit = factor(visit, levels = 1:7)) %>%
        filter(totalReadCount < rarefactionDepthWGS) %>%
        tally(), aes(x = rarefactionDepthWGS + 5000, y = 5, label = n, color = 'red')) +
    facet_grid(visit ~ batch) +
    ylab("Number of\nsamples") +
    xlab("mOTU counts") +
    theme_presentation() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(breaks = c(1, 3, 5, 7, 9, 11))

ggsave(plot = depth_histo, filename = here("plots/KLGPG_221206/depth_histogram_WGS.pdf"), width = 10, height = 6)

highDepthSamples <- depths %>%
    mutate(visit = as.factor(visit)) %>%
    group_by(sampleID, batch) %>%
    filter(totalReadCount >= rarefactionDepthWGS) %>%
    select(sampleID)

print("Keeping samples with reasonably depth...")
print(str_c("Rarefaction depth: ", rarefactionDepthWGS))
#for (B in c("modellingBatchA", "modellingBatchB", "interim")) {
for (B in unique(depths$batch)) {
    n <- profiles %>% ungroup() %>% filter(batch == B) %>% anti_join(highDepthSamples, by = c("sampleID", "batch")) %>% select(sampleID) %>% distinct() %>% nrow()
    N <- profiles %>% ungroup() %>% filter(batch == B) %>% select(sampleID) %>% distinct() %>% nrow()
    f <- n / N
    print(str_c("Removing samples with low depth in batch ", B, ": ", n, ' corresponding to ', round(f, 3), '% of all samples'))
}

profiles <- profiles %>%
    inner_join(highDepthSamples, by = c("sampleID", "batch"))

print("Rarefying...")
set.seed(112312)
profiles <- profiles %>%
    pivot_wider(id_cols = c(sampleID, batch), names_from = motu_raw, values_from = count) %>%
    mutate(tmp = str_c(sampleID, batch, sep = "___")) %>%
    relocate(tmp) %>%
    column_to_rownames("tmp") %>%
    select(-sampleID, -batch) %>%
    rrarefy(sample = rarefactionDepthWGS) %>%
    as.data.frame() %>%
    rownames_to_column('tmp') %>%
    mutate(sampleID = str_split_fixed(tmp, "___", n = 2)[, 1]) %>%
    mutate(batch = str_split_fixed(tmp, "___", n = 2)[, 2]) %>%
    select(-tmp) %>%
    pivot_longer(-c(sampleID, batch)) %>%
    rename(motu_raw = name, count = value)

profiles <- profiles %>%
    group_by(sampleID, batch) %>%
    mutate(relAb = count / sum(count)) %>%
    select(-count) %>%
    filter(motu_raw != "unassigned")

print("Computing pairwise distances and PCOA object...")

profiles <- profiles %>%
    mutate(relAbOrig = relAb) %>%
    mutate(relAb = log10(relAb + pseudoCount)) %>%
    inner_join(meta, by = c("sampleID", "batch")) %>%
    #left_join(fullTax %>% mutate(genus = str_replace(genus, "^g__", "")), by = 'genus')
    identity() # fix later

profiles <- profiles %>%
    group_by(PSN, visit) %>%
    nest() %>%
    group_by(PSN) %>%
    nest() %>%
    mutate(data = map(data, function(x) {
        x <- x %>%
            mutate(visit = as.numeric(as.character(visit))) %>%
            arrange(visit)
        return(x)
    })) %>%
    unnest() %>%
    unnest() %>%
    relocate(sampleID, motu_raw, relAb, relAbOrig, PSN, visit)

pairwiseDistances <- pivot_wider(profiles, id_cols = motu_raw, names_from = c(sampleID, batch), values_from = relAb, names_sep = "___") %>%
    column_to_rownames("motu_raw") %>%
    as.data.frame() %>%
    as.matrix() %>%
    t() %>%
    vegdist(method = "euclidean", k = 2)

pairwiseDistancesIdentityEuclidean <- pivot_wider(profiles, id_cols = motu_raw, names_from = c(sampleID, batch), values_from = relAbOrig, names_sep = "___") %>%
    column_to_rownames("motu_raw") %>%
    as.data.frame() %>%
    as.matrix() %>%
    t() %>%
    vegdist(method = "euclidean", k = 2)

set.seed(1)
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
# model_covariates <- read_excel('/g/scb/zeller/karcher/PRISMA/data/16S_metadata/covariable_columns.xlsx') %>%
model_covariates <- read_excel('/g/scb/zeller/karcher/PRISMA/data/16S_metadata/covariable_columns_v2.xlsx') %>%
    # ATTENTION: This will have to change at a later stage
    filter(C_D_Ratio_Relevance_prio)
tmp <- model_covariates$Column_newName
model_covariates <- model_covariates$Column
names(model_covariates) <- tmp

outcomeInformation <- outcomeInformation %>%
    rename(
        patientID = v65_pat_id,
        visitNumber = v62_visit_number,
        height = v92_height, # height in cm
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
        tac_concentration = v254_Tacrolimus,
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
stopifnot(all(profiles$PSN %in% outcomeInformation$patientID))

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

##########################################################
############## This keeps coming to bite you in the ass...
##########################################################
outcomeInformation <- outcomeInformation[, !str_detect(colnames(outcomeInformation), "v[0-9]+_") | colnames(outcomeInformation) %in% names(model_covariates)]

# Merge dose columns to have only one meaningful one and then calc CD ratio
outcomeInformation <- outcomeInformation %>%
    mutate(fin_tac_dose = pmap_dbl(list(tacDosePrograf, tacDoseEnvarsus, tacDoseAdvagraf, tacDoseModigraf), function(a, b, c, d) {
        tmp <- c(a, b, c, d)
        if (all(is.na(tmp))) {
            return(NA)
        }
        # stopifnot(sum(!is.na(tmp)) == 1)
        return(tmp[!is.na(tmp)])
    })) %>%
    select(-all_of(colnames(.)[str_detect(colnames(.), 'tacDose')])) %>%
    mutate(CD = tac_concentration / fin_tac_dose) %>%
    # height in cm, weight in kg, Haycock formula
    mutate(bsa_haycock = 0.024265 * (height^0.3964) * weight^0.5378) %>%
    # mutate(bsa_duboisdubois = 0.20247 * (height / 100)^0.725 * weight^0.425) %>%
    # mutate(bsa_mosteller = sqrt((height * weight) / 3600)) %>%
    mutate(CD_corrected = tac_concentration / (fin_tac_dose / bsa_haycock))
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

if (FALSE) {

    compare_CD_with_CD_corrected(
        outcomeInformation,
        CD_corrected_for_body_surface_area_midpoint = 1,
        CD_midpoint = 1)

    compare_CD_with_CD_corrected(
        outcomeInformation,
        CD_corrected_for_body_surface_area_midpoint = 'median',
        CD_midpoint = 1)


    compare_CD_with_CD_corrected(
        outcomeInformation,
        CD_corrected_for_body_surface_area_midpoint = 'median',
        CD_midpoint = 'median')

}

outcomeInformation <- outcomeInformation %>%
    # Finally, set 'others' and 'none' to NA
    mutate(ABxSubClass = ifelse(ABxSubClass %in% c("none"), NA, ABxSubClass))

outcomeInformation <- outcomeInformation %>%
    select(-all_of(abxInfo$allAbx))

outcomeInformation <- outcomeInformation %>%
    mutate(postTransplant = visitNumber >= 4)

outcomeInformation <- outcomeInformation %>%
    mutate(across(c(hospitalization, rejection, changeImmunosuppRegimen), \(x) ifelse(is.na(x), FALSE, x)))

outcomeInformation <- outcomeInformation %>%
    mutate(CDbinary = factor(ifelse(CD >= 1, "high", "low"), levels = c('low', 'high'))) %>%
    mutate(CDbinary_corrected = factor(ifelse(CD_corrected >= median(outcomeInformation$CD_corrected, na.rm = TRUE), "high", "low"), levels = c('low', 'high')))

##################################
##################################
############ IMPORTANT ###########
##################################
##################################

# clinicalMetadata summarizes, clinical metadata. i.e. clinical model
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
            mutate(visitPlusOne = c(visit[2:(length(visit))], NA))
        return(x)
    })) %>%
    unnest() %>%
    unnest() %>%
    mutate(visit = factor(visit, levels = 1:7)) %>%
    mutate(visitPlusOne = factor(visitPlusOne, levels = 2:7)) %>%
    # I'm removing data nested data as I'm fiddling with the visit column and I don't want to overwrite things by mistake...
    select(-data)

# For reason that will become clear later on (essentially, in certain situations, I want to predict outcome at T from microbiome + metadata at T-1),
# I here distinguish between outcomeInformation and clinicalMetadata

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
        CD_corrected,
        CDbinary,
        CDbinary_corrected
    ) %>%
    rename(visit = visitNumber) %>%
    group_by(patientID, visit) %>%
    nest() %>%
    group_by(patientID) %>%
    nest() %>%

    mutate(data = map(data, function(x) {
        x <- x %>%
            mutate(visit = as.numeric(as.character(visit))) %>%
            arrange(visit) %>%
            # mutate(visitPlusOne = c(visit[2:(length(visit))], NA))
            # mutate(visitPlusOne = c(NA, visit[1:(length(visit) - 1)]))
            mutate(visitPlusOne = visit - 1)
        return(x)
    })) %>%
    unnest() %>%
    unnest() %>%
    mutate(visit = factor(visit, levels = 1:7)) %>%
    mutate(visitPlusOne = factor(visitPlusOne, levels = 1:6)) %>%
    mutate(CD = ifelse(CD > 5, 5, CD))

# Transform patientID into factor where ordering corresponds to a meaningful ordering
orderDFPatientID <- clinicalMetadata %>%
    rename(PSN = patientID) %>%
    group_by(PSN) %>%
    left_join(outcomeInformation %>% select(patientID, visit, CD, rejection, changeImmunosuppRegimen, hospitalization), by = c("PSN" = 'patientID', "visit" = 'visit')) %>%
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

importantTaxaMotuRaw <- profiles %>%
    mutate(relAb = (10^relAb) - pseudoCount) %>%
    mutate(taxa = motu_raw) %>%
    mutate(taxa = as.character(taxa)) %>%
    group_by(taxa, PSN, visit) %>%
    summarize(relAb = sum(relAb)) %>%
    group_by(taxa) %>%
    summarize(m = mean(relAb > pseudoCount) > 0.2, mm = any(relAb > 0.01)) %>%
    filter(m & mm) %>%
    select(taxa) %>%
    # filter(taxa %in% c("Roseburia", "Coprococcus", "Anaerostipes", "Enterococcus"))
    identity()

tax_levels_prof <- str_split_fixed(profiles$motu_raw, "[|]", n = 8)
colnames(tax_levels_prof) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species", "motu")
profiles <- cbind(profiles, tax_levels_prof) %>%
    relocate(kingdom, phylum, class, order, family, genus, species, motu) %>%
    filter(!str_detect(genus, "incertae")) %>%
    mutate(
        genus = str_replace(genus, ".*gen. ", ""),
        genus = str_replace(genus, "\\[", ""),
        genus = str_replace(genus, "\\]", ""),
        genus = str_replace_all(genus, "/", "_")
    )

profiles_family <- profiles %>%
    group_by(sampleID, PSN, visit, family, batch) %>%
    summarize(relAbOrig = sum(relAbOrig)) %>%
    mutate(relAb = log10(relAbOrig + pseudoCount)) %>%
    mutate(family = str_replace(family, "f__", ""))

profiles_genus <- profiles %>%
    group_by(sampleID, PSN, visit, genus, batch) %>%
    summarize(relAbOrig = sum(relAbOrig)) %>%
    mutate(relAb = log10(relAbOrig + pseudoCount)) %>%
    mutate(genus = str_replace(genus, "g__", "")) %>%
    #mutate(family = str_replace(family, "f__", "")) %>%
    ungroup()

importantTaxaGenus <- profiles_genus %>%
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

pairwiseDistancesGenus <- pivot_wider(profiles_genus, id_cols = genus, names_from = c(sampleID, batch), values_from = relAb, names_sep = "___") %>%
    column_to_rownames("genus") %>%
    as.data.frame() %>%
    as.matrix() %>%
    t() %>%
    vegdist(method = "euclidean", k = 2)

set.seed(1)
pcoaGenus <- cmdscale(pairwiseDistancesGenus) %>%
    as.data.frame() %>%
    rownames_to_column("tmp") %>%
    mutate(sampleID = str_split_fixed(tmp, "___", n = 2)[, 1]) %>%
    mutate(batch = str_split_fixed(tmp, "___", n = 2)[, 2]) %>%
    select(-tmp) %>%
    left_join(meta, by = c("sampleID", 'batch')) %>%
    as_tibble() %>%
    mutate(visit = as.character(visit)) %>%
    mutate(visit = factor(visit, levels = 1:7))

profiles_wgs <- profiles
profiles_wgs_genus <- profiles_genus
profiles_wgs_family <- profiles_family

print(str_c("Saving all objects for downstream analysis to object", obj_path))

dataList <- c(
    'meta', 
    'profiles_wgs', 
    'profiles_wgs_family', 
    "profiles_wgs_genus", 
    "pcoa", 
    "pcoaGenus",
    'pairwiseDistances', 
    'pairwiseDistancesGenus',
    'pairwiseDistancesIdentityEuclidean', 
    "outcomeInformation", 
    "clinicalMetadata", 
    "orderDFPatientID", 
    "importantTaxaMotuRaw",
    "importantTaxaGenus")
save(list = dataList, file = obj_path)
