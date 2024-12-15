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
library(randomForest)
library(pROC)
library(ggrepel)

# source('/home/karcher/utils/utils.r')
source(here('scripts/utils.r'))

microbiome_confounders <- c(
    "weight",
    "sex",
    "ageCategorical",
    "v502_cakut_dis",
    "v69_Vasculitis_dis",
    "v70_Diabetic_Nephropathy_dis",
    "v71_Glomerulonephritides_dis",
    "v72_Hypertensive_Nephropathy_dis",
    "v73_Nephrolithiasis_dis",
    "v74_Polycystic_Kidney__dis",
    "v75_Coronary_artery_comor",
    "v76_Vesicoureteral_Reflux_dis",
    "v78_Diabetes_comor",
    "v79_Epilepsy_comor",
    "v80_Hypertension_comor",
    "v81_Inflammatory_bowel_comor",
    "v82_Irritable_bowel_comor",
    "v83_Parkinson_Disease_comor",
    "v84_Psoriasis_comor",
    "v85_Rheumatoid_Arthritis_comor",
    "v12a_smoking",
    "v14_alcohol",
    "v15_diet",
    "v400_previous_tx",
    "v501_pretransplant",
    "v66a_renal_prior",
    "v66c_renal_prior_type"
)

# resamp_n_model <- 5
resamp_n_model <- 1 # 1 for debugging,testing. 5 for production.

model_type <- "RF"
# model_type <- "logreg"

cd_type_to_use <- "CDbinary"
# cd_type_to_use <- "CDbinary_corrected"

# TODO: Make this relevant for which profiles are being read around line 100!
tax_and_profiler_choice <- "ncbi_motus"
# tax_and_profiler_choice <- "ncbi_mapseq"

# These are manually selected taxa that are being used to pre-select the features before training the models
# candidate_taxa_for_prediction is only used/relevant when predefine_features == TRUE
predefine_features <- TRUE
# predefine_features <- FALSE

if (tax_and_profiler_choice == "ncbi_mapseq") {
    candidate_taxa_for_prediction <- c(
        "Tyzzerella",
        "Anaerosporobacter",
        "Coprococcus",
        "Roseburia",
        "Dorea",
        "Faecalibacterium"
    )
} else if (tax_and_profiler_choice == "ncbi_motus") {
    candidate_taxa_for_prediction <- c(
        "ref_mOTU_v31_03702",
        "ref_mOTU_v31_03674",
        "ref_mOTU_v31_03668",
        "ref_mOTU_v31_03667",
        "ref_mOTU_v31_03570",
        "ref_mOTU_v31_00719",
        "ref_mOTU_v31_03690",
        "ref_mOTU_v31_00856",
        "ref_mOTU_v31_05137",
        "ref_mOTU_v31_04300"
    )
}

if (!predefine_features) {
    candidate_taxa_for_prediction <- NULL
}

if (!tax_and_profiler_choice %in% c("ncbi_mapseq", "ncbi_motus")) {
    stop("Unknown taxonomy annotation")
}

obj_path <- here(str_c('objects/PRISMA_', tax_and_profiler_choice, '.rdata'))
load_data(obj_path)

if (tax_and_profiler_choice == "ncbi_mapseq") {
    profiles <- profiles
} else if (tax_and_profiler_choice == "ncbi_motus") {
    motus_species_map <- profiles_wgs %>%
        ungroup() %>% 
        select(species, motu) %>%
        distinct()
    profiles <- profiles_wgs %>%
        mutate(genus = motu) # ATTENTION: I'm naming this 'genus' here but this is just for historical reasons. This should have been naemd 'taxa' to be non-confusing.
    importantTaxaGenus <- importantTaxaMotuRaw %>%
        mutate(taxa = str_split_fixed(taxa, '[|]', n = 8)[, 8])
    profiles_family <- profiles_wgs_family
}

preTransplantProfiles <- profiles %>%
    mutate(genus = str_replace_all(genus, "-", "_")) %>%
    inner_join(data.frame(visit = c(1, 2)), by = 'visit') %>%
    group_by(PSN) %>%
    nest() %>%
    mutate(data = map(data, \(x) {
        if (1 %in% x$visit) {
            return(x %>% filter(visit == 1))
        } else {
            return(x %>% filter(visit == 2))
        }
    })) %>%
    unnest(data) %>%
    inner_join(importantTaxaGenus %>% rename(genus = taxa), by = 'genus')

# preTransplantProfilesFamily <- profiles_family %>%
#     inner_join(data.frame(visit = c(1, 2)), by = 'visit') %>%
#     group_by(PSN) %>%
#     nest() %>%
#     mutate(data = map(data, \(x) {
#         if (1 %in% x$visit) {
#             return(x %>% filter(visit == 1))
#         } else {
#             return(x %>% filter(visit == 2))
#         }
#     })) %>%
#     unnest(data) 

abundant_and_prevalent_taxa <- unique(c(unique(preTransplantProfiles$genus)))
abundant_and_prevalent_taxa <- abundant_and_prevalent_taxa[!str_detect(abundant_and_prevalent_taxa, "\\[")]
clinical_covars <- c("cyp3a5star3", "firstAlbuminMeasurement", "ageCategorical", "firstHematocritMeasurement", "sex", "weight")

##############################################################################
#### Primary endpoint prediction: Predict CD at baseline from microbiome ####
##############################################################################

tpFilterLow <- 5
tpFilterHigh <- 5
allowDifference <- 1

flipper <- list(
    low = "high",
    high = "low"
)

if (abs(tpFilterHigh - tpFilterLow) <= 1) {
    print("Setting allowDifference variable to 0...")
    allowDifference <- 0
}

data <- outcomeInformation %>%
    mutate(visit = factor(visit, levels = 1:7, ordered = TRUE)) %>%
    group_by(patientID) %>%
    filter(!is.na(CD)) %>%
    arrange(visit) %>%
    filter(visit >= tpFilterLow) %>%
    filter(visit <= tpFilterHigh) %>%
    mutate(visit = factor(visit, levels = 4:7)) %>%
    nest() %>%
    mutate(`cdMetabolism` = map_chr(data, \(x) {
        samples <- dim(x)[1]
        return(flipper[[x[[cd_type_to_use]]]])
    })) %>%
    identity() %T>%
    write_tsv(here("results/CD_metabolism_map.tsv")) %>%
    unnest(data) %>%
    identity()

if (cd_type_to_use == "CDbinary") {
    cd_type_to_use_what <- "CD"
} else if (cd_type_to_use == "CDbinary_corrected") {
    cd_type_to_use_what <- "CD_corrected"
} else {
    stop[str_c("Unknown CD variable: ", cd_type_to_use)]
}

p <- ggplot(data = data
    , aes_string(x = "visit", y = cd_type_to_use_what)) +
    geom_hline(yintercept = 1, linetype = 'dotted') +
    geom_boxplot(outlier.color = NA) +
    {
        if (abs(tpFilterHigh - tpFilterLow) <= 1) {
            geom_jitter(aes(color = `cdMetabolism`), width = 0.05, height = 0, alpha = 0.3)
        } else {
            geom_jitter(, width = 0.05, height = 0, alpha = 0.3)
        }
    } +
    geom_path(aes(group = patientID, color = `cdMetabolism`), alpha = 0.5) +
    theme_presentation() +
    scale_color_manual(values = cdMetabColors
    ) +
    scale_alpha(range = c(0.2, 1)) +

    ggtitle(str_c("Regarding cdMetabolism classification:\nAllowing, ",
        allowDifference,
        "sample(s) to disagree with label")) +
    scale_x_discrete_prisma(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    NULL

ggsave(plot = p, filename = str_c(here("plots/KLGPG_221206/CDOverTime_allowDifference_"), allowDifference, ".pdf"), width = 5, height = 5.5)
ggsave(plot = p + scale_x_discrete_prisma(drop = TRUE), filename = str_c(here("plots/KLGPG_221206/CDOverTime_allowDifference_"), allowDifference, "_dropped.pdf"), width = 3, height = 5.5)

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
    theme_presentation() +
    ylab("Number of patients")) %>%
    ggsave(filename = here("plots/KLGPG_221206/cyp_genotypes.pdf"), width = 4, height = 2.5)

(read_tsv(here("results/CD_metabolism_map.tsv")) %>%
    select(-data) %>%
    left_join(
        clinicalMetadata %>%
            select(patientID, cyp3a5star3, cyp3a4star22) %>%
            distinct(),
            by = 'patientID'
    ) %>%
    pivot_longer(-c(patientID, `cdMetabolism`)) %>%
    rename(cyp_genotype = name) %>%
    identity() %>%
    select(-patientID) %>%
    relocate(`cdMetabolism`, value) %>%
    table() %>%
    as.data.frame() %>%
    ggplot() + geom_bar(aes(x = value, y = Freq, fill = cdMetabolism), stat = 'identity') +
    theme_presentation() +
    facet_grid(. ~ cyp_genotype) +
    scale_fill_manual(values = cdMetabColors) +
    ylab("Number of patients")) %>%
    ggsave(filename = here("plots/KLGPG_221206/cyp_genotype_CD_metabolism.pdf"))

###############################################################################
## Fit univariate log. regression models (adjusted and unaadjusted) to get an idea of the association of the microbiome with CD
###############################################################################

cdModelDataSmall <- read_tsv(here("results/CD_metabolism_map.tsv"), show_col_types = FALSE) %>%
    select(-data) %>%
    inner_join(
        clinicalMetadata %>%
            select(patientID, visit, cyp3a5star3, cyp3a4star22, firstAlbuminMeasurement, firstHematocritMeasurement, all_of(microbiome_confounders)) %>%
            # for weight
            filter(!is.na(cyp3a5star3)) %>%
            filter(!is.na(cyp3a4star22)) %>%
            filter(visit == 1) %>%
            select(-visit) %>%
            distinct(),
            by = 'patientID'
    ) %>%
    filter(cdMetabolism != 'mixed') %>%
    mutate(cdMetabolism = factor(cdMetabolism, levels = c('low', 'high'))) %>%
    mutate(`CD-ratio` = factor(ifelse(cdMetabolism == 'low', "high", 'low'), levels = c('low', "high"))) %>%
    mutate(sex = as.factor(sex))

# This doesn't really make much sense anymore, and is also overfitting (see later)
modelDataAll <- list()
res <- list()
resUnadjusted <- list()
resAdjusted <- list()
print("Getting single-variable assocations...")
for (g in unique(preTransplantProfiles$genus)) {
    
    cdModelData <- cdModelDataSmall %>%
        left_join(preTransplantProfiles %>%
            filter(genus == g) %>%
            select(genus, relAb, PSN) %>%
            rename(patientID = PSN),
            by = 'patientID')
    resUnadjustedSmall <- list()
    for (covar in c('naive_model', "all_covariates", microbiome_confounders)) {
        if (covar == 'naive_model') {
            cdModel <- glm(data = cdModelData,
                # relAb here is already log10-scaled...
                formula = as.formula(str_c("cdMetabolism ~ relAb")), family = 'binomial')
        } else if (covar == 'all_covariates') {
            cdModel <- glm(data = cdModelData,
                # relAb here is already log10-scaled...
                formula = as.formula(str_c("cdMetabolism ~ relAb + ", str_c(microbiome_confounders, collapse = " + "))), family = 'binomial')
        } else {
            cdModel <- glm(data = cdModelData,
                # relAb here is already log10-scaled...
                formula = as.formula(str_c("cdMetabolism ~ relAb + ", covar)), family = 'binomial')
        }

        resUnadjustedSmall[[length(resUnadjustedSmall) + 1]] <- cdModel
        names(resUnadjustedSmall)[length(resUnadjustedSmall)] <- covar
    }
    resAdjusted[[length(resAdjusted) + 1]] <- tibble(covar = names(resUnadjustedSmall), models = resUnadjustedSmall) %>%
        mutate(summary = map(models, summary)) %>%
        mutate(taxon_pvalue = map_dbl(summary, \(x) {
            x$coefficients[rownames(x$coefficients) == "relAb", 4]
        })) %>%
        mutate(taxon_pvalue = as.numeric(taxon_pvalue))
    names(resAdjusted)[length(resAdjusted)] <- g

    cdModelUnadjusted <- glm(data = cdModelData,
        # relAb here is already log10-scaled...
        formula = cdMetabolism ~ relAb, family = 'binomial')

    res[[length(res) + 1]] <- cdModel
    names(res)[length(res)] <- g



    resUnadjusted[[length(resUnadjusted) + 1]] <- cdModelUnadjusted
    names(resUnadjusted)[length(resUnadjusted)] <- g

    modelDataAll[[length(modelDataAll) + 1]] <- cdModelData
    names(modelDataAll)[length(modelDataAll)] <- g
}

resTibbleAdjusted <- enframe(resAdjusted) %>%
    unnest(value) %>%
    rename(genus = name, covariate = covar) %>%
    select(genus, covariate, taxon_pvalue)

genus_order <- resTibbleAdjusted %>%
    filter(covariate == 'naive_model') %>%
    arrange(taxon_pvalue) %>%
    pull(genus)

resTibbleAdjusted <- resTibbleAdjusted %>%
    mutate(genus = factor(genus, levels = genus_order)) %>%
    mutate(covariate = factor(covariate, levels = rev(
        c(
            'naive_model',
            'weight',
            'sex',
            'ageCategorical',
            'v15_diet',
            "v12a_smoking",
            "v14_alcohol",
            microbiome_confounders[!microbiome_confounders %in% c(
                'naive_model',
                'weight',
                'sex',
                'ageCategorical',
                'v15_diet',
                "v12a_smoking",
                "v14_alcohol")],
            'all_covariates')
    )
    )
    )

p <- ggplot() +
    geom_tile(data = resTibbleAdjusted, aes(x = genus, y = covariate, fill = -log10(taxon_pvalue)), color = 'white') +
    geom_text(data = resTibbleAdjusted %>%
        mutate(label = ifelse(taxon_pvalue < 0.05, "*", "")), aes(x = genus, y = covariate, label = label), color = '#d14481', nudge_y = -0.25, size = 3) +
    theme_publication() +
    scale_fill_viridis_c() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot = p, filename = here("plots/KLGPG_221206/glm_cd_cyp_tax_profiles_heatmap_single_covariate_adjusted.pdf"), width = 11, height = 4)

resTibbleUnadjusted <- tibble(genus = names(resUnadjusted), models = resUnadjusted) %>%
    mutate(summary = map(models, summary)) %>%
    mutate(taxon_pvalue = map(summary, \(x) {
        x$coefficients[rownames(x$coefficients) == "relAb", 4]
    })) %>%
    mutate(taxon_estimate = map(summary, \(x) {
        x$coefficients[rownames(x$coefficients) == "relAb", 1]
    })) %>%
    mutate(taxon_estimate_na = map_lgl(taxon_estimate, \(x) is.na(x) || length(x) == 0)) %>%
    mutate(taxon_pvalue_na = map_lgl(taxon_pvalue, \(x) is.na(x) || length(x) == 0)) %>%
    filter(!taxon_estimate_na) %>%
    filter(!taxon_pvalue_na) %>%
    mutate(taxon_estimate = as.numeric(taxon_estimate)) %>%
    mutate(taxon_pvalue = as.numeric(taxon_pvalue)) %>%
    left_join(profiles %>% 
        ungroup() %>% 
        select(genus, family, phylum) %>%
        mutate(genus = str_replace(genus, "g__", "")) %>%
        mutate(family = str_replace(family, "f__", "")) %>%
        mutate(phylum = str_replace(phylum, "p__", "")) %>%
        distinct(genus, .keep_all = TRUE) %>%
        ungroup() %>%
        distinct() %>%
        identity()
        , by = c('genus' = 'genus')) %>%
    relocate(genus, family, phylum) %>%
    arrange(taxon_pvalue) %>%
    mutate(taxon_estimate = ifelse(taxon_estimate < -5, -5, taxon_estimate)) %>%
    mutate(taxon_estimate = ifelse(taxon_estimate > 5, 5, taxon_estimate)) %>%
    distinct(genus, .keep_all = TRUE) %>%
    arrange(taxon_pvalue)

resTibbleUnadjusted$taxon_pvalue_adjusted_BH <- p.adjust(resTibbleUnadjusted$taxon_pvalue, method = 'BH')

pUnadjusted <- ggplot(data = resTibbleUnadjusted) +
    geom_vline(xintercept = 0, linetype = 'dotted') +
    geom_point(aes(x = taxon_estimate, y = -log10(taxon_pvalue)), alpha = 0.5) +
    geom_text_repel(data = resTibbleUnadjusted %>%
        arrange(taxon_pvalue) %>%
        head(10)
    , aes(x = taxon_estimate, y = -log10(taxon_pvalue), label = genus), max.overlaps = Inf) +
    theme_presentation() +
    xlab("Effect size [odds ratio]") +
    ylab("-log10(p-value)") +
    ggtitle("UNADJUSTED log. regression model\n predicting CD metabolism\nfrom baseline information") +
    NULL

ggsave(pUnadjusted + plot_layout(guides = 'collect'), filename = here("plots/KLGPG_221206/glm_cd_cyp_tax_profiles_volcano_plots.pdf"), width = 7, height = 7)

tmp <- resTibbleUnadjusted %>%
    head(50) %>% 
    group_by(genus) %>% 
    mutate(family = str_c(family, 1:length(family))) %>%
    arrange(taxon_pvalue) %>%
    {
        if(tax_and_profiler_choice == "ncbi_motus") {
            (.) %>% left_join(motus_species_map, by = c('genus' = 'motu'))
        } else {
            (.)
        }
    }
tmp$genus <- factor(tmp$genus, levels = tmp$genus)

truncate_string <- function(string, max_length = 50) {
    # if (nchar(string) > max_length) {
    #     return(substr(string, 1, max_length))
    # } else {
    #     return(string)
    # }
    return(str_replace(str_c(str_split(string, " ")[[1]][1:2], collapse = " "), "s__", ""))
}

(ggplot(data = tmp, aes(x = genus, y = taxon_pvalue, fill = phylum)) +
theme_presentation() +
theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
geom_bar(stat = 'identity') +
{
    if(tax_and_profiler_choice == "ncbi_motus") {
        scale_x_discrete(labels = map_chr(tmp$species, truncate_string))
    } else {
        NULL
    }
} +
theme(
    plot.margin = unit(c(1, 1, 1, 2), "cm")
)) %>%
ggsave(filename = here("plots/KLGPG_221206/glm_cd_cyp_tax_profiles_phylum.pdf"), width = 12, height = 4.75)

plots <- list()
for (g in candidate_taxa_for_prediction) {
    if(tax_and_profiler_choice == "ncbi_motus") {
        g_title <- motus_species_map$species[motus_species_map$motu == g][1]
        g_title <- truncate_string(g_title)
        g_title <- str_replace(g_title, "s__", "")
    } else {
        g_title <- g
    }
    plots[[length(plots) + 1]] <- illustrate_taxon_hit(do.call('rbind', modelDataAll), g, meta, by_batch = FALSE) + ggtitle(g_title) + theme(plot.title = element_text(size = 8, face = "bold"))
}

ggsave(plot = wrap_plots(plots, guides = 'collect', nrow = 3),
    filename = here("plots/KLGPG_221206/cd_metabolism_hits.pdf"), width = 6.25, height = 6)

plots <- list()
for (g in candidate_taxa_for_prediction) {
    if(tax_and_profiler_choice == "ncbi_motus") {
        g_title <- motus_species_map$species[motus_species_map$motu == g][1]
        g_title <- truncate_string(g_title)
        g_title <- str_replace(g_title, "s__", "")
    } else {
        g_title <- g
    }    
    plots[[length(plots) + 1]] <- illustrate_taxon_hit(do.call('rbind', modelDataAll), g, meta, by_batch = TRUE) + ggtitle(g_title) + theme(plot.title = element_text(size = 8, face = "bold"))
}

ggsave(plot = wrap_plots(plots, guides = 'collect', nrow = 3),
    filename = here("plots/KLGPG_221206/cd_metabolism_hits_by_batch.pdf"), width = 8, height = 5)

###############################################################################
##  train RF models to predict CD bracket based on clinical meta + microbiome
###############################################################################

cdModelDataSmall$firstAlbuminMeasurement[is.na(cdModelDataSmall$firstAlbuminMeasurement)] <- mean(cdModelDataSmall$firstAlbuminMeasurement[!is.na(cdModelDataSmall$firstAlbuminMeasurement)])
cdModelDataSmall$weight[is.na(cdModelDataSmall$weight)] <- mean(cdModelDataSmall$weight[!is.na(cdModelDataSmall$weight)])

vals_cyp3a5star3 <- compute_tpr_fpr_from_variable_and_ground_truths(
    ground_truths_boolean = cdModelDataSmall$cdMetabolism == 'high',
    predictions_boolean = cdModelDataSmall$cyp3a5star3
)
vals_cyp3a4star22 <- compute_tpr_fpr_from_variable_and_ground_truths(
    ground_truths_boolean = cdModelDataSmall$cdMetabolism == 'high',
    predictions_boolean = cdModelDataSmall$cyp3a4star22
)

rocObjectModelSmallAll <- get_model_performances(
    model_data = cdModelDataSmall,
    model_feature_string = clinical_covars,
    resamp_n_model = resamp_n_model,
    microbial_feature_selection_internal = NULL, # Only fit on clinical data
    model_type = model_type
    )

cdModelDataBig <- cdModelDataSmall %>%
    inner_join(preTransplantProfiles %>%
        filter(genus %in% abundant_and_prevalent_taxa) %>%
        select(genus, relAb, PSN) %>%
        rename(patientID = PSN) %>%
        pivot_wider(id_cols = patientID, names_from = genus, values_from = relAb)) %>%
    # left_join(
    #     preTransplantProfilesFamily %>%
    #         select(family, relAb, PSN) %>%
    #         rename(patientID = PSN) %>%
    #         inner_join(data.frame(family = candidate_taxa_for_prediction)) %>% pivot_wider(id_cols = patientID, names_from = family, values_from = relAb)
    # )
    identity()

rocObjectModelBigAll <- get_model_performances(
    model_data = cdModelDataBig,
    model_feature_string = c(clinical_covars, abundant_and_prevalent_taxa),
    resamp_n_model = resamp_n_model,
    microbial_feature_selection_internal = candidate_taxa_for_prediction,
    model_type = model_type,
    taxa_to_use = abundant_and_prevalent_taxa)

cdModelDataOnlyTax <- cdModelDataSmall %>%
    inner_join(preTransplantProfiles %>%
        filter(genus %in% abundant_and_prevalent_taxa) %>%
        select(genus, relAb, PSN) %>%
        rename(patientID = PSN) %>%
        pivot_wider(id_cols = patientID, names_from = genus, values_from = relAb)) %>%
    # left_join(
    #     preTransplantProfilesFamily %>%
    #         select(family, relAb, PSN) %>%
    #         rename(patientID = PSN) %>%
    #         inner_join(data.frame(family = candidate_taxa_for_prediction)) %>% pivot_wider(id_cols = patientID, names_from = family, values_from = relAb)
    # )
    identity()

rocObjectModelOnlyTaxAll <- get_model_performances(
    model_data = cdModelDataOnlyTax,
    model_feature_string = abundant_and_prevalent_taxa,
    resamp_n_model = resamp_n_model,
    microbial_feature_selection_internal = candidate_taxa_for_prediction,
    model_type = model_type,
    taxa_to_use = abundant_and_prevalent_taxa)

cdModels <- tibble(
    resamp = 1:resamp_n_model,
    cm_and_microbiome_roc = map(rocObjectModelBigAll, \(x) x[[1]]),
    clinical_model_roc = map(rocObjectModelSmallAll, \(x) x[[1]]),
    microbiome_roc = map(rocObjectModelOnlyTaxAll, \(x) x[[1]])
) %>%
    pivot_longer(-resamp) %>%
    rename(model_type = name, roc = value) %>%
    mutate(specs = map(roc, \(x) {
        return(data.frame(TPR = x$specificities, FPR = 1 - x$sensitivities))
    })) %>%
    mutate(auc = map_dbl(roc, \(x) x$auc)) %>%
    mutate(group = case_when(
        model_type == "clinical_model_roc" ~ "clinical model",
        model_type == "cm_and_microbiome_roc" ~ "CM + microbiome",
        model_type == "microbiome_roc" ~ "microbiome"
    )) %>%
    mutate(group = factor(group, levels = rev(c(
        'CM + microbiome',
        'clinical model',
        "microbiome")), ordered = TRUE)) %>%
    arrange(group) %>%
    rename(Features = group) %>%
    group_by(Features) %>%
    nest() %>%
    ungroup() %>%
    mutate(y = seq(0.15, 0.025, length.out = length(levels(Features)))) %>%
    unnest(data) %>%
    identity()

blue_color <- "#3498db" 
red_color <- "#e74c3c"
green_color <- "#2ecc71"
purple_color <- "#9b59b6" 
orange_color <- "#F39C12"

# Display the colors
colors <- c(blue_color, red_color, green_color, purple_color, orange_color)
names(colors) <- c(levels(cdModels$Features), "cyp3a5star3", "cyp3a4star22")

pAll <- ggplot() +
    geom_line(data = cdModels %>%
        select(resamp, Features, specs) %>%
        unnest(), aes(x = FPR, y = TPR, group = interaction(Features, resamp), color = Features), alpha = 1) +
    theme_presentation() +
    scale_color_manual(values = colors, breaks = names(colors)) +
    xlab("False Positive Rate") +
    ylab("True Positive Rate") +
    geom_point(data = rbind(
        data.frame(
            FPR = vals_cyp3a5star3$FPR,
            TPR = vals_cyp3a5star3$TPR,
            Features = "cyp3a5star3"
        ),
        data.frame(
            FPR = vals_cyp3a4star22$FPR,
            TPR = vals_cyp3a4star22$TPR,
            Features = "cyp3a4star22"
    )), aes(x = FPR, y = TPR, color = Features), size = 4, shape = 4) +
    geom_text(data = cdModels %>%
        group_by(Features) %>%
        summarize(label = round(median(auc), 3), y = y[1]), aes(x = 0.275, y = y, label = str_c(Features, ": ", label)), inherit.aes = FALSE, hjust = 0) +
    NULL

ggsave(
    plot = pAll,
    filename = here(str_c("plots/KLGPG_221206/cdMetabolismPrediction", model_type, ".pdf")), width = 5, height = 3.25)
