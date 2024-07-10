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

resamp_n_model <- 5

# model_type <- "RF"
model_type <- "logreg"

cd_what <- "CDbinary"
# cd_what <- "CDbinary_corrected"

taxonomy_annot <- "ncbi_mapseq"
# taxonomy_annot <- "gtdb_idtaxa"

if (taxonomy_annot == "ncbi_mapseq") {
    ### For NCBI
    candidate_taxa_for_prediction <- c(
        "Tyzzerella",
        "Anaerosporobacter",
        "Coprococcus",
        "Roseburia",
        "Dorea",
        "Faecalibacterium"
        # "Lachnospiraceae"
        # "Leuconostoc" # Super lowly abundant and heavily dependent on rarefaction seed...
    )
} else if (taxonomy_annot == "gtdb_idtaxa") {
    # For GTDB
    candidate_taxa_for_prediction <- c(
        "Dorea_A",
        # "Coprobacter",
        "Bariatricus",
        "Roseburia",
        "Dorea"
    )
}


## NCBI_ENTRY <-> GTDB_ENTRY
# Tyzzerella <-> Faecalimonas, Anaerotignum
# Anaerosporobacter <-> Anaerosporobacter
# Coprococcos <-> Coprococcus_A, Faecalimonas, Batriatricus
# Load data

if (!taxonomy_annot %in% c("ncbi_mapseq", "gtdb_idtaxa")) {
    stop("Unknown taxonomy annotation")
}

obj_path <- here(str_c('objects/PRISMA_', taxonomy_annot, '.rdata'))
load_data(obj_path)

preTransplantProfiles <- profiles %>%
    mutate(genus = str_replace_all(genus, "-", "_")) %>%
    inner_join(data.frame(visit = c(1, 2))) %>%
    group_by(PSN) %>%
    nest() %>%
    mutate(data = map(data, \(x) {
        if (1 %in% x$visit) {
            return(x %>% filter(visit == 1))
        } else {
            return(x %>% filter(visit == 2))
        }
    })) %>%
    unnest() %>%
    # filter(mean(relAb > 0.01) > 0.1) %>%
    inner_join(importantTaxaGenus %>% rename(genus = taxa))

preTransplantProfilesFamily <- profiles_family %>%
    inner_join(data.frame(visit = c(1, 2))) %>%
    group_by(PSN) %>%
    nest() %>%
    mutate(data = map(data, \(x) {
        if (1 %in% x$visit) {
            return(x %>% filter(visit == 1))
        } else {
            return(x %>% filter(visit == 2))
        }
    })) %>%
    unnest()

candidateGenera <- unique(c(unique(preTransplantProfiles$genus), unique((preTransplantProfilesFamily$family))))
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
    # mutate(visit = factor(visit, levels = levels(visit)[(tpFilterLow):length(levels(visit))])) %>%
    mutate(visit = factor(visit, levels = 4:7)) %>%
    nest() %>%
    # mutate(varianceCD = map_dbl(data, \(x) {
    #     return(var(x$CD))
    # })) %>%
    mutate(`cdMetabolism` = map_chr(data, \(x) {
        samples <- dim(x)[1]
        # browser()
        # return(ifelse(as.character(x[[cd_what]]) == "high", "low", "high"))
        return(flipper[[x[[cd_what]]]])
    })) %>%
    identity() %T>%
    # mutate(`cdRatio` = map_dbl(data, \(x) {
    #     if (dim(x)[1] != 1) {
    #         dsaadsadsds
    #     }
    #     return(as.numeric(x$CD))
    # }
    # )) %T>%
    write_tsv(here("results/CD_metabolism_map.tsv")) %>%
    unnest() %>%
    identity()

if (cd_what == "CDbinary") {
    cd_what_what <- "CD"
} else if (cd_what == "CDbinary_corrected") {
    cd_what_what <- "CD_corrected"
} else {
    stop[str_c("Unknown CD variable: ", cd_what)]
}

p <- ggplot(data = data
    , aes_string(x = "visit", y = cd_what_what)) +
    geom_hline(yintercept = 1, linetype = 'dotted') +
    geom_boxplot(outlier.color = NA) +
    {
        if (abs(tpFilterHigh - tpFilterLow) <= 1) {
            geom_jitter(aes(color = `cdMetabolism`), width = 0.05, height = 0, alpha = 0.3)
        } else {
            geom_jitter(, width = 0.05, height = 0, alpha = 0.3)
        }
    } +
    # geom_path(aes(group = patientID, color = sqrt(varianceCD), alpha = sqrt(varianceCD))) +
    geom_path(aes(group = patientID, color = `cdMetabolism`), alpha = 0.5) +
    theme_presentation() +
    # scale_colour_gradient(low = "grey", high = "red") +
    scale_color_manual(values = cdMetabColors
    ) +
    scale_alpha(range = c(0.2, 1)) +
    #        annotate(geom = "text", x = 0.5, y = 1.5, label = "high\nmetab", hjust = 0) +
    #        annotate(geom = "text", x = 0.5, y = 0.5, label = "low\nmetab", hjust = 0) +
    ggtitle(str_c("Regarding cdMetabolism classification:\nAllowing, ",
        allowDifference,
        "sample(s) to disagree with label")) +
    # scale_x_discrete(drop = FALSE) +
    scale_x_discrete_prisma(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    NULL

ggsave(plot = p, filename = str_c(here("plots/KLGPG_221206/CDOverTime_allowDifference_"), allowDifference, ".pdf"), width = 5, height = 5.5)
ggsave(plot = p + scale_x_discrete_prisma(drop = TRUE), filename = str_c(here("plots/KLGPG_221206/CDOverTime_allowDifference_"), allowDifference, "_dropped.pdf"), width = 3, height = 5.5)
# ggsave(plot = p, filename = here("plots/KLGPG_221206/CDOverTime.png"), width = 5, height = 5.5)

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
    theme_presentation() +
    facet_grid(. ~ cyp_genotype) +
    scale_fill_manual(values = cdMetabColors) +
    ylab("Number of patients")) %>%
    ggsave(filename = here("plots/KLGPG_221206/cyp_genotype_CD_metabolism.pdf"))

###############################################################################
## Fit univariate log. regression models (adjusted and unaadjusted) to get an idea of the association of the microbiome with CD
###############################################################################

cdModelDataSmall <- read_tsv(here("results/CD_metabolism_map.tsv")) %>%
    select(-data) %>%
    inner_join(
        clinicalMetadata %>%
            # select(patientID, visit, cyp3a5star3, cyp3a4star22, firstAlbuminMeasurement, ageCategorical, sex, weight, firstHematocritMeasurement) %>%
            select(patientID, visit, cyp3a5star3, cyp3a4star22, firstAlbuminMeasurement, firstHematocritMeasurement, all_of(microbiome_confounders)) %>%
            # for weight
            filter(!is.na(cyp3a5star3)) %>%
            filter(!is.na(cyp3a4star22)) %>%
            filter(visit == 1) %>%
            select(-visit) %>%
            distinct()
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
for (g in unique(preTransplantProfiles$genus)) {
    print(g)
    cdModelData <- cdModelDataSmall %>%
        left_join(preTransplantProfiles %>%
            filter(genus == g) %>%
            select(genus, relAb, PSN) %>%
            rename(patientID = PSN))
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
        # mutate(cyp3a5star3_pvalue = map_dbl(summary, \(x) {
        #     x$coefficients[rownames(x$coefficients) == "cyp3a5star3TRUE", 4]
        # })) %>%
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
    unnest() %>%
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
    # mutate(cyp3a5star3_pvalue = map_dbl(summary, \(x) {
    #     x$coefficients[rownames(x$coefficients) == "cyp3a5star3TRUE", 4]
    # })) %>%
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
    left_join(profiles %>% ungroup() %>% select(genus, family, phylum) %>% distinct(), by = c('genus' = 'genus')) %>%
    relocate(genus, family, phylum) %>%
    arrange(taxon_pvalue) %>%
    mutate(taxon_estimate = ifelse(taxon_estimate < -5, -5, taxon_estimate)) %>%
    mutate(taxon_estimate = ifelse(taxon_estimate > 5, 5, taxon_estimate))

lab_unadjusted <- resTibbleUnadjusted %>% filter(taxon_pvalue < 0.1)
pUnadjusted <- ggplot(data = resTibbleUnadjusted) +
    geom_vline(xintercept = 0, linetype = 'dotted') +
    geom_point(aes(x = taxon_estimate, y = -log10(taxon_pvalue)), alpha = 0.5) +
    # geom_text_repel(data = resTibbleUnadjusted %>% filter(taxon_pvalue < 0.1), aes(x = taxon_estimate, y = -log10(taxon_pvalue), label = genus)) +
    geom_text_repel(data = resTibbleUnadjusted %>%
        # filter(genus %in% candidate_taxa_for_prediction)
        arrange(taxon_pvalue) %>%
        head(10)
    , aes(x = taxon_estimate, y = -log10(taxon_pvalue), label = genus), max.overlaps = Inf) +
    theme_presentation() +
    xlab("Effect size [odds ratio]") +
    ylab("-log10(p-value)") +
    ggtitle("UNADJUSTED log. regression model\n predicting CD metabolism\nfrom baseline information") +
    NULL

# ggsave(pAdjusted + pUnadjusted + scatter_plot + plot_layout(guides = 'collect'), filename = here("plots/KLGPG_221206/glm_cd_cyp_tax_profiles_volcano_plots.pdf"), width = 12, height = 5)
ggsave(pUnadjusted + plot_layout(guides = 'collect'), filename = here("plots/KLGPG_221206/glm_cd_cyp_tax_profiles_volcano_plots.pdf"), width = 4, height = 4.5)

(resTibbleUnadjusted %>%
    arrange(taxon_pvalue) %>%
    head(50) %>%
    mutate(genus = factor(genus, levels = genus)) %>%
    ggplot(aes(x = genus, y = taxon_pvalue, fill = phylum)) +
    theme_presentation() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    geom_bar(stat = 'identity')) %>%
    ggsave(filename = here("plots/KLGPG_221206/glm_cd_cyp_tax_profiles_phylum.pdf"), width = 12, height = 3.5)


top_fam <- resTibbleUnadjusted %>%
    arrange(taxon_pvalue) %>%
    head(50) %>%
    group_by(family) %>%
    tally() %>%
    arrange(desc(n)) %>%
    head(3) %>%
    pull(family)
resTibbleUnadjusted$top_family <- resTibbleUnadjusted$family
resTibbleUnadjusted <- resTibbleUnadjusted %>%
    mutate(top_family = ifelse(!top_family %in% top_fam, "Other", top_family)) %>%
    mutate(top_family = factor(top_family, levels = c(top_fam, "Other")))
(resTibbleUnadjusted %>%
    arrange(taxon_pvalue) %>%
    head(50) %>%
    mutate(genus = factor(genus, levels = genus)) %>%
    ggplot(aes(x = genus, y = taxon_pvalue, fill = top_family)) +
    theme_presentation() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    geom_bar(stat = 'identity')) %>%
    ggsave(filename = here("plots/KLGPG_221206/glm_cd_cyp_tax_profiles_family.pdf"), width = 12, height = 3.5)

plots <- list()
for (g in candidate_taxa_for_prediction) {
    plots[[length(plots) + 1]] <- illustrate_taxon_hit(do.call('rbind', modelDataAll), g, meta, by_batch = FALSE) + ggtitle(g) + theme(plot.title = element_text(size = 8, face = "bold"))
}
for (g in c("Lachnospiraceae")) {
    plots[[length(plots) + 1]] <- illustrate_taxon_hit(cdModelDataSmall %>%
        left_join(preTransplantProfilesFamily %>%
            filter(family == g) %>%
            select(family, relAb, PSN) %>%
            rename(patientID = PSN)), g, meta, by_batch = FALSE, tax_level = "family") + ggtitle(g) + theme(plot.title = element_text(size = 8, face = "bold"))
}

ggsave(plot = wrap_plots(plots, guides = 'collect', nrow = 3),
    filename = here("plots/KLGPG_221206/cd_metabolism_hits.pdf"), width = 6.25, height = 6)

plots <- list()
for (g in candidate_taxa_for_prediction) {
    plots[[length(plots) + 1]] <- illustrate_taxon_hit(do.call('rbind', modelDataAll), g, meta, by_batch = TRUE) + ggtitle(g) + theme(plot.title = element_text(size = 8, face = "bold"))
}
for (g in c("Lachnospiraceae")) {
    plots[[length(plots) + 1]] <- illustrate_taxon_hit(cdModelDataSmall %>%
        left_join(preTransplantProfilesFamily %>%
            filter(family == g) %>%
            select(family, relAb, PSN) %>%
            rename(patientID = PSN)), g, meta, tax_level = "family", by_batch = TRUE) + ggtitle(g) + theme(plot.title = element_text(size = 8, face = "bold"))
}

ggsave(plot = wrap_plots(plots, guides = 'collect', nrow = 3),
    filename = here("plots/KLGPG_221206/cd_metabolism_hits_by_batch.pdf"), width = 8, height = 5)

###############################################################################
##  train RF models to predict CD bracket based on clinical meta + microbiome
###############################################################################

cdModelDataSmall$firstAlbuminMeasurement[is.na(cdModelDataSmall$firstAlbuminMeasurement)] <- mean(cdModelDataSmall$firstAlbuminMeasurement[!is.na(cdModelDataSmall$firstAlbuminMeasurement)])
cdModelDataSmall$weight[is.na(cdModelDataSmall$weight)] <- mean(cdModelDataSmall$weight[!is.na(cdModelDataSmall$weight)])

# rocObjectModelSmallcyp3a5star3 <- get_model_performances(
#     model_data = cdModelDataSmall,
#     # model_feature_string = c("cyp3a5star3", "cyp3a4star22", "firstAlbuminMeasurement", "ageCategorical", "firstHematocritMeasurement", "sex", "weight"),
#     model_feature_string = c("cyp3a5star3"),
#     resamp_n_model = resamp_n_model,
#     microbial_feature_selection_internal = FALSE,
#     # model_type = "logreg")
#     model_type = model_type)

# rocObjectModelSmallcyp3a4star22 <- get_model_performances(
#     model_data = cdModelDataSmall,
#     # model_feature_string = c("cyp3a5star3", "cyp3a4star22", "firstAlbuminMeasurement", "ageCategorical", "firstHematocritMeasurement", "sex", "weight"),
#     model_feature_string = c("cyp3a4star22"),
#     resamp_n_model = resamp_n_model,
#     microbial_feature_selection_internal = FALSE,
#     # model_type = "logreg")
#     model_type = model_type)

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
    # model_feature_string = c("cyp3a5star3", "cyp3a4star22"),
    resamp_n_model = resamp_n_model,
    microbial_feature_selection_internal = FALSE,
    # model_type = "logreg")
    model_type = model_type)

cdModelDataBig <- cdModelDataSmall %>%
    inner_join(preTransplantProfiles %>%
        filter(genus %in% candidateGenera) %>%
        select(genus, relAb, PSN) %>%
        rename(patientID = PSN) %>%
        pivot_wider(id_cols = patientID, names_from = genus, values_from = relAb)) %>%
    left_join(
        preTransplantProfilesFamily %>%
            select(family, relAb, PSN) %>%
            rename(patientID = PSN) %>%
            inner_join(data.frame(family = candidate_taxa_for_prediction)) %>% pivot_wider(id_cols = patientID, names_from = family, values_from = relAb)
    )

rocObjectModelBig <- get_model_performances(
    model_data = cdModelDataBig,
    model_feature_string = c("cyp3a5star3", candidateGenera),
    resamp_n_model = resamp_n_model,
    microbial_feature_selection_internal = candidate_taxa_for_prediction,
    # microbial_feature_selection_internal = TRUE,
    # microbial_feature_selection_internal = FALSE,
    # model_type = "logreg")
    model_type = model_type)

rocObjectModelBigAll <- get_model_performances(
    model_data = cdModelDataBig,
    model_feature_string = c(clinical_covars, candidateGenera),
    resamp_n_model = resamp_n_model,
    microbial_feature_selection_internal = candidate_taxa_for_prediction,
    # microbial_feature_selection_internal = TRUE,
    # microbial_feature_selection_internal = FALSE,
    # model_type = "logreg")
    model_type = model_type)

cdModelDataOnlyTax <- cdModelDataSmall %>%
    inner_join(preTransplantProfiles %>%
        filter(genus %in% candidateGenera) %>%
        select(genus, relAb, PSN) %>%
        rename(patientID = PSN) %>%
        pivot_wider(id_cols = patientID, names_from = genus, values_from = relAb)) %>%
    left_join(
        preTransplantProfilesFamily %>%
            select(family, relAb, PSN) %>%
            rename(patientID = PSN) %>%
            inner_join(data.frame(family = candidate_taxa_for_prediction)) %>% pivot_wider(id_cols = patientID, names_from = family, values_from = relAb)
    )

rocObjectModelOnlyTaxAll <- get_model_performances(
    model_data = cdModelDataOnlyTax,
    model_feature_string = candidateGenera,
    resamp_n_model = resamp_n_model,
    microbial_feature_selection_internal = candidate_taxa_for_prediction,
    # microbial_feature_selection_internal = TRUE,
    # microbial_feature_selection_internal = FALSE,
    # model_type = "logreg")
    model_type = model_type)

cdModels <- tibble(
    resamp = 1:resamp_n_model,
    # small_roc = map(rocObjectModelSmall, \(x) x[[1]]),
    # big_roc = map(rocObjectModelBig, \(x) x[[1]]),
    big_all_roc = map(rocObjectModelBigAll, \(x) x[[1]]),
    small_all_roc = map(rocObjectModelSmallAll, \(x) x[[1]]),
    onlytax_roc = map(rocObjectModelOnlyTaxAll, \(x) x[[1]])
) %>%
    pivot_longer(-resamp) %>%
    rename(model_type = name, roc = value) %>%
    mutate(specs = map(roc, \(x) {
        return(data.frame(TPR = x$specificities, FPR = 1 - x$sensitivities))
    })) %>%
    mutate(auc = map_dbl(roc, \(x) x$auc)) %>%
    mutate(group = case_when(
        # model_type == "small_roc" ~ "cyp genotype",
        model_type == "small_all_roc" ~ "clinical model",
        # model_type == "big_roc" ~ "cyp + microbiome",
        model_type == "big_all_roc" ~ "CM + microbiome",
        model_type == "onlytax_roc" ~ "microbiome"
    )) %>%
    mutate(group = factor(group, levels = rev(c(
        #' cyp genotype',
        # "cyp + microbiome",
        'CM + microbiome',
        'clinical model',
        "microbiome")), ordered = TRUE)) %>%
    arrange(group) %>%
    rename(Features = group) %>%
    group_by(Features) %>%
    nest() %>%
    ungroup() %>%
    mutate(y = seq(0.15, 0.025, length.out = length(levels(Features)))) %>%
    unnest() %>%
    identity()


# colVec <- # Define colors
grey_color <- "#888888" # Grey
blue_color <- "#3498db" # Blue
green_color <- "#2ecc71" # Green
purple_color <- "#9b59b6" # Purple
red_color <- "#e74c3c" # Red
orange_color <- "#F39C12" # Orange
teal_color <- "#1ABC9C" # Teal
pink_color <- "#E84393" # Pink
brown_color <- "#8B4513" # Brown
navy_color <- "#2C3E50" # Navy

# Display the colors
colors <- c(blue_color, red_color, green_color, purple_color, orange_color)
names(colors) <- c(levels(cdModels$Features), "cyp3a5star3", "cyp3a4star22")

pClinical <- ggplot() +
    geom_line(data = cdModels %>%
        select(resamp, Features, specs) %>%
        unnest() %>%
        filter(Features == 'cyp genotype'), aes(x = FPR, y = TPR, group = interaction(Features, resamp), color = Features), alpha = 1) +
    theme_presentation() +
    scale_color_manual(values = colors) +
    xlab("False Positive Rate") +
    ylab("True Positive Rate") +
    geom_text(data = cdModels %>%
        filter(Features == 'cyp genotype') %>%
        group_by(Features) %>%
        summarize(label = round(median(auc), 3), y = y[1]), aes(x = 0.275, y = y, label = str_c(Features, ": ", label)), inherit.aes = FALSE, hjust = 0) +
    NULL

ggsave(
    plot = pClinical,
    filename = here(str_c("plots/KLGPG_221206/cdMetabolismPredictionOnlyClinical", model_type, ".pdf")), width = 5, height = 3.25)

pC <- ggplot() +
    geom_line(data = cdModels %>%
        select(resamp, Features, specs) %>%
        unnest() %>%
        filter(Features == 'cyp + microbiome'), aes(x = FPR, y = TPR, group = interaction(Features, resamp), color = Features), alpha = 1) +
    theme_presentation() +
    scale_color_manual(values = colors) +
    xlab("False Positive Rate") +
    ylab("True Positive Rate") +
    geom_text(data = cdModels %>%
        filter(Features == 'cyp + microbiome') %>%
        group_by(Features) %>%
        summarize(label = round(median(auc), 3), y = y[1]), aes(x = 0.275, y = y, label = str_c(Features, ": ", label)), inherit.aes = FALSE, hjust = 0) +
    NULL

ggsave(
    plot = pC,
    filename = here(str_c("plots/KLGPG_221206/cdMetabolismPredictionClinicalPlusMicrobiome", model_type, ".pdf")), width = 5, height = 3.25)

pC <- ggplot() +
    geom_line(data = cdModels %>%
        select(resamp, Features, specs) %>%
        unnest() %>%
        filter(Features == "microbiome"), aes(x = FPR, y = TPR, group = interaction(Features, resamp), color = Features), alpha = 1) +
    theme_presentation() +
    scale_color_manual(values = colors) +
    xlab("False Positive Rate") +
    ylab("True Positive Rate") +
    geom_text(data = cdModels %>%
        filter(Features == "microbiome") %>%
        group_by(Features) %>%
        summarize(label = round(median(auc), 3), y = y[1]), aes(x = 0.275, y = y, label = str_c(Features, ": ", label)), inherit.aes = FALSE, hjust = 0) +
    NULL

ggsave(
    plot = pC,
    filename = here(str_c("plots/KLGPG_221206/cdMetabolismPredictionClinicalOnlyMicrobiome", model_type, ".pdf")), width = 5, height = 3.25)

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
    # geom_point(data = data.frame(FPR = vals_cyp3a4star22$FPR, TPR = vals_cyp3a4star22$TPR), aes(x = FPR, y = TPR), color = colors[5], size = 4, shape = 4) +
    geom_text(data = cdModels %>%
        group_by(Features) %>%
        summarize(label = round(median(auc), 3), y = y[1]), aes(x = 0.275, y = y, label = str_c(Features, ": ", label)), inherit.aes = FALSE, hjust = 0) +
    NULL

ggsave(
    plot = pAll,
    filename = here(str_c("plots/KLGPG_221206/cdMetabolismPrediction", model_type, ".pdf")), width = 5, height = 3.25)
