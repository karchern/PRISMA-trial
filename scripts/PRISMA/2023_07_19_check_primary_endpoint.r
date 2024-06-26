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

# Load data
obj_path <- here('objects/PRISMA.rdata')
load_data(obj_path)

preTransplantProfiles <- profiles %>%
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

# candidateGenera <- c("Enterococcus", "Roseburia", "Coprococcus")
candidateGenera <- unique(preTransplantProfiles$genus)

##############################################################################
#### Primary endpoint prediction: Predict CD at baseline from microbiome ####
##############################################################################

metabCDThreshold <- 1
tpFilterLow <- 4
tpFilterHigh <- 4
allowDifference <- 1

if (abs(tpFilterHigh - tpFilterLow) <= 1) {
    print("Setting allowDifference variable to 0...")
    allowDifference <- 0
}

p <- ggplot(outcomeInformation %>%
    mutate(visitNumber = factor(visitNumber, levels = 1:7, ordered = TRUE)) %>%
    group_by(patientID) %>%
    filter(!is.na(CD)) %>%
    arrange(visitNumber) %>%
    filter(visitNumber >= tpFilterLow) %>%
    filter(visitNumber <= tpFilterHigh) %>%
    mutate(visitNumber = factor(visitNumber, levels = levels(visitNumber)[(tpFilterLow):length(levels(visitNumber))])) %>%
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
    write_tsv(here("results/CD_metabolism_map.tsv")) %>%
    unnest() %>%
    identity()
, aes(x = visitNumber, y = CD)) +
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

cdModelDataSmall <- read_tsv(here("results/CD_metabolism_map.tsv")) %>%
    select(-data) %>%
    inner_join(
        clinicalMetadata %>%
            select(patientID, visit, cyp3a5star3, cyp3a4star22, firstAlbuminMeasurement, ageCategorical, sex, weight, firstHematocritMeasurement) %>%
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


# This doesn't really make much sense anymore, and is also overfitting (see later)
modelDataAll <- list()
res <- list()
resUnadjusted <- list()
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
        formula = cdMetabolism ~ cyp3a5star3 + cyp3a4star22 + firstAlbuminMeasurement + ageCategorical + firstHematocritMeasurement + sex + weight + relAb, family = 'binomial')
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

resTibble <- tibble(genus = names(res), models = res) %>%
    mutate(summary = map(models, summary)) %>%
    mutate(cyp3a5star3_pvalue = map_dbl(summary, \(x) {
        x$coefficients[rownames(x$coefficients) == "cyp3a5star3TRUE", 4]
    })) %>%
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
    left_join(profiles %>% ungroup() %>% select(genus, phylum) %>% distinct(), by = c('genus' = 'genus')) %>%
    relocate(genus, phylum) %>%
    arrange(taxon_pvalue) %>%
    mutate(taxon_estimate = ifelse(taxon_estimate < -5, -5, taxon_estimate)) %>%
    mutate(taxon_estimate = ifelse(taxon_estimate > 5, 5, taxon_estimate))

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
    left_join(profiles %>% ungroup() %>% select(genus, phylum) %>% distinct(), by = c('genus' = 'genus')) %>%
    relocate(genus, phylum) %>%
    arrange(taxon_pvalue) %>%
    mutate(taxon_estimate = ifelse(taxon_estimate < -5, -5, taxon_estimate)) %>%
    mutate(taxon_estimate = ifelse(taxon_estimate > 5, 5, taxon_estimate))

lab_adjusted <- resTibble %>% filter(taxon_pvalue < 0.1)
pAdjusted <- ggplot(data = resTibble) +
    geom_vline(xintercept = 0, linetype = 'dotted') +
    geom_point(aes(x = taxon_estimate, y = -log10(taxon_pvalue)), alpha = 0.5) +
    geom_text_repel(data = lab_adjusted, aes(x = taxon_estimate, y = -log10(taxon_pvalue), label = genus)) +
    theme_presentation() +
    ggtitle("ADJUSTED log. regression model\npredicting CD metabolism\nfrom baseline information") +
    NULL

lab_unadjusted <- resTibbleUnadjusted %>% filter(taxon_pvalue < 0.1)
pUnadjusted <- ggplot(data = resTibbleUnadjusted) +
    geom_vline(xintercept = 0, linetype = 'dotted') +
    geom_point(aes(x = taxon_estimate, y = -log10(taxon_pvalue)), alpha = 0.5) +
    geom_text_repel(data = resTibbleUnadjusted %>% filter(taxon_pvalue < 0.1), aes(x = taxon_estimate, y = -log10(taxon_pvalue), label = genus)) +
    theme_presentation() +
    ggtitle("UNADJUSTED log. regression model\n predicting CD metabolism\nfrom baseline information") +
    NULL

scatter_data <- full_join(
    resTibble %>% select(genus, taxon_pvalue),
    resTibbleUnadjusted %>% select(genus, taxon_pvalue),
    by = 'genus',
    suffix = c("_adjusted", "_unadjusted")
)
scatter_plot <- ggplot(scatter_data) +
    geom_abline(intercept = 0, slope = 1, linetype = 'dotted') +
    geom_point(aes(x = -log10(taxon_pvalue_adjusted), y = -log10(taxon_pvalue_unadjusted)), alpha = 0.5) +
    geom_text_repel(data =
        scatter_data %>%
            inner_join(
                rbind(
                    lab_adjusted %>% select(genus),
                    lab_unadjusted %>% select(genus)
                ) %>%
                    distinct()
            )
    , aes(x = -log10(taxon_pvalue_adjusted), y = -log10(taxon_pvalue_unadjusted), label = genus)) +
    theme_presentation() +
    xlim((c(0, 2))) +
    ylim(c(0, 2)) +
    NULL

ggsave(pAdjusted + pUnadjusted + scatter_plot + plot_layout(guides = 'collect'), filename = here("plots/KLGPG_221206/glm_cd_cyp_tax_profiles_volcano_plots.pdf"), width = 12, height = 5)

(resTibble %>%
    arrange(taxon_pvalue) %>%
    head(50) %>%
    mutate(genus = factor(genus, levels = genus)) %>%
    ggplot(aes(x = genus, y = taxon_pvalue, fill = phylum)) +
    theme_presentation() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    geom_bar(stat = 'identity')) %>%
    ggsave(filename = here("plots/KLGPG_221206/glm_cd_cyp_tax_profiles.pdf"), width = 12, height = 3.5)

plots <- list()
for (g in c("Kopriimonas", "Erysipelatoclostridium", "Enterococcus", "Roseburia", "Coprococcus", "Parasutterella", 'Sutterella')) {
    # ggsave(filename = str_c(here("plots/KLGPG_221206/"), g, ".pdf"), width = 5, height = 4)
    plots[[length(plots) + 1]] <- illustrate_taxon_hit(do.call('rbind', modelDataAll), g) + ggtitle(g) + theme(plot.title = element_text(size = 8, face = "bold"))
}

ggsave(plot = wrap_plots(plots, guides = 'collect', nrow = 2),
    filename = here("plots/KLGPG_221206/cd_metabolism_hits.pdf"), width = 9, height = 5)

###############################################################################
##  train RF models to predict CD bracket based on clinical meta + microbiome
###############################################################################



# ATTENTION: I have to artificially include some noise cause otherwise the logistic model won't fit...
# cdModelDataSmall$cyp3a5star3[c(1)] <- TRUE
# ATTENTION: I impute albumin/hematocrit/weight with the mean
cdModelDataSmall$firstAlbuminMeasurement[is.na(cdModelDataSmall$firstAlbuminMeasurement)] <- mean(cdModelDataSmall$firstAlbuminMeasurement[!is.na(cdModelDataSmall$firstAlbuminMeasurement)])
cdModelDataSmall$weight[is.na(cdModelDataSmall$weight)] <- mean(cdModelDataSmall$weight[!is.na(cdModelDataSmall$weight)])

resamp_n <- 5

rocObjectModelSmallAll <- list()
for (seed in 1:resamp_n) {
    print(str_c("Seed: ", seed))
    ps <- list()
    set.seed(seed)
    for (patientID in cdModelDataSmall$patientID) {
        test <- cdModelDataSmall[cdModelDataSmall$patientID == patientID, ]
        train <- cdModelDataSmall[cdModelDataSmall$patientID != patientID, ]
        cdModelSmall <- randomForest(cdMetabolism ~ cyp3a5star3 + cyp3a4star22 + firstAlbuminMeasurement + ageCategorical + firstHematocritMeasurement + sex + weight, data = train, proximity = TRUE)
        p <- predict(cdModelSmall, test, type = 'prob')[, 1]
        ps[[length(ps) + 1]] <- p
    }
    rocObjectModelSmall <- roc(predictor = unlist(ps), response = as.numeric(cdModelDataSmall$cdMetabolism))
    rocObjectModelSmall
    rocObjectModelSmallAll[[seed]] <- list(rocObjectModelSmall, ps)
}

cdModelDataBig <- cdModelDataSmall %>%
    inner_join(preTransplantProfiles %>%
        filter(genus %in% candidateGenera) %>%
        select(genus, relAb, PSN) %>%
        rename(patientID = PSN) %>%
        pivot_wider(id_cols = patientID, names_from = genus, values_from = relAb))


rocObjectModelBigAll <- list()
for (seed in 1:resamp_n) {
    print(str_c("Seed: ", seed))
    ps <- list()
    set.seed(seed)
    for (patientID in cdModelDataBig$patientID) {
        test <- cdModelDataBig[cdModelDataBig$patientID == patientID, ]
        train <- cdModelDataBig[cdModelDataBig$patientID != patientID, ]
        cdModelBig <- randomForest(formula = as.formula(str_c("cdMetabolism ~ cyp3a5star3 + firstAlbuminMeasurement + ageCategorical + firstHematocritMeasurement + sex + weight + ", str_c(candidateGenera, collapse = " + "))), data = train, proximity = FALSE)
        p <- predict(cdModelBig, test, type = 'prob')[, 1]
        ps[[length(ps) + 1]] <- p
    }
    rocObjectModelBig <- roc(predictor = unlist(ps), response = as.numeric(cdModelDataBig$cdMetabolism))
    rocObjectModelBigAll[[seed]] <- list(rocObjectModelBig, ps)
}

cdModelDataOnlyTax <- cdModelDataSmall %>%
    inner_join(preTransplantProfiles %>%
        filter(genus %in% candidateGenera) %>%
        select(genus, relAb, PSN) %>%
        rename(patientID = PSN) %>%
        pivot_wider(id_cols = patientID, names_from = genus, values_from = relAb))

rocObjectModelOnlyTaxAll <- list()
for (seed in 1:resamp_n) {
    print(str_c("Seed: ", seed))
    ps <- list()
    set.seed(seed)
    for (patientID in cdModelDataOnlyTax$patientID) {
        test <- cdModelDataOnlyTax[cdModelDataOnlyTax$patientID == patientID, ]
        train <- cdModelDataOnlyTax[cdModelDataOnlyTax$patientID != patientID, ]
        cdModelonlyTax <- randomForest(formula = as.formula(str_c("cdMetabolism ~ ", str_c(candidateGenera, collapse = " + "))), data = train, proximity = FALSE)
        p <- predict(cdModelonlyTax, test, type = 'prob')[, 1]
        ps[[length(ps) + 1]] <- p
    }
    rocObjectModelOnlyTax <- roc(predictor = unlist(ps), response = as.numeric(cdModelDataOnlyTax$cdMetabolism))
    rocObjectModelOnlyTaxAll[[seed]] <- list(rocObjectModelOnlyTax, ps)
}

cdModels <- tibble(
    resamp = 1:resamp_n,
    small_roc = map(rocObjectModelSmallAll, \(x) x[[1]]),
    big_roc = map(rocObjectModelBigAll, \(x) x[[1]]),
    onlytax_roc = map(rocObjectModelOnlyTaxAll, \(x) x[[1]]),
) %>%
    pivot_longer(-resamp) %>%
    rename(model_type = name, roc = value) %>%
    mutate(specs = map(roc, \(x) {
        return(data.frame(TPR = x$specificities, FPR = 1 - x$sensitivities))
    })) %>%
    mutate(auc = map_dbl(roc, \(x) x$auc)) %>%
    mutate(group = case_when(
        model_type == "small_roc" ~ "Clinical model",
        model_type == "big_roc" ~ "CM + microbiome",
        model_type == "onlytax_roc" ~ "microbiome"
    )) %>%
    mutate(group = factor(group, levels = rev(c('Clinical model', "microbiome", "CM + microbiome")), ordered = TRUE)) %>%
    arrange(group) %>%
    rename(Features = group) %>%
    group_by(Features) %>%
    nest() %>%
    ungroup() %>%
    mutate(y = seq(0.175, 0.05, length.out = length(levels(Features)))) %>%
    unnest() %>%
    identity()

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
        select(resamp, Features, specs) %>%
        unnest() %>%
        filter(Features == 'Clinical model'), aes(x = FPR, y = TPR, group = interaction(Features, resamp), color = Features), alpha = 0.5) +
    theme_presentation() +
    scale_color_manual(values = colors) +
    xlab("False Positive Rate") +
    ylab("True Positive Rate") +
    geom_text(data = cdModels %>%
        filter(Features == 'Clinical model') %>%
        group_by(Features) %>%
        summarize(label = round(median(auc), 3), y = y[1]), aes(x = 0.4, y = y, label = str_c(Features, ": ", label)), inherit.aes = FALSE, hjust = 0) +
    NULL

ggsave(
    plot = pClinical,
    filename = here("plots/KLGPG_221206/cdMetabolismPredictionOnlyClinical.pdf"), width = 5, height = 3.25)

pC <- ggplot() +
    geom_line(data = cdModels %>%
        select(resamp, Features, specs) %>%
        unnest() %>%
        filter(Features == 'CM + microbiome'), aes(x = FPR, y = TPR, group = interaction(Features, resamp), color = Features), alpha = 0.5) +
    theme_presentation() +
    scale_color_manual(values = colors) +
    xlab("False Positive Rate") +
    ylab("True Positive Rate") +
    geom_text(data = cdModels %>%
        filter(Features == 'CM + microbiome') %>%
        group_by(Features) %>%
        summarize(label = round(median(auc), 3), y = y[1]), aes(x = 0.4, y = y, label = str_c(Features, ": ", label)), inherit.aes = FALSE, hjust = 0) +
    NULL

ggsave(
    plot = pC,
    filename = here("plots/KLGPG_221206/cdMetabolismPredictionClinicalPlusMicrobiome.pdf"), width = 5, height = 3.25)

pC <- ggplot() +
    geom_line(data = cdModels %>%
        select(resamp, Features, specs) %>%
        unnest() %>%
        filter(Features == "microbiome"), aes(x = FPR, y = TPR, group = interaction(Features, resamp), color = Features), alpha = 0.5) +
    theme_presentation() +
    scale_color_manual(values = colors) +
    xlab("False Positive Rate") +
    ylab("True Positive Rate") +
    geom_text(data = cdModels %>%
        filter(Features == "microbiome") %>%
        group_by(Features) %>%
        summarize(label = round(median(auc), 3), y = y[1]), aes(x = 0.4, y = y, label = str_c(Features, ": ", label)), inherit.aes = FALSE, hjust = 0) +
    NULL

ggsave(
    plot = pC,
    filename = here("plots/KLGPG_221206/cdMetabolismPredictionClinicalOnlyMicrobiome.pdf"), width = 5, height = 3.25)

pAll <- ggplot() +
    geom_line(data = cdModels %>%
        select(resamp, Features, specs) %>%
        unnest(), aes(x = FPR, y = TPR, group = interaction(Features, resamp), color = Features), alpha = 0.5) +
    theme_presentation() +
    scale_color_manual(values = colors) +
    xlab("False Positive Rate") +
    ylab("True Positive Rate") +
    geom_text(data = cdModels %>%
        group_by(Features) %>%
        summarize(label = round(median(auc), 3), y = y[1]), aes(x = 0.4, y = y, label = str_c(Features, ": ", label)), inherit.aes = FALSE, hjust = 0) +
    NULL

ggsave(
    plot = pAll,
    filename = here("plots/KLGPG_221206/cdMetabolismPrediction.pdf"), width = 5, height = 3.25)
