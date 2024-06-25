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
# source('/home/karcher/utils/utils.r')
source(here('scripts/utils.r'))

# Load data
obj_path <- here('objects/PRISMA.rdata')
load_data(obj_path)


########################################################################################
# LMMs to find associations between taxa and (primary and secondary) outcomes over time!
########################################################################################

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
            inner_join(importantTaxaGenus)
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
            # join outcome metadata (will always remain at T)
            left_join(outcomeInformation, by = c("PSN" = 'patientID', "visit" = 'visitNumber')) %>%
            mutate(visit = factor(visit, levels = 1:7))
        # for (outcomeString in c("anyComplication", "rejection", "changeImmunosuppRegimen", "hospitalization", "CD")) {
        # for (outcomeString in c("anyComplication", "rejection", "hospitalization", "CD", "CDbinary")) {
        for (outcomeString in c("CD", "CDbinary")) {
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
                lmms[[length(lmms) + 1]] <- list(lmm, tmp2)
                names(lmms)[length(lmms)] <- str_c(outcomeString, taxon, visitType, sep = "___", collapse = "___")
            }
        }
    }
}

for (lmmType in c("lmms", "lmmsAdjusted")) {
    # for (compType in c('primary', "secondary")) {
    for (compType in c('primary')) {
        lmmsTibble <- tibble(outcomeTaxa = names(lmms), lmms_list = lmms, lmmsAdjusted = lmmsAdjusted) %>%
            mutate(lmms = map(lmms_list, \(x) x[[1]])) %>%
            mutate(data = map(lmms_list, \(x) x[[2]])) %>%
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
            select(taxa, lmms, data, outcome, pValTaxon, pValTaxonAdjusted, effectSizeTaxon, visitType) %>%
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
            theme_presentation() +
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
            theme_presentation() +
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
            ggsave(plot = lmmPhylumBar + lmmHeatmap + plot_layout(widths = c(1, 10), guides = 'collect'), filename = str_c(here("plots/KLGPG_221206/LmmHeatmap_"), lmmType, "__", compType, ".pdf"), width = 7.5, height = 7)
        } else if (compType == "secondary") {
            ggsave(plot = lmmPhylumBar + lmmHeatmap + plot_layout(widths = c(1, 10), guides = 'collect'), filename = str_c(here("plots/KLGPG_221206/LmmHeatmap_"), lmmType, "__", compType, ".pdf"), width = 8, height = 7)
        }
    }
}

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

### CAREFUL: THis code is non-functional in that you cannot distinguish between visit and visitMinusOne here...
plots <- list()
for (g in list(
    c("Roseburia", "CDbinary"),
    c("Anaerostipes", "CDbinary"),
    c("Tyzzerella", "CDbinary"),
    c("Bifidobacterium", "CDbinary"),
    c("Eubacterium", "CDbinary")
)) {
    taxon <- g[1]
    complication <- g[2]
    dir.create(str_c(here("plots/KLGPG_221206/complication_LMMs_results/"), complication), showWarnings = FALSE)
    yLabel <- ifelse(taxon != "richness", "Relative Abundance (log10)", "richness (log10)")
    if (complication == "CD") {
        plots[[length(plots) + 1]] <- compareTaxAssocsScatter(
            taxon = taxon,
            outcomeMeasure = complication,
            lmmObject = lmmsTibble %>% filter(taxa == taxon) %>% filter(outcome == "CD") %>% pull(lmms) %>% magrittr::extract2(1)) + ylab("CD ratio")
    } else {
        plots[[length(plots) + 1]] <- compareTaxAssocsQuantilePlots(taxon = taxon, complication = complication) + ylab(yLabel) #            ggsave(plot = compareTaxAssocsQuantilePlots(taxon = taxon, complication = complication) + ylab(yLabel), filename = str_c(here("plots/KLGPG_221206/complication_LMMs_results/',complication, ")/", taxon, "_LMM_quantilePlots.pdf"), width = 5, height = 5)
    }
}

ggsave(plot = wrap_plots(plots, ncol = 2, nrow = 3, guides = 'collect'), filename = here("plots/KLGPG_221206/complicationAssocsSeleceted.pdf"), width = 8, height = 10)
