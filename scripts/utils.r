# Define some convenience functions

time_point_colors <- c(
    "#3a7bfc", # light blue
    "#3a7bfc", # slightly less intense blue
    "#3a7bfc", # original blue
    "#ff3636",
    "#f87676",
    "#f8aeae",
    "#fcebeb")

labelLink <- c(
    "1" = "pre-transplant (1)",
    "2" = "pre-transplant (2)",
    "3" = "post-immunosuppression\npre-transplant",
    "4" = "week 1 post-transplant",
    "5" = "week 4 post-transplant",
    "6" = "Month 3 post-transplant",
    "7" = "Month 6 post-transplant"
)

cdMetabColors <- c(
    "high" = "#f78480",
    "low" = "#6cd076",
    'mixed' = "#6c6c6c"
)

cdRatioColors <- c(
    "low" = "#f78480",
    "high" = "#6cd076",
    'mixed' = "#6c6c6c"
)

quantilesP <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

batch_colors <- c(
    interim = "#85C1E9",
    modellingBatchA = "#F88379",
    modellingBatchB = "#A9DFBF")

rarefactionDepth <- 1E4
rarefactionDepthWGS <- 5E3
pseudoCount <- 1E-4

load_data <- function(obj_path) {
    print(str_c("Loading data from ", obj_path))

    # Load the data into the global environment
    loaded_objects <- load(obj_path, envir = .GlobalEnv)
    cat("Loaded following objects\n########################")
    for (obj in loaded_objects) {
        print(obj)
    }
}

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
                select(patientID, visit, rejection, changeImmunosuppRegimen, hospitalization) %>%
                mutate(visit = as.factor(visit)),
            by = c('PSN' = 'patientID', "visit" = 'visit'))

    dataB$PSN <- factor(dataB$PSN, levels = orderDFPatientID$PSN)

    dataB$visit <- factor(dataB$visit, levels = 1:7)

    dataB$visit <- factor(dataB$visit, levels = 1:7)

    pObj <- pObj +
        geom_bar(data = dataB,
            aes(x = PSN, y = relAb, fill = taxa), position = 'stack', stat = 'identity') +
        geom_text(data = outcomeInformation %>%
            rename(PSN = patientID, visit = visit) %>%
            # select(PSN, visit, rejection) %>%
            mutate(rejection = ifelse(rejection, "R", "")) %>%
            mutate(PSN = factor(PSN, levels = orderDFPatientID$PSN)),
        aes(x = PSN, y = 1.525, label = rejection), size = 2.5) +
        geom_text(data = outcomeInformation %>%
            rename(PSN = patientID, visit = visit) %>%
            # select(PSN, visit, changeImmunosuppRegimen) %>%
            mutate(changeImmunosuppRegimen = ifelse(changeImmunosuppRegimen, "C", "")) %>%
            mutate(PSN = factor(PSN, levels = orderDFPatientID$PSN)),
        aes(x = PSN, y = 1.375, label = changeImmunosuppRegimen), size = 2.5) +
        geom_text(data = outcomeInformation %>%
            rename(PSN = patientID, visit = visit) %>%
            # select(PSN, visit, changeImmunosuppRegimen) %>%
            mutate(hospitalization = ifelse(hospitalization, "H", "")) %>%
            mutate(PSN = factor(PSN, levels = orderDFPatientID$PSN)),
        aes(x = PSN, y = 1.225, label = hospitalization), size = 2.5) +
        # geom_point(data = dataB %>% ungroup() %>% filter(anyABx) %>% select(PSN, ABxSubClass) %>% distinct()  ,
        # aes(x = PSN, y = 1.1, color = ABxSubClass)) +
        theme_presentation() +
        # facet_grid(anyComplicationEver2 ~ visit)
        facet_grid(visit ~ ., scales = "fixed") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
        scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1.725)) +
        # ylim(c(0, 1.725)) +
        ylab("relative bacterial abundance")

    return(pObj)
}


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

.f_resolve_taxonomy <- function(collated_mat, taxLevel = "genus") {
    # resolves taxonomy for collated vknight results of MAPseq and mOTUs3.
    # assumes full taxonomic annotation in row names, separated by a pipe ("|")

    taxLevel_vec <- c("kingdom", "phylum", "class", "order", "family", "genus", "species", "mOTU")
    stopifnot("taxLevel must be either: kingdom,phylum,class,order,family,genus,species,mOTU" = taxLevel %in% taxLevel_vec)

    # Split by pipe
    # split taxonomy and assign taxonomy names
    taxSeparated_df <-
        suppressWarnings(collated_mat %>%
            as_tibble(rownames = "tax") %>%
            separate(tax, into = taxLevel_vec, sep = "\\|")) %>%
        select(-c(taxLevel_vec[taxLevel_vec != taxLevel]))

    # Group by the selected tax level and summarise counts
    # convert "NA" to "not_resolved" since they represent bacterial (and archaeal reads) that are not resolved at the selected tax level
    taxLevel_mat <-
        taxSeparated_df %>%
        gather(-taxLevel, key = "Sample_ID", value = "count") %>%
        group_by(!!as.symbol(taxLevel), Sample_ID) %>%
        summarise(count = sum(count)) %>%
        pivot_wider(names_from = Sample_ID, values_from = count, values_fill = 0) %>%
        mutate(!!as.symbol(taxLevel) := case_when(is.na(!!as.symbol(taxLevel)) ~ "not_resolved",
            TRUE ~ !!as.symbol(taxLevel)),
        !!as.symbol(taxLevel) := str_remove(!!as.symbol(taxLevel), pattern = "^[a-z]__")) %>%
        column_to_rownames(taxLevel) %>%
        as.matrix()

    return(taxLevel_mat)

}

scale_x_discrete_prisma <- function(labelMap = labelLink, how = 'discrete', ...) {
    if (how == 'discrete') {
        scale_x_discrete(labels = labelMap, ...)
    } else if (how == "continuous") {
        scale_x_continuous(labels = labelMap, ...)
    }
}

make_first_and_third_letter_uppercase <- function(str) {
    if (length(str_split(str, "")[[1]]) == 4) {
        return(paste0(toupper(substr(str, 1, 1)), substr(str, 2, 2), toupper(substr(str, 3, 3)), substr(str, 4, 4)))
    } else if (length(str_split(str, "")[[1]]) == 5) {
        if (str_detect(substr(str, 1, 3), "ae") || str_detect(substr(str, 1, 3), "oe") || str_detect(substr(str, 1, 3), "ue")) {
            return(paste0(toupper(substr(str, 1, 1)), substr(str, 2, 3), toupper(substr(str, 4, 4)), substr(str, 5, 5)))
        } else if (str_detect(substr(str, 3, 5), "ae") || str_detect(substr(str, 3, 5), "oe") || str_detect(substr(str, 3, 5), "ue")) {
            return(paste0(toupper(substr(str, 1, 1)), substr(str, 2, 2), toupper(substr(str, 3, 3)), substr(str, 4, 5)))
        } else {
            dasdafsadsdfasfsfdafsasf
        }

    } else {
        dasdassdaad
    }
}

clean_patient_clinical_metadata <- function(df_all, how = "normal") {

    print(how)
    if (how != 'normal') {
        # I have to clean the patient IDs...
        df_all$v65_pat_id <- map(df_all$v65_pat_id, \(x) {
            a <- str_split(x, "-", n = 3)[[1]][1]
            a <- make_first_and_third_letter_uppercase(a)
            b <- str_split(x, "-", n = 3)[[1]][2]
            b <- toupper(b)
            c <- str_split(x, "-", n = 3)[[1]][3]
            return(str_c(a, b, c, sep = "-"))
        })
        # Also clean the god damn fucking date to be cosnistent with the previous format... shit
        df_all$v13_dob <- map_chr(df_all$v13_dob, \(x) {
            parts <- str_split(x, '[.]')[[1]]
            year <- parts[3]
            month <- parts[2]
            day <- parts[1]
            # This is super hacky but works... the oldest person was born 1948 so this way we can separate the years
            year_prefix <- ifelse(year > 0 && year < 47, "20", "19")
            return(str_c(str_c(year_prefix, year), month, day, sep = '-'))
        })
    }

    #### CLEAN PATIENT IDs####
    # Replace umlauts
    df_all$v65_pat_id <- gsub("ä", "ae", df_all$v65_pat_id)
    df_all$v65_pat_id <- gsub("ö", "oe", df_all$v65_pat_id)
    df_all$v65_pat_id <- gsub("ü", "ue", df_all$v65_pat_id)
    df_all$v65_pat_id <- gsub("Ä", "AE", df_all$v65_pat_id)
    df_all$v65_pat_id <- gsub("Ö", "OE", df_all$v65_pat_id)
    df_all$v65_pat_id <- gsub("Ü", "UE", df_all$v65_pat_id)

    # #Transform to lower characters and trim potential whitespace
    # df_all$v65_pat_id <- tolower(df_all$v65_pat_id)
    # df_all$v65_pat_id <- trimws(df_all$v65_pat_id,which = "both")

    # Homogenize ID lenght
    df_all$v65_pat_id <- gsub("NTXMUE", "NZMU", df_all$v65_pat_id)
    # df_all$v65_pat_id <- gsub("-(\\d)$", "-0\\1", df_all$v65_pat_id)
    df_all$v65_pat_id <- gsub("-0", "-", df_all$v65_pat_id)
    return(df_all)
}


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


compareTaxAssocsQuantilePlots <- function(taxon = NULL, complication = "anyComplication", simpleLegend = TRUE, plot_kind = 'boxplot') {
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

    # browser()
    p <- ggplot(tmp2, aes(x = as.factor(anyABx), y = relAb)) +
        geom_boxplot(aes_string(x = "anyABx", y = "relAb", fill = complication))
    # geom_quantileplot(aes(fill = CDbinary), quantilesP = quantilesP) +
    # scale_fill_quantile(cdMetabColors, quantilesP) +
    theme_presentation()
    # browser()
    # # Remove NAs
    # tmp2 <- tmp2 %>%
    #     filter(!if_any(all_of(c('anyABx', complication, "relAb")), is.na))
    # p <- get_quantile_plot(tmp2, "anyABx", complication, "relAb", plotGroup = NA, expectedNumLevels = 2)
    # return(p)
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

compareTaxAssocsScatter <- function(taxon = NULL, outcomeMeasure = "CD", lmmObject) {
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
    p <- ggplot(data = tmp2, aes(x = relAb, y = outcome)) +
        geom_point(alpha = 0.2) +
        theme_presentation() +
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

illustrate_taxon_hit <- function(modelData = NULL, taxon = NULL, meta = NULL, by_batch = FALSE, tax_level = "genus") {
    modelData <- modelData %>%
        filter(.data[[tax_level]] == taxon) %>%
        rename(`Tacrolimus\nmetabolism` = cdMetabolism) %>%
        mutate(relAb = (10^(relAb) * 100)) %>%
        inner_join(
            meta %>%
                select(PSN, batch) %>%
                distinct(),
            by = c('patientID' = "PSN")) %>%
        mutate(batch = ifelse(batch == "interim", "interim", "rest\nmodel. cohort"))
    p <- ggplot() +
        geom_boxplot(
            data = modelData,
            aes(x = `Tacrolimus\nmetabolism`, y = relAb, fill = `Tacrolimus\nmetabolism`), outlier.color = NA) +
        geom_jitter(
            data = modelData,
            aes(x = `Tacrolimus\nmetabolism`, y = relAb, fill = `Tacrolimus\nmetabolism`), position = position_jitter(), alpha = 0.3) +
        theme_publication() +
        theme(plot.title = element_text(size = 14, face = "bold")) +
        ylab("Bacterial\nrelative abundance [%]") +
        # scale_fill_manual(values = c('low' = "#4a5dca", "high" = "#d43e3e")) +
        scale_fill_manual(values = cdMetabColors) +
        scale_y_continuous(trans = 'log10', limits = c(0.005, max(modelData$relAb) * 1.05)) +
        {
            if (by_batch) {
                facet_grid(~batch)
            } else {
                NULL
            }
        } +
        NULL
    return(p)
}


diagnose_time_shift <- function(df, PSN = NULL) {
    df <- ungroup(df)
    if (is.null(PSN)) {
        set.seed(13213)
        PS <- sample(df$PSN, 1)
        print(PS)
    }
    return(df %>%
        ungroup() %>%
        select(PSN, visit, anyABx, anyComplication) %>%
        rename(`co-medication` = anyABx, complication = anyComplication) %>%
        distinct() %>%
        mutate(visit = as.numeric(as.character(visit))) %>%
        filter(PSN == PS) %>%
        select(-PSN) %>%
        pivot_longer(-visit) %>%
        ggplot() +
        geom_tile(aes(x = visit, y = name, fill = value)) +
        # ylim(c(1,7)) +
        # scale_y_continuous(breaks = 1:7) +
        theme_presentation() +
        ylab("Visit") +
        ggtitle(PS) +
        scale_x_continuous(breaks = 1:7) +
        scale_fill_manual(values = c("TRUE" = "darkgreen", "FALSE" = "darkblue")) +
        NULL)
}

compute_tpr_fpr_from_variable_and_ground_truths <- function(ground_truths_boolean, predictions_boolean) {
    stopifnot(is.logical(ground_truths_boolean) && is.logical(predictions_boolean) && length(ground_truths_boolean) == length(predictions_boolean))
    # Create a confusion matrix
    cm <- table(Predicted = factor(predictions_boolean, levels = c("FALSE", "TRUE")), Actual = factor(ground_truths_boolean, levels = c("FALSE", "TRUE")))
    # Compute TPR and FPR
    TPR <- cm["TRUE", "TRUE"] / (cm["TRUE", "TRUE"] + cm["FALSE", "TRUE"])
    FPR <- cm["TRUE", "FALSE"] / (cm["TRUE", "FALSE"] + cm["FALSE", "FALSE"])

    # Return a list with TPR and FPR
    return(list(TPR = TPR, FPR = FPR))
}

get_wilcox_results_for_internal_filtering <- function(
    tab,
    ground_truth
    ) {
    wilcox_test_results <- list()
    ground_truth_true <- ground_truth == 'low'
    ground_truth_false <- ground_truth == 'high'
    for (genus in colnames(tab)) {
        wilcox_test_results[[genus]] <- wilcox.test(
            x = tab[[genus]][ground_truth_true],
            y = tab[[genus]][ground_truth_false],
        )
    }
    return(wilcox_test_results)
}

get_model_performances <- function(
    model_data,
    model_feature_string = None,
    resamp_n_model = 1,
    microbial_feature_selection_internal = TRUE,
    top_microbial_features_if_microbial_feature_selection_internal = 10,
    model_type = "RF"
    ) {
    model_feature_string_original <- model_feature_string
    rocObjectsAll <- list()
    for (seed in 1:resamp_n_model) {
        print(str_c("Seed: ", seed))
        ps <- list()
        set.seed(seed)
        for (patientID in model_data$patientID) {
            model_feature_string <- model_feature_string_original
            model_feature_string_non_microbial <- model_feature_string_original[!model_feature_string_original %in% candidateGenera]
            test <- model_data[model_data$patientID == patientID, ]
            train <- model_data[model_data$patientID != patientID, ]
            if (!is_logical(microbial_feature_selection_internal) || microbial_feature_selection_internal) {
                all_microbial_features <- c(candidateGenera)
                train_only_microbial <- train[, colnames(train) %in% all_microbial_features]
                train_rest <- train[, !colnames(train) %in% all_microbial_features]
                if (is.logical(microbial_feature_selection_internal)) {
                    print("Running training-set internal feature slection")
                    wilcox_test_results <- get_wilcox_results_for_internal_filtering(
                        train_only_microbial,
                        train$cdMetabolism)
                    top_microbial_features <- enframe(wilcox_test_results) %>%
                        rename(genus = name) %>%
                        mutate(p_val = map_dbl(value, \(x) x$p.value)) %>%
                        arrange(p_val) %>%
                        head(top_microbial_features_if_microbial_feature_selection_internal) %>%
                        select(genus)
                } else {
                    print("Taking predifined top features")
                    top_microbial_features <- data.frame(genus = microbial_feature_selection_internal)
                }
                # Caution: For testing only, since overfitting
                # top_microbial_features <- data.frame(genus = c("Coprococcus"))
                if (!all(top_microbial_features$genus %in% colnames(train_only_microbial))) {
                    print(top_microbial_features$genus[!top_microbial_features$genus %in% colnames(train_only_microbial)])
                    stop("Not all top microbial features you supplied are in the training data.")
                }
                train <- cbind(train_rest, train_only_microbial[, colnames(train_only_microbial) %in% top_microbial_features$genus])
                model_feature_string <- c(model_feature_string_non_microbial, model_feature_string_original[model_feature_string_original %in% top_microbial_features$genus])
            }

            input_formula <- as.formula(str_c("cdMetabolism ~ ", str_c(model_feature_string, collapse = " + ")))
            print(input_formula)
            if (model_type == "RF") {
                cdModel <- randomForest(formula = input_formula, data = train, proximity = TRUE)
            } else if (model_type == "logreg") {
                cdModel <- glm(formula = input_formula, data = train, family = 'binomial')
            } else {
                asdadsd
            }
            if (model_type == "RF") {
                p <- predict(cdModel, test, type = 'prob')[, 1]
            } else if (model_type == "logreg") {
                p <- predict(cdModel, test, type = 'response')
            } else {
                asdadsds
            }

            ps[[length(ps) + 1]] <- p
        }
        rocObject <- roc(predictor = unlist(ps), response = as.numeric(model_data$cdMetabolism))
        rocObject
        rocObjectsAll[[seed]] <- list(rocObject, ps)
    }
    return(rocObjectsAll)
}

compare_CD_with_CD_corrected <- function(
    outcomeInformation = NULL,
    bsa_methods = c(
        "bsa_haycock",
        "bsa_duboisdubois",
        "bsa_mosteller"
    ),
    CD_corrected_for_body_surface_area_midpoint = NULL,
    CD_midpoint = NULL) {

    if (is.numeric(CD_corrected_for_body_surface_area_midpoint)) {

    } else if (is.character((CD_corrected_for_body_surface_area_midpoint)) && CD_corrected_for_body_surface_area_midpoint == 'median') {
        CD_corrected_for_body_surface_area_midpoint <- median(outcomeInformation$CD_corrected, na.rm = TRUE)
    } else {
        stop('CD_corrected_for_body_surface_area_midpoint should be either numeric or character equalling median')
    }

    if (is.numeric(CD_midpoint)) {

    } else if (is.character((CD_midpoint)) && CD_midpoint == 'median') {
        CD_midpoint <- median(outcomeInformation$CD, na.rm = TRUE)
    } else {
        stop('CD_corrected_for_body_surface_area_midpoint should be either numeric or character equalling median')
    }

    for (bsa_method in bsa_methods) {
        outcomeInformation <- outcomeInformation %>%
            filter(!is.na(CD)) %>% # pre-transplant samples
            mutate(CD_corrected = tac_concentration / (fin_tac_dose / .data[[bsa_method]]))
        p1 <- outcomeInformation %>%
            filter(!is.na(CD)) %>% # pre-transplant samples
            select(CD, CD_corrected, age, ageCategorical) %>%
            ggplot(aes(x = CD, y = CD_corrected, color = ageCategorical)) +
            geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
            geom_point(alpha = 0.5) +
            theme_presentation()

        p2 <- outcomeInformation %>%
            filter(!is.na(CD)) %>% # pre-transplant samples
            select(CD, CD_corrected, age, ageCategorical) %>%
            ggplot(aes(x = CD, y = CD_corrected, color = age)) +
            geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
            geom_point(alpha = 0.3) +
            theme_presentation() +
            scale_color_gradient(low = '#6ccbe1', high = 'red')

        p3 <- outcomeInformation %>%
            mutate(CD = ifelse(CD > 4, 4, CD)) %>%
            mutate(CD_corrected = ifelse(CD_corrected > 4, 4, CD_corrected)) %>%
            filter(!is.na(CD)) %>% # pre-transplant samples
            select(patientID, visitNumber, CD, CD_corrected) %>%
            pivot_longer(-c(patientID, visitNumber)) %>%
            mutate(group = str_c(patientID, visitNumber)) %>%
            ggplot(aes(x = name, y = value)) +
            geom_boxplot() +
            theme_presentation() +
            xlab("CD type") +
            ylab('value') +
            geom_point(data = data.frame(x = "CD", y = CD_midpoint), aes(x = x, y = y), color = 'red', size = 3, inherit.aes = FALSE) +
            geom_point(data = data.frame(x = "CD_corrected", y = CD_corrected_for_body_surface_area_midpoint), aes(x = x, y = y), color = 'red', size = 3, inherit.aes = FALSE)


        p4_data <- outcomeInformation %>%
            filter(!is.na(CD)) %>%
            mutate(CD = ifelse(CD > 4, 4, CD)) %>%
            mutate(CD_corrected = ifelse(CD_corrected > 4, 4, CD_corrected)) %>%
            mutate(CD_bin = ifelse(CD > CD_midpoint, 'high', 'low')) %>%
            mutate(CD_corrected_bin = ifelse(CD_corrected > CD_corrected_for_body_surface_area_midpoint, 'high', 'low')) %>%
            select(patientID, visitNumber, CD_bin, CD_corrected_bin) %>%
            # pivot_longer(-c(patientID, visitNumber)) %>%
            # group_by(name, value) %>% tally() %>%
            # rename(CD_type = name, bracket = value) %>%
            # select(CD_type, bracket) %>%
            select(CD_bin, CD_corrected_bin) %>%
            table() %>%
            as.data.frame() %>%
            rename(count = Freq)
        p4 <- ggplot(data = p4_data, aes(x = CD_bin, y = CD_corrected_bin)) +
            # geom_bar(stat = 'identity') +
            geom_tile(aes(fill = count)) +
            geom_text(aes(label = count), color = '#888888', size = 5) +
            theme_presentation() +
            scale_fill_continuous(limits = c(0, 200), low = '#4c4c4c', high ='#f8f8f8')
        # scale_fill_manual(values = c('low' = 'darkgrey', 'high' = 'lightgray')) +
        NULL

        ggsave(plot = (p1 | p2) / (p3 | p4) + plot_layout(guides = 'collect'), filename = here(str_c('plots/KLGPG_221206/CD_vs_CD_corrected_', bsa_method, '.pdf')), width = 6.5, height = 5)

        rank_cor <- outcomeInformation %>%
            filter(!is.na(CD)) %>%
            select(CD, CD_corrected) %>%
            cor(method = 'spearman') %>%
            as.numeric()

        print(rank_cor)
    }
}