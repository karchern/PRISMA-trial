# Define some convenience functions

time_point_colors <- c(
    "#a5c2fc", # light blue
    "#82aafa", # slightly less intense blue
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

    p <- ggplot(tmp2, aes(x = as.factor(anyABx), y = relAb)) +
        # geom_boxplot(aes_string(x = "anyABx", y = "relAb", fill = complication)) +
        geom_quantileplot(aes(fill = CDbinary), quantilesP = quantilesP) +
        scale_fill_quantile(cdMetabColors, quantilesP) +
        theme_presentation()
    # # Remove NAs
    # tmp2 <- tmp2 %>%
    #     filter(!if_any(all_of(c('anyABx', complication, "relAb")), is.na))
    # p <- get_quantile_plot(tmp2, "anyABx", complication, "relAb", plotGroup = NA, expectedNumLevels = 2)
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
    p <- ggplot(data = tmp2, aes(x = relAb, y = outcome, color = PSN)) +
        geom_point() +
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

illustrate_taxon_hit <- function(modelData = NULL, taxon = NULL) {
    modelData <- modelData %>%
        filter(genus == taxon) %>%
        rename(`Tacrolimus\nmetabolism` = cdMetabolism) %>%
        mutate(relAb = (10^relAb) * 100)
    p <- ggplot() +
        geom_boxplot(
            data = modelData,
            aes(x = `Tacrolimus\nmetabolism`, y = relAb, fill = `Tacrolimus\nmetabolism`), outlier.color = NA) +
        geom_jitter(
            data = modelData,
            aes(x = `Tacrolimus\nmetabolism`, y = relAb, fill = `Tacrolimus\nmetabolism`), position = position_jitter()) +
        theme_presentation() +
        ylab("Bacterial\nrelative abundance [%]") +
        # scale_fill_manual(values = c('low' = "#4a5dca", "high" = "#d43e3e")) +
        scale_fill_manual(values = cdMetabColors) +
        scale_y_continuous(trans = 'log10', limits = c(0.005, 5))
    return(p)
}
