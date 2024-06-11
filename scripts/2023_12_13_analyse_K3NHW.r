library(tidyverse)
library(readxl)
library(ggsignif)
library(vegan)
library(patchwork)

# Just some convenience functions that I've gotten used to
source('/home/karcher/utils/utils.r')

ssDepth <- 10000
pseudoCount <- 1/ssDepth

.f_resolve_taxonomy <- function(collated_mat, taxLevel = "genus") {
    # resolves taxonomy for collated vknight results of MAPseq and mOTUs3.(post-v0.11.3 patch)
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

meta <- read_excel('/g/scb/zeller/karcher/PRISMA/data/16S_metadata/220310_K3NHW/16S_Sequencing_Overview.xlsx', sheet = 1) %>%
    select(Sample_Type, Genecore_ID, Samples, oxygen_condition, `Drug Pool`) %>%
    mutate(Genecore_ID = str_replace(Genecore_ID, "_[0]+", "")) %>%
    mutate(Genecore_ID = str_replace(Genecore_ID, "_", "")) %>%
    rename(originalCommunity = Samples, sampleType = Sample_Type, sampleID = Genecore_ID, drugPool = `Drug Pool`) %>%
    # I'm trying to fix the weird entries for the drug incubations...
    mutate(originalCommunity = case_when(
        originalCommunity == 'MB0015 080721' ~ "MB015",
        originalCommunity == 'MB0016 290721' ~ "MB016",
        originalCommunity == 'MB0017 290721' ~ "MB017",
        originalCommunity == 'MB014 010721' ~ "MB014",
        originalCommunity == 'MB013 010721' ~ "MB013",
        .default = originalCommunity
    )) %>%
    mutate(originalCommunity = case_when(
        originalCommunity == "MB006A" & oxygen_condition == "AA" ~ "MB006",
        .default = originalCommunity  
    )) %>%
    filter(originalCommunity != "MB006S") %>%
        filter(originalCommunity != "MB010 160321") %>%
        mutate(originalCommunity = ifelse(originalCommunity == "MB010 290720", "MB010", originalCommunity)) %>%
        anti_join(data.frame(originalCommunity = c("MB012", "MB021"))) %>%
        # mutate(originalcommunity = ifelse(originalCommunity == "MB012", "MB021", originalCommunity)) %>%
        # filter(drugPool != "Pool_DMSO")
        identity() %>%
        group_by(originalCommunity) %>%
        nest() %>%
        mutate(data = map(data, \(x) {
            if (!any(str_detect(x$drugPool, "DMSO"))) {
                return(x)
            } else {
                x$drugPool[str_detect(x$drugPool,"DMSO")] <- str_c(x$drugPool[str_detect(x$drugPool,"DMSO")], 1:length(x$drugPool[str_detect(x$drugPool,"DMSO")]), sep = "_")
                return(x)
            }
        })) %>%
        unnest()
        
profiles <- readRDS('/g/scb/zeller/karcher/PRISMA/profiles/16S/220310_K3NHW/res_mapseq.rds')
colnames(profiles) <- str_replace(colnames(profiles), ".*_lane1", "")
profilesGenus <- .f_resolve_taxonomy(profiles, "genus")
profilesGenusLong <- profilesGenus %>%
    as.data.frame() %>%
    rownames_to_column('genus') %>%
    pivot_longer(-genus) %>%
    rename(sampleID = name, count = value) %>%
    inner_join(meta, by = 'sampleID')

(profilesGenusLong %>%
    group_by(sampleID, sampleType) %>%
    summarize(totalDepth = sum(count)) %>%
    ggplot(aes(x = totalDepth, fill = sampleType), alpha = 0.3) +
    geom_histogram() +
    theme_classic() +
    xlab("Sequencing depth")) %>%
    ggsave(filename = "/g/scb/zeller/karcher/PRISMA/plots/220310_K3NHW/depth_histo.pdf", width = 4, height = 4)

(profilesGenusLong %>%
    group_by(sampleID, sampleType) %>%
        summarize(totalDepth = sum(count)) %>%
        ggplot(aes(x = sampleType, y = totalDepth, fill = sampleType)) +
        geom_boxplot() +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        geom_hline(yintercept = ssDepth) +
        geom_text(data = profilesGenusLong %>%
            group_by(sampleID, sampleType) %>%
            summarize(totalDepth = sum(count)) %>%
            filter(totalDepth < ssDepth) %>%
            group_by(sampleType) %>%
            tally(), aes(x = sampleType, y = ssDepth * 0.75, label = n), inherit.aes = F, nudge_x = -0.25) +
        geom_text(data = profilesGenusLong %>%
            distinct(sampleID, sampleType) %>%
            group_by(sampleType) %>%
            tally(), aes(x = sampleType, y = 60000, label = n), inherit.aes = F) +        
    xlab("Sequencing depth")) %>%
    ggsave(filename = "/g/scb/zeller/karcher/PRISMA/plots/220310_K3NHW/depth_boxplot_with_test.pdf", width = 4, height = 4)

print('These are the low depth samples Im removing, mostly seem to be (negative?) controls')
profilesGenusLong %>%
            group_by(sampleID, sampleType, originalCommunity) %>%
            summarize(totalDepth = sum(count)) %>%
            filter(totalDepth < ssDepth)

# Remove low depth samples
profilesGenusLong <- profilesGenusLong %>%
    inner_join(profilesGenusLong %>%
        group_by(sampleID, sampleType, originalCommunity) %>%
        summarize(totalDepth = sum(count)) %>%
        filter(totalDepth >= ssDepth) %>%
        select(sampleID) %>%
        distinct())

# I'm removing MB019, MB020 from this experiment since glycerol stock sequencing failed, and hence we have nothing to compare to.MB019
profilesGenusLong <- profilesGenusLong %>%
    anti_join(data.frame(originalCommunity = c("MB019", "MB020")))

# rarefy
set.seed(1)
profilesGenusLong <- profilesGenusLong %>%
    pivot_wider(id_cols = sampleID, names_from = genus, values_from = count) %>%
    column_to_rownames("sampleID") %>%
    rrarefy(sample = ssDepth) %>%
    as.data.frame() %>%
    rownames_to_column('sampleID') %>%
    pivot_longer(-sampleID) %>%
    rename(genus = name, count = value) %>%
    inner_join(meta, by = 'sampleID')

profilesGenusLongFinal <- profilesGenusLong %>%
    group_by(sampleID, sampleType, originalCommunity) %>%
    mutate(relAb = count / sum(count)) %>%
    filter(genus != "not_resolved")

pcoa <- profilesGenusLongFinal %>%
    mutate(relAb = log10(relAb + pseudoCount)) %>%
    pivot_wider(id_cols = sampleID, names_from = genus, values_from = relAb) %>%
    column_to_rownames('sampleID') %>%
    as.matrix() %>%
    vegdist(method = "euclidean") %>%
    cmdscale(k = 2) %>%
    as.data.frame() %>%
    rownames_to_column('sampleID') %>%
    as_tibble() %>%
    rename(`PCo 1` = V1, `PCo 2` = V2) %>%
    left_join(meta, by = 'sampleID')



ed <- function(c1, c2) {
    return(sqrt((c1[1] - c2[1]) ^ 2  + (c1[2] - c2[2])  ^ 2))
}

data <- pcoa %>%
    filter(sampleType %in% c('overnight_culture', "glycerol_stock")) %>%
    rbind(pcoa %>% filter(sampleType == "glycerol_stock") %>% mutate(oxygen_condition = "MA")) %>%
    # group_by(originalCommunity) %>%
    # nest() %>%
    # mutate(data = map(data, \(x) {
    #     gsI <- which(x$sampleType == "glycerol_stock")
    #     if (length(gsI) == 0) {
    #         x$lineGroup <- NA
    #         return(x)
    #     }
    #     print(gsI)
    #     eds <- list()
    #     for (i in 1:dim(x)[1]) {
    #         if (i != gsI) {
    #             eds[[length(eds) + 1]] <- ed(c(x$`PCo 1`[i], x$`PCo 2`[i]), c(x$`PCo 1`[gsI], x$`PCo 2`[gsI]))
    #         }
    #     }
    #     eds <- unlist(eds)
    #     i <- which.min(eds)
    #     x$lineGroup <- ifelse(1:dim(x)[1] %in% c(i, gsI), originalCommunity, NA)
    #     return(x)
    # })) %>%
    # unnest()
    identity()

    (ggplot() +
    geom_point(data = data, aes(x = `PCo 1`, y = `PCo 2`, size = sampleType, color = originalCommunity)) +
        geom_line(data = data %>% filter(!is.na(originalCommunity)), aes(x = `PCo 1`, y = `PCo 2`, group = originalCommunity)) +
        scale_size_manual(values = c("glycerol_stock" = 3, 'overnight_culture' = 1)) +
    theme_classic() +
    facet_grid(oxygen_condition~.)
    ) %>%
    ggsave(filename = "/g/scb/zeller/karcher/PRISMA/plots/220310_K3NHW/pcoa_overview_overnight.pdf", width = 6.5, height = 6.5 )

(pcoa %>%
    filter(sampleType %in% c('drug_plate_10', "glycerol_stock")) %>%
        inner_join(pcoa %>%
            filter(sampleType %in% c('drug_plate_10', "glycerol_stock")) %>%
            group_by(originalCommunity) %>%
            tally() %>%
            filter(n > 10) %>%
            select(originalCommunity)) %>%
        #filter(drugPool != "Pool_DMSO") %>%
        filter(!str_detect(drugPool, "Pool_DMSO")) %>%
    ggplot() +
    geom_point(aes(x = `PCo 1`, y = `PCo 2`, size = sampleType)) +
    geom_line(aes(x = `PCo 1`, y = `PCo 2`, group = originalCommunity, color = originalCommunity)) +
    scale_size_manual(values = c("glycerol_stock" = 3, 'drug_plate_10'= 1)) +
    theme_classic()) %>%
    ggsave(filename = "/g/scb/zeller/karcher/PRISMA/plots/220310_K3NHW/pcoa_overview_drug_incubation.pdf", width = 6, height = 5 )

overnight_culture_scatters <- profilesGenusLongFinal %>%
    filter(sampleType %in% c('overnight_culture', "glycerol_stock")) %>%
    rbind(profilesGenusLongFinal %>% filter(sampleType == "glycerol_stock") %>% mutate(oxygen_condition = "MA")) %>%
    group_by(originalCommunity, oxygen_condition) %>%
    nest() %>%
    arrange(originalCommunity) %>%
    filter(originalCommunity != "Control") %>%
    mutate(data = map(data, \(x) {
        return(x %>%
            mutate(sampleType = factor(sampleType, levels = c("glycerol_stock", unique(sampleType)[unique(sampleType) != "glycerol_stock"]))) %>%
            pivot_wider(id_cols = genus, names_from = sampleType, values_from = relAb) %>%
            filter(overnight_culture != 0 | glycerol_stock != 0) %>%
            # mutate(across(all_of(c("overnight_culture", "glycerol_stock")), function(x) log10(x + pseudoCount))))
            identity())
    })) %>%
    mutate(pearsonCor = map_dbl(data, \(x) {
        # print(x)
        return(cor(x = log10(x$glycerol_stock + pseudoCount), y = log10(x$overnight_culture + pseudoCount), method = 'pearson'))
    }))


p <- overnight_culture_scatters %>%
    select(originalCommunity, oxygen_condition, data) %>%
    unnest() %>%
    ggplot(aes(x = glycerol_stock + pseudoCount, y = overnight_culture + pseudoCount)) +
    geom_point(alpha = 0.3) +
    scale_x_continuous(trans = 'log10') +
    scale_y_continuous(trans = 'log10') +
    theme_classic() +
    facet_wrap(originalCommunity ~ oxygen_condition) +
    xlab("Glycerol stock") +
    ylab("Overnight culture") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(data = overnight_culture_scatters, aes(x = 5E-4, y = 0.5, label = round(pearsonCor, 3)), inherit.aes = FALSE) +
    #annotate("text", x = -3, y = -1, label = round(pc, 3)) +
    #ggtitle(str_c(o, ox, sep = ","))
    NULL

#ggsave(plot = wrap_plots(overnight_culture_scatters$plots, nrow = 5), filename = "/g/scb/zeller/karcher/PRISMA/plots/220310_K3NHW/overnight_culture_scatters.pdf", width = 10, height = 8.25)
ggsave(plot = p, filename = "/g/scb/zeller/karcher/PRISMA/plots/220310_K3NHW/overnight_culture_scatters.pdf", width = 9.25, height = 10.5)

drug_incubation_scatters <- profilesGenusLongFinal %>%
    filter(sampleType %in% c('drug_plate_10', "glycerol_stock")) %>%
    inner_join(profilesGenusLongFinal %>%
        filter(sampleType %in% c('drug_plate_10', "glycerol_stock")) %>%
        group_by(originalCommunity) %>%
        distinct(sampleID) %>%
        tally() %>%
        filter(n > 10) %>%
        select(originalCommunity)) %>%
        #filter(drugPool != "Pool_DMSO") %>%
        filter(!str_detect(drugPool, "Pool_DMSO")) %>%
    mutate(sampleType = ifelse(sampleType == "glycerol_stock", sampleType, drugPool)) %>%
    group_by(originalCommunity, oxygen_condition) %>%
    nest() %>%
    arrange(originalCommunity) %>%
    filter(originalCommunity != "Control") %>%
    mutate(data = map(data, \(x) {
        return(x %>%
            mutate(sampleType = factor(sampleType, levels = c("glycerol_stock", unique(sampleType)[unique(sampleType) != "glycerol_stock"]))) %>%
            pivot_wider(id_cols = genus, names_from = sampleType, values_from = relAb) %>%
            pivot_longer(-c(genus, glycerol_stock)) %>%
            rename(pool = name, relAb = value) %>%
            filter(glycerol_stock != 0 | relAb != 0) %>%
            # mutate(
            #     glycerol_stock = log10(glycerol_stock + pseudoCount),
            #     relAb = log10(relAb + pseudoCount),
            # ))
            identity())
    })) %>%
    mutate(pearsonCor = map_dbl(data, \(x) {
        return(cor(x = log10(x$glycerol_stock + pseudoCount), y = log10(x$relAb + pseudoCount), method = 'pearson'))
    }))

p <- drug_incubation_scatters %>%
    select(originalCommunity, oxygen_condition, data) %>%
    unnest() %>%
    ggplot(aes(x = glycerol_stock + pseudoCount, y = relAb + pseudoCount, color = pool)) +
    geom_point(alpha = 0.2) +
    scale_x_continuous(trans = 'log10') +
    scale_y_continuous(trans = 'log10') +
    theme_classic() +
    # annotate("text", x = -3, y = -1, label = round(pc, 3)) +
    geom_text(data = drug_incubation_scatters, aes(x = 5E-4, y = 1, label = round(pearsonCor, 3)), inherit.aes = FALSE) +
    facet_wrap(originalCommunity ~ .) +
    xlab("Glycerol stock") +
    ylab("Drug incubation") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +    
    NULL

#ggsave(plot = wrap_plots(drug_incubation_scatters$plots, nrow = 2) + plot_layout(guides = 'collect'), filename = "/g/scb/zeller/karcher/PRISMA/plots/220310_K3NHW/drug_incubation_scatters.pdf", width = 5.75, height = 4)
ggsave(plot = p, filename = "/g/scb/zeller/karcher/PRISMA/plots/220310_K3NHW/drug_incubation_scatters.pdf", width = 5.75, height = 4)

dmso_vs_glycerol <- profilesGenusLongFinal %>%
    filter(sampleType %in% c('drug_plate_10', "glycerol_stock")) %>%
    inner_join(profilesGenusLongFinal %>%
        filter(sampleType %in% c('drug_plate_10', "glycerol_stock")) %>%
        group_by(originalCommunity) %>%
        distinct(sampleID) %>%
        tally() %>%
        filter(n > 10) %>%
        select(originalCommunity)) %>%
    filter(str_detect(drugPool, "Pool_DMSO") | sampleType == "glycerol_stock") %>%
    mutate(sampleType = ifelse(sampleType == "glycerol_stock", sampleType, drugPool)) %>%
    group_by(originalCommunity, oxygen_condition) %>%
    nest() %>%
    arrange(originalCommunity) %>%
    filter(originalCommunity != "Control") %>%
    mutate(data = map(data, \(x) {
        return(x %>%
            mutate(sampleType = factor(sampleType, levels = c("glycerol_stock", unique(sampleType)[unique(sampleType) != "glycerol_stock"]))) %>%
            pivot_wider(id_cols = genus, names_from = sampleType, values_from = relAb) %>%
            pivot_longer(-c(genus, glycerol_stock)) %>%
            rename(pool = name, relAb = value) %>%
            filter(glycerol_stock != 0 | relAb != 0) %>%
            # mutate(
            #     glycerol_stock = log10(glycerol_stock + pseudoCount),
            #     relAb = log10(relAb + pseudoCount),
            # ))
            identity())
    })) %>%
    mutate(pearsonCor = map_dbl(data, \(x) {
        return(cor(x = log10(x$glycerol_stock + pseudoCount), y = log10(x$relAb + pseudoCount), method = 'pearson'))
    }))

p <- dmso_vs_glycerol %>%
    select(originalCommunity, oxygen_condition, data) %>%
    unnest() %>%
    ggplot(aes(x = glycerol_stock + pseudoCount, y = relAb + pseudoCount, color = pool)) +
    geom_point(alpha = 0.2) +
    scale_x_continuous(trans = 'log10') +
    scale_y_continuous(trans = 'log10') +
    theme_classic() +
    # annotate("text", x = -3, y = -1, label = round(pc, 3)) +
    geom_text(data = dmso_vs_glycerol, aes(x = 5E-4, y = 1, label = round(pearsonCor, 3)), inherit.aes = FALSE) +
    facet_wrap(originalCommunity ~ .) +
    xlab("Glycerol stock") +
    ylab("DMSO-incubated") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +    
    NULL

#ggsave(plot = wrap_plots(drug_incubation_scatters$plots, nrow = 2) + plot_layout(guides = 'collect'), filename = "/g/scb/zeller/karcher/PRISMA/plots/220310_K3NHW/drug_incubation_scatters.pdf", width = 5.75, height = 4)
ggsave(plot = p, filename = "/g/scb/zeller/karcher/PRISMA/plots/220310_K3NHW/dmso_vs_glycerol_stock_scatters.pdf", width = 5.75, height = 4)


dmso_vs_drugs <- profilesGenusLongFinal %>%
    filter(sampleType %in% c('drug_plate_10')) %>%
    inner_join(profilesGenusLongFinal %>%
        filter(sampleType %in% c('drug_plate_10')) %>%
        group_by(originalCommunity) %>%
        distinct(sampleID) %>%
        tally() %>%
        filter(n > 10) %>%
        select(originalCommunity)) %>%
    filter(str_detect(drugPool, "Pool")) %>%
    mutate(sampleType = ifelse(sampleType == "glycerol_stock", sampleType, drugPool)) %>%
    group_by(originalCommunity, oxygen_condition) %>%
    nest() %>%
    arrange(originalCommunity) %>%
    filter(originalCommunity != "Control") %>%
    mutate(data = map(data, \(x) {
        return(x %>%
            mutate(sampleType = factor(sampleType, levels = c("glycerol_stock", unique(sampleType)[unique(sampleType) != "glycerol_stock"]))) %>%
                pivot_wider(id_cols = genus, names_from = sampleType, values_from = relAb) %>%
                pivot_longer(-c(genus, Pool_DMSO_1, Pool_DMSO_2, Pool_DMSO_3, Pool_DMSO_4)) %>%
                rename(pool = name, relAb = value) %>%
                # select(-pool) %>%
                pivot_longer(-c(genus, pool, relAb)) %>%
            filter(value != 0 | relAb != 0) %>%
            # mutate(
            #     glycerol_stock = log10(glycerol_stock + pseudoCount),
            #     relAb = log10(relAb + pseudoCount),
            # ))
            identity())
    })) %>%
    mutate(pearsonCor = map_dbl(data, \(x) {
        return(cor(x = log10(x$value + pseudoCount), y = log10(x$relAb + pseudoCount), method = 'pearson'))
    }))

p <- dmso_vs_drugs %>%
    select(originalCommunity, oxygen_condition, data) %>%
    unnest() %>%
    ggplot(aes(x = value + pseudoCount, y = relAb + pseudoCount)) +
    geom_point(alpha = 0.2) +
    scale_x_continuous(trans = 'log10') +
    scale_y_continuous(trans = 'log10') +
    theme_classic() +
    # annotate("text", x = -3, y = -1, label = round(pc, 3)) +
    geom_text(data = dmso_vs_drugs, aes(x = 5E-4, y = 1, label = round(pearsonCor, 3)), inherit.aes = FALSE) +
    facet_wrap(originalCommunity ~ .) +
    xlab("Glycerol stock") +
    ylab("DMSO-incubated") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +    
    NULL

#ggsave(plot = wrap_plots(drug_incubation_scatters$plots, nrow = 2) + plot_layout(guides = 'collect'), filename = "/g/scb/zeller/karcher/PRISMA/plots/220310_K3NHW/drug_incubation_scatters.pdf", width = 5.75, height = 4)
ggsave(plot = p, filename = "/g/scb/zeller/karcher/PRISMA/plots/220310_K3NHW/dmso_vs_drugs_stock_scatters.pdf", width = 5.75, height = 4)


(profilesGenusLongFinal %>%
    filter(sampleType %in% c('overnight_culture', "glycerol_stock")) %>%
    rbind(profilesGenusLongFinal %>% filter(sampleType == "glycerol_stock") %>% mutate(oxygen_condition = "MA")) %>%
    group_by(originalCommunity, oxygen_condition) %>%
    nest() %>%
    arrange(originalCommunity) %>%
    filter(originalCommunity != "Control") %>%
    mutate(data = map(data, \(x) {
        return(x %>%
            mutate(sampleType = factor(sampleType, levels = c("glycerol_stock", unique(sampleType)[unique(sampleType) != "glycerol_stock"]))) %>%
            pivot_wider(id_cols = genus, names_from = sampleType, values_from = relAb) %>%
            filter(overnight_culture != 0 | glycerol_stock != 0) %>%
            mutate(across(all_of(c("overnight_culture", "glycerol_stock")), function(x) x > pseudoCount))) %>%
            filter(overnight_culture)
    })) %>%
    mutate(fractionGlycStockGeneraRetained = map_dbl(data, \(x) {
        return(mean(x$glycerol_stock))
    })) %>%
        ggplot(aes(x = originalCommunity, y = fractionGlycStockGeneraRetained, fill = oxygen_condition)) +
            geom_bar(stat = 'identity', position = 'dodge') +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            ylab("Frac. of glycerol stock genera\nretained in overnight cultures"))  %>%
    ggsave(filename = "/g/scb/zeller/karcher/PRISMA/plots/220310_K3NHW/fraction_retained_genera_overnight_cultures.pdf", width = 5, height = 3)

(profilesGenusLongFinal %>%
    filter(sampleType %in% c('drug_plate_10', "glycerol_stock")) %>%
    inner_join(profilesGenusLongFinal %>%
    filter(sampleType %in% c('drug_plate_10', "glycerol_stock")) %>%
        group_by(originalCommunity) %>%
        distinct(sampleID) %>%
        tally() %>%
        filter(n>10) %>%
    select(originalCommunity)) %>%   
    mutate(sampleType = ifelse(sampleType == "glycerol_stock", sampleType, drugPool)) %>%
    group_by(originalCommunity, oxygen_condition) %>%
    nest() %>%
    arrange(originalCommunity) %>%
    filter(originalCommunity != "Control") %>%
    mutate(data = map(data, \(x) {
        return(x %>%
            mutate(sampleType = factor(sampleType, levels = c("glycerol_stock", unique(sampleType)[unique(sampleType) != "glycerol_stock"]))) %>%
            pivot_wider(id_cols = genus, names_from = sampleType, values_from = relAb) %>%
            pivot_longer(-c(genus, glycerol_stock)) %>%
            rename(pool = name, relAb = value) %>%
            filter(glycerol_stock != 0 | relAb != 0) %>%
            mutate(
                glycerol_stock = glycerol_stock > pseudoCount,
                relAb = relAb + pseudoCount,
            )) %>%
            filter(glycerol_stock)
    })) %>%
    mutate(fractionGlycStockGeneraRetained = map_dbl(data, \(x) {
        return(mean(x$glycerol_stock))
    })) %>% ggplot(aes(x = originalCommunity, y = fractionGlycStockGeneraRetained, fill = oxygen_condition)) +
            geom_bar(stat = 'identity', position = 'dodge') +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            ylab("Frac. of glycerol stock genera\nretained after drug incubation"))  %>%
    ggsave(filename = "/g/scb/zeller/karcher/PRISMA/plots/220310_K3NHW/fraction_retained_genera_drug_incubation.pdf", width = 3.5, height = 3)
