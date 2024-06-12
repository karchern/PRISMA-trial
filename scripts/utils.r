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

clean_patient_clinical_metadata <- function(df_all) {
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
