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
