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
