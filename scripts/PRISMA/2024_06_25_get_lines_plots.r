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


##########################################
##########################################
# Disclaimer: This code has not been tested since the big reshuffling. Should you want to revive this, you need to do so carefully
##########################################
##########################################

library(pROC)

# cdModelDataBig %>%
#     select(patientID, Roseburia, Enterococcus, `Erysipelatoc.`, Coprococcus) %>%
#     left_join(clinicalMetadata)
data <- rbind(profiles %>%
    select(PSN, visit, genus, relAb) %>%
    rename(patientID = PSN) %>%
    mutate(relAb = (10^relAb) * 100),
outcomeInformation %>%
    select(patientID, visitNumber, CD) %>%
    filter(!is.na(CD)) %>%
    rename(visit = visitNumber) %>%
    mutate(genus = "CD") %>%
    rename(relAb = CD)) %>%
    arrange(patientID, visit) %>%
    # filter(genus == "CD") %>%
    inner_join(data.frame(genus = c("CD", candidateGenera, "Erysipelatoclostridium"))) %>%
    group_by(genus) %>%
    # mutate(relAb = scale(relAb)[, 1]) %>%
    identity()
# colnames(data)[colnames(data) == "relAb[,1]"] <- 'relAb'

get_line_plot <- function(PSN = NULL) {

    tmp <- c(
        "4" = "week 1",
        "5" = "week 4",
        "6" = "Month 3",
        "7" = "Month 6"
    )

    data$measurement <- data$genus

    data <- data %>%
        filter(patientID == PSN) %>%
        filter(visit > 3)

    data <- data %>%
        filter(measurement %in% c("CD", "Roseburia", "Coprococcus", "Erysipelatoclostridium"))

    scaleFactor <- max(data %>% filter(measurement != "CD") %>% pull(relAb)) / max(data %>% filter(measurement == "CD") %>% pull(relAb))

    p <- ggplot() +
        geom_line(data = data %>% filter(measurement != "CD"), aes(x = visit, y = relAb, color = measurement, group = measurement)) +
        geom_line(data = data %>% filter(measurement == "CD"), aes(x = visit, y = relAb * scaleFactor, color = measurement, group = measurement), color = 'grey', linetype = 2) +
        theme_presentation() +
        ylab("Bacterial\nrelative abundances [%]") +
        # ylab("bacterial abundance/\nCD ratio") +
        scale_y_continuous(name = "Bacterial\nrelative abundances [%]", sec.axis = sec_axis(~ . / scaleFactor, name = "CD-ratio")) +
        ggtitle(str_c(PSN)) +
        scale_x_continuous(labels = tmp, breaks = as.numeric(names(tmp))) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1))

    return(p)
}

plots <- list()
for (PSN in unique(data$patientID)) {
    p <- get_line_plot(PSN = PSN)
    plots[[length(plots) + 1]] <- p
    llAppend(plots, p, PSN)
}
names(plots) <- unique(data$patientID)

map2(plots, names(plots), \(p, n){
    ggsave(plot = p, filename = str_c(here("plots/KLGPG_221206/cd_lineplots_per_patient/", n, ".pdf")), width = 5.5, height = 3.5)
})
