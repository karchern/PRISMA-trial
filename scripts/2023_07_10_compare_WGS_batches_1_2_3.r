library(tidyverse)

#setwd("/Volumes/zeller/karcher/PRISMA/")
setwd("/g/scb/zeller/karcher/PRISMA")

################# CAREFUL ##################
# When comparing WGS and 16S reads, make sure to not join by sampleID!
# You will find overlap but the matching is incorrect!
# Instead always join by c(PSN, Visit) to ensure correct matching!!
############################################
