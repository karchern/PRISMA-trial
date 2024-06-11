#!/bin/bash -e
#cd WGS 
#> wgs_depths
#find . | grep gz | /g/scb/zeller/karcher/anaconda3/envs/parallel-fastq-dump/bin/parallel -j 32 bash ../get_depths.sh {} >> wgs_number_lines_times_4.tab
batch=$(echo $1 | cut -d "/" -f2); sampleID=$(echo $1  | cut -d "/" -f 3); echo -e "${batch}\t${sampleID}\t$(zcat ${1} | wc -l)" 

# ... run this later from R.
# a <- read_tsv('/g/scb/zeller/karcher/PRISMA/data/wgs_number_lines_times_4.tab', col_names = F) %>% mutate(batch = X1, gigaBases = ((X3/4) * 150)/1E9) %>% ggplot(aes(x = batch, y = gigaBases)) + geom_boxplot() + theme_classic() + scale_y_log10() + theme(axis.text.x = element_text(angle = 45 , hjust = 1))
