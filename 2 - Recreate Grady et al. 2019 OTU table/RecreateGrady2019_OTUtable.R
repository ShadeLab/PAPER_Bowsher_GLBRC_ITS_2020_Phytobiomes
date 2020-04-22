
###
#Download directory "R" from Grady et al. 2019 Github repository:
##https://github.com/ShadeLab/PAPER_GradySorensenStopnisek_NatComm_2019
#Portions of the below code were copied directly from "Workflow_GLBRC16S.R" to replicate 
#the OTU table used by Grady et al. 2019
###


library(dplyr)
library(tidyr)
library(reshape2)
library(RSQLite)
library(stringr)
# Read in OTU table
otu <- read.table("R/InputFiles/table_combined_merged_trimmed_otus.txt",sep="\t", header=TRUE, stringsAsFactors = FALSE, row.names=1)

#Prepping the environmental metadata
glbrc <- dbConnect(SQLite(), dbname="R/InputFiles/GLBRC_bioenergy_db.db" )

#get tables into R
glbrc_NA <- dbGetQuery(glbrc, "select * from nucleic_acids") 
glbrc_soil <- dbGetQuery(glbrc, 'select * from soil')
glbrc_plant <- dbGetQuery(glbrc, 'select * from plant')
glbrc_plot <- dbGetQuery(glbrc, 'select * from plot')
glbrc_sampling <- dbGetQuery(glbrc, 'select * from sampling')
glbrc_sequncing <- dbGetQuery(glbrc, 'select * from sequencing')
glbrc_sequncing <- glbrc_sequncing %>% 
  mutate(nucleic_acid_name = str_trim(glbrc_sequncing$nucleic_acid_name, side = "both"))

#joining tables to create complete map file
metadata <- full_join(glbrc_plot, glbrc_sampling, by='plotID')
metadata <- full_join(metadata, glbrc_soil, by='sampleID')
metadata <- full_join(metadata, glbrc_plant, by='sampleID')
metadata <- full_join(metadata, glbrc_NA, by='sampleID')
metadata <- full_join(metadata, glbrc_sequncing, by='nucleic_acid_name')

map <- metadata

map_16 <- subset(map, map$primers == 'EMP V4')
map_16$help_name = as.character(lapply(strsplit(as.character(map_16$nucleic_acid_name), split="D"), "[", 1))

#finding duplicates
n_occur <- data.frame(table(map_16$help_name))
n_occur[n_occur$Freq > 1,]
duplicate_df <- map_16[map_16$help_name %in% n_occur$Var1[n_occur$Freq > 1],]
list_dupli <- duplicate_df$sequence_name #list of duplicate samples

D1 <- duplicate_df[grep('D1', duplicate_df$sequence_name),]
D1$removing <- 'remove'
dim(D1)

map_full <- full_join(map, D1) 

#creating numeric time column
map_full$sampling_date <- paste0(map_full$month,'-', map_full$day,'-',map_full$year) 
map_full$sampling_date <- as.POSIXct(map_full$sampling_date, format='%m-%d-%Y')
time <- as.POSIXct(map_full$sampling_date, format='%m-%d-%Y')
time_numeric <- as.numeric(time)
map_time <- cbind(map_full, time_numeric)

#adding weather data
weather <- read.csv('R/InputFiles/kbs_weather_09212017.csv', encoding = 'UTF-8', na.strings= c("NA", " ", ""))
dim(weather)
head(weather)
weather$sampling_date <- as.POSIXct(weather$date, format='%d.%m.%y')

sub_weather <- weather[weather$sampling_date %in% map_time$sampling_date,] #subsetting weather file for sample dates

#merging dataframes - map file and weather
map_complete <- full_join(map_time, weather)

#Subset to 16S Samples
map_16S <- subset(map_complete, map_complete$primers=="EMP V4")
# Get the name of the samples we have data for
# Remove any duplicates that we don't want to use
map_small <- map_16S[map_16S$exclude_from_analysis=="N",]

#map_small_v1 <- subset(map_16S_small, map_16S_small$removing == 'remove') 
#map_small <- map_16S_small[!(rownames(map_16S_small) %in% c(33,108, 185,259,37,418,467,487,571)),]

samples <- colnames(otu)
# put taxonomy into its own variable
taxonomy <- read.csv('R/InputFiles/taxonomy_combined_merged_trimmed_otus.csv', header = T, row.names = 1, na.strings= c("NA", " ", ""))
# subset the map to include only those samples we have sequence data
map_small <- map_small[map_small$sequence_name %in% samples,]


# Remove rows that have N/A for the sequence name
map_small<- map_small[complete.cases(map_small$sequence_name),]

# Subset the samples to those we want to analyse (IE remove the duplicates)
samples <- samples[samples %in% map_small$sequence_name]
#Subset the OTU table to only the samples we want to analyze(IE remove the duplicates)
otu_sub <- otu[,colnames(otu) %in% samples]
# Order the samples
otu_sub <- otu_sub[,order(colnames(otu_sub))]
# Order the samples of the map the same way
map_small <- map_small[order(map_small$sequence_name),]
# Check to make sure they all match with each other
colnames(otu_sub) == map_small$sequence_name

# Map file that has duplicates of samples removed
map_16S <- map_small
# OTU table that removes duplicates of single samples
otu <- otu_sub
otu_CM <- otu
#removing Eukaryota from the OTU table
tax_short <- taxonomy[!grepl("Mitochondria", taxonomy$Family),]
tax_short <- tax_short[!grepl("Chloroplast", tax_short$Class),]

otu <- otu[rownames(otu) %in% rownames(tax_short),]

taxonomy_full <- taxonomy
taxonomy <- taxonomy[rowSums(otu_sub)>0,]
otu <- otu[rowSums(otu)>0,]
otu_soil <- otu[,map_16S$source=="soil"]

tax_filtered <- tax_short%>%
  mutate(otu = rownames(tax_short)) %>%
  filter(!is.na(Phylum))

silva_bact_only <- read.csv('R/InputFiles/silva_bacteria_only_glbrc.csv', header=T)

keep_otus <- silva_bact_only %>%
  filter(lca_tax_slv != 'Unclassified;') %>%
  separate(lca_tax_slv, into=c("Kingdom", "Phylum", "Class", 
                               "Order", "Family", "Genus", "Species"), sep=";", remove=F) %>%
  mutate(taxonomy = lca_tax_slv) %>%
  filter(Kingdom!= 'Eukaryota') %>%
  select(-lca_tax_slv) %>%
  bind_rows(tax_filtered) %>%
  select(-taxonomy)

otu_filtered <- otu[rownames(otu) %in% keep_otus$otu,]

library(vegan)
set.seed(13)
#otu_rare <- t(rrarefy(t(otu), min(colSums(otu))))
otu_rare<- t(rrarefy(t(otu_filtered), 1000))












#next, get the 2017 switchgrass phyllosphere coreOTUs.


#Core analysis V2
map_16S$sampling_week <- 0
map_16S$sampling_week[map_16S$sampling_date == '2017-05-15 EDT'] <- 1 
map_16S$sampling_week[map_16S$sampling_date == '2017-06-05 EDT'] <- 2 
map_16S$sampling_week[map_16S$sampling_date == '2017-06-26 EDT'] <- 3
map_16S$sampling_week[map_16S$sampling_date == '2017-07-17 EDT'] <- 4
map_16S$sampling_week[map_16S$sampling_date == '2017-08-07 EDT'] <- 5 
map_16S$sampling_week[map_16S$sampling_date == '2017-08-28 EDT'] <- 6 
map_16S$sampling_week[map_16S$sampling_date == '2017-09-18 EDT'] <- 7 

swit17_otu <- otu_rare[,map_16S$plant=="switchgrass" & (map_16S$source=='phyllosphere') & (map_16S$year==2017)]

rel_otu_rare <- decostand(otu_rare, method="total", MARGIN=2)
selected_otus <- data.frame(otu = as.factor(row.names(rel_otu_rare)), rel_otu_rare) %>% gather(sequence_name, abun, -otu) %>% 
  left_join(map_16S[, c('sequence_name','rep','time_numeric', 'treatment' ,'source', 'plant', 'month', 'sampling_date', 
                        'carbon_per_nitrogen', 'nitrogen_percent', 'carbon_percent', 'year')], by = 'sequence_name') %>%
  left_join(tax_filtered, by='otu')


source='switchgrass 2017'
map_file <- map_16S %>%
  filter(plant == 'switchgrass',
         source == 'phyllosphere',
         year == '2017')
input <- swit17_otu
weeks <- c(1:10)
result <- NULL

for(i in weeks) {
  if(i %in% unique(map_file$sampling_week)) {
    name_sample <- map_file$sequence_name[map_file$sampling_week == i]
    timed_matrix <- input[,colnames(input) %in% name_sample]
    timed_matrix <- timed_matrix[rowSums(timed_matrix)>0,]
    timed_matrix_PA <- 1*((timed_matrix>0)==1)
    timed_matrix_PA <- timed_matrix_PA[rowSums(timed_matrix_PA)>0,]
    Occ <- rowSums(timed_matrix_PA)/ncol(timed_matrix_PA)
    rel_abun <- decostand(timed_matrix_PA, method="total", MARGIN=2)
    Mean_rel_abund <- apply(rel_abun, 1, mean)
    df_o <- data.frame(otu=names(Occ), occ=Occ) 
    df_a <- data.frame(otu=names(Mean_rel_abund), abun=log10(Mean_rel_abund))
    table <- left_join(df_a, df_o, by='otu') %>% mutate(week = i, source = source)
    result <- rbind(result, table)
  }
  else {
  }
}
occ1_s17 <- result[result$occ==1,]
tempCore_sw17 <- as.character(unique(occ1_s17$otu))

colSums(otu_rare)
#save rarefied OTU table
write.table(otu_rare, "otu_rare_16s.txt", quote=FALSE, sep='\t')
#save core OTU list
write.table(tempCore_sw17, "core_16s.txt", quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
#save taxonomy
write.table(tax_filtered, "tax_16s.txt", quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
