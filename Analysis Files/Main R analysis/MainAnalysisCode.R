
library(dplyr)##dplyr_0.8.3
library(tidyr)##tidyr_0.8.3
library(RSQLite)##RSQLite_2.1.2
library(stringr)##stringr_1.4.0
library(ggplot2)##ggplot2_3.2.0
library(gridExtra)##gridExtra_2.3
library(vegan)##vegan_2.5-5
library(scales)##scales_1.0.0
library(cowplot)##cowplot_1.0.0
library(purrr) ##purrr_0.3.2
library(car) ##car_3.0-3
library(tibble) ##tibble_2.1.3 
sessionInfo()


# Read in OTU table
otu <- read.table("otu_table_ITS_UPARSE_filtered_180.txt",sep="\t",header=TRUE, row.names=1)


#Prepping the environmental metadata
glbrc <- dbConnect(SQLite(), dbname="R/InputFiles/GLBRC_bioenergy_db.db" )

#Content of the DB
dbListTables(glbrc)

#Content of selected tables in the DB
dbListFields(glbrc, 'plant')
dbListFields(glbrc, 'sequencing')
dbListFields(glbrc, 'soil')
dbListFields(glbrc, 'nucleic_acids')

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

unique(metadata$nucleic_acid_name)
map <- metadata
dim(map)

#subset metadata to only include ITS2 data.
map_ITSprimers <- subset(map, map$primers == 'ITS2')
dim(map_ITSprimers)
head(map_ITSprimers)
map_ITSprimers$help_name = as.character(lapply(strsplit(as.character(map_ITSprimers$nucleic_acid_name), split="D"), "[", 1))
unique(map_ITSprimers$help_name)

#creating numeric time column
map_ITSprimers$sampling_date <- paste0(map_ITSprimers$month,'-', map_ITSprimers$day,'-',map_ITSprimers$year) 
map_ITSprimers$sampling_date <- as.POSIXct(map_ITSprimers$sampling_date, format='%m-%d-%Y')
time <- as.POSIXct(map_ITSprimers$sampling_date, format='%m-%d-%Y')
time_numeric <- as.numeric(time)
map_time <- cbind(map_ITSprimers, time_numeric)

#adding weather data
weather <- read.csv('R/InputFiles/kbs_weather_09212017.csv', encoding = 'UTF-8', na.strings= c("NA", " ", ""))
dim(weather)
head(weather)
weather$sampling_date <- as.POSIXct(weather$date, format='%d.%m.%y')

sub_weather <- weather[weather$sampling_date %in% map_time$sampling_date,] #subsetting weather file for sample dates

#merging dataframes - map file and weather
map_complete <- full_join(map_time, weather)

#Subset to ITS Samples
map_ITS <- subset(map_complete, map_complete$primers=="ITS2")


#edit column names of otu table to match correspnding sample in the GLBRC metadata
samples <- colnames(otu)
samples
samples = gsub(pattern = "G5R1", replacement = "G5R1_", x = samples)
samples = gsub(pattern = "G5R2", replacement = "G5R2_", x = samples)
samples = gsub(pattern = "G5R3", replacement = "G5R3_", x = samples)
samples = gsub(pattern = "G5R4", replacement = "G5R4_", x = samples)
samples = gsub(pattern = "NF", replacement = "NF_", x = samples)
samples = gsub(pattern = "MAIN", replacement = "MAIN_", x = samples)
samples = gsub(pattern = "2017", replacement = "2017_", x = samples)
samples = gsub(pattern = "LD1", replacement = "LD1_", x = samples)
samples = gsub(pattern = "SD1", replacement = "SD1_", x = samples)
samples
colnames(otu)<- samples

#remove samples from the GLBRC metadata that are not in the OTU table
map_ITS <- subset(map_ITS, sequence_name %in% colnames(otu))

# Order the samples in the otu table
otu <- otu[,order(colnames(otu))]
# Order the samples of the metadata the same way
map_ITS <- map_ITS[order(map_ITS$sequence_name),]
# Check to make sure they all match with each other
colnames(otu) == map_ITS$sequence_name

#read in taxonomy 
taxonomy <- read.table('constax_taxonomy/consensus_taxonomy.txt', sep='\t',header = T, row.names = 1)

#Check for non-target taxa
unique(taxonomy$Kingdom)
unique(taxonomy$Phylum)

##see what OTUs were in the two non-fungal phyla detected above.
OTUs_nonfungal <- rownames(subset(taxonomy, Phylum=='Cercozoa' | Phylum=="Chlorophyta"))

#removing non-fungi from the OTU table
otu = otu[!row.names(otu) %in% OTUs_nonfungal, ]
#removing non-fungi from the taxonomy table
taxonomy = taxonomy[!row.names(taxonomy) %in% OTUs_nonfungal,]

#check whether rowSums =1 or 0 for any OTU
row.names(otu[which(rowSums(otu) <= 1),])

#make new taxonomy dataframe, but add OTUs as new column
keep_otus <- tibble::rownames_to_column(taxonomy, "otu")



#rarefy dataset.
set.seed(13)
otu_rare <- t(rrarefy(t(otu), min(colSums(otu))))
colSums(otu_rare)



# calculate number of OTUs per sample
rare_OTU_pa <- as.data.frame(1*((otu_rare>0)==1))
min(colSums(rare_OTU_pa))
max(colSums(rare_OTU_pa))
mean(colSums(rare_OTU_pa))


#Prepare Figure 1 (rarefaction curves)
curve_colors <- rep("darkgreen", ncol(otu))
curve_colors[map_ITS$source=="phyllosphere"] <- "green4"
curve_colors[map_ITS$source=="soil"] <- "chocolate4"

out<-rarecurve(t(otu), step=1000, sample=max(colSums(otu)), label=FALSE, col = curve_colors)

rare <- lapply(out, function(x){
 b <- as.data.frame(x)
 b <- data.frame(OTU = b[,1], raw.read = rownames(b))
 b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
 return(b)
})

sample_names <- rownames(t(otu))
names(rare) <- sample_names

library(purrr)
rare <- map_dfr(rare, function(x){
 z <- data.frame(x)
 return(z)
}, .id = "sample")

rare$source <- ifelse(grepl("SD", rare$sample), "soil", "phyllosphere")

#Figure 1
Fig1 <-ggplot(data = rare, aes(x=raw.read,y=OTU, color=source))+
geom_point(size=.1)+
 scale_color_manual(values=c("green4","chocolate4"))+
 ylab("Richness (No. OTUs)")+
 xlab("No. reads per sample")+
 scale_x_continuous(labels = comma)+
 geom_vline(xintercept=min(colSums(otu)))+
 guides(color = guide_legend(override.aes = list(size=3)))+
 theme_bw()

Fig1
 
#ggsave("Figure1.tiff",Fig1,width=6,height=4,units=c("in"))




#Adding sampling week property to metadata
map_ITS$sampling_week[map_ITS$sampling_date == '2017-04-24 EDT'] <-0
map_ITS$sampling_week[map_ITS$sampling_date == '2017-05-15 EDT'] <- 1 
map_ITS$sampling_week[map_ITS$sampling_date == '2017-06-05 EDT'] <- 2 
map_ITS$sampling_week[map_ITS$sampling_date == '2017-06-26 EDT'] <- 3
map_ITS$sampling_week[map_ITS$sampling_date == '2017-07-17 EDT'] <- 4
map_ITS$sampling_week[map_ITS$sampling_date == '2017-08-07 EDT'] <- 5 
map_ITS$sampling_week[map_ITS$sampling_date == '2017-08-28 EDT'] <- 6 
map_ITS$sampling_week[map_ITS$sampling_date == '2017-09-18 EDT'] <- 7 

#############################
#Occupancy abundance analysis
#############################

#' thanks to https://rdrr.io/github/jerryzhujian9/ezR/src/R/basic.R
blank2na = function(x,na.strings=c('','.','NA','na','N/A','n/a','NaN','nan')) {
  if (is.factor(x)) {
    lab = attr(x, 'label', exact = T)
    labs1 <- attr(x, 'labels', exact = T)
    labs2 <- attr(x, 'value.labels', exact = T)
    
    # trimws will convert factor to character
    x = trimws(x,'both')
    if (! is.null(lab)) lab = trimws(lab,'both')
    if (! is.null(labs1)) labs1 = trimws(labs1,'both')
    if (! is.null(labs2)) labs2 = trimws(labs2,'both')
    
    if (!is.null(na.strings)) {
      # convert to NA
      x[x %in% na.strings] = NA
      # also remember to remove na.strings from value labels 
      labs1 = labs1[! labs1 %in% na.strings]
      labs2 = labs2[! labs2 %in% na.strings]
    }
    
    # the levels will be reset here
    x = factor(x)
    
    if (! is.null(lab)) attr(x, 'label') <- lab
    if (! is.null(labs1)) attr(x, 'labels') <- labs1
    if (! is.null(labs2)) attr(x, 'value.labels') <- labs2
  } else if (is.character(x)) {
    lab = attr(x, 'label', exact = T)
    labs1 <- attr(x, 'labels', exact = T)
    labs2 <- attr(x, 'value.labels', exact = T)
    
    # trimws will convert factor to character
    x = trimws(x,'both')
    if (! is.null(lab)) lab = trimws(lab,'both')
    if (! is.null(labs1)) labs1 = trimws(labs1,'both')
    if (! is.null(labs2)) labs2 = trimws(labs2,'both')
    
    if (!is.null(na.strings)) {
      # convert to NA
      x[x %in% na.strings] = NA
      # also remember to remove na.strings from value labels 
      labs1 = labs1[! labs1 %in% na.strings]
      labs2 = labs2[! labs2 %in% na.strings]
    }
    
    if (! is.null(lab)) attr(x, 'label') <- lab
    if (! is.null(labs1)) attr(x, 'labels') <- labs1
    if (! is.null(labs2)) attr(x, 'value.labels') <- labs2
  } else {
    x = x
  }
  return(x)
}

lastValue <- function(x) tail(x[!is.na(x)], 1)

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)}


#in taxonomy table, add column naming the highest resolution taxonomy achieved for each OTU
rownames(keep_otus) <- keep_otus$otu
keep_otus[] = lapply(keep_otus, blank2na, na.strings=c('','NA','na','N/A','n/a','NaN','nan'))
last_taxons<- apply(keep_otus, 1, lastValue)
keep_otus$last_taxon <- last_taxons
keep_otus$final_names <- paste(keep_otus$last_taxon, keep_otus$otu, sep=' - ')



#Core analysis 
#First, convert rarefied OTU table abundances to relative abundances
rel_otu_rare <- decostand(otu_rare, method="total", MARGIN=2)

#Using relabund table (above), give each OTU in each sample a row.
#Add relative abundance, metadata, and taxonomy
selected_otus <- data.frame(otu = as.factor(row.names(rel_otu_rare)), rel_otu_rare) %>% gather(sequence_name, abun, -otu) %>% 
  left_join(map_ITS[, c('sequence_name','rep','time_numeric', 'treatment' ,'source', 'plant', 'month', 'sampling_date', 
                        'carbon_per_nitrogen', 'nitrogen_percent', 'carbon_percent', 'year')], by = 'sequence_name') %>%
  left_join(keep_otus, by='otu')

#Subset rarefied OTU table to only leaf samples
leaf_otu <- otu_rare[,map_ITS$source=="phyllosphere"]

#subset mapping data to only leaf samples
map_leaf <- map_ITS %>%
  filter(source == 'phyllosphere')


#generate abundance/occupancy data for each otu for each week
weeks <- c(0:7)
result <- NULL
source='switchgrass 2017'

for(i in weeks) {
  if(i %in% unique(map_leaf$sampling_week)) {
    name_sample <- map_leaf$sequence_name[map_leaf$sampling_week == i]
    timed_matrix <- leaf_otu[,colnames(leaf_otu) %in% name_sample]
    timed_matrix <- timed_matrix[rowSums(timed_matrix)>0,]
    timed_matrix_PA <- 1*((timed_matrix>0)==1)
    timed_matrix_PA <- timed_matrix_PA[rowSums(timed_matrix_PA)>0,]
    Occ <- rowSums(timed_matrix_PA)/ncol(timed_matrix_PA)
    rel_abun <- decostand(timed_matrix, method="total", MARGIN=2)
    Mean_rel_abund <- apply(rel_abun, 1, mean)
    df_o <- data.frame(otu=names(Occ), occ=Occ) 
    df_a <- data.frame(otu=names(Mean_rel_abund), abun=log10(Mean_rel_abund))
    table <- left_join(df_a, df_o, by='otu') %>% mutate(week = i, source = source)
    result <- rbind(result, table)
  }
  else {
  }
}

#determine otus at occupancy of 1 for any week)
occ1_s17 <- result[result$occ==1,]
fungalcore_sw17 <- as.character(unique(occ1_s17$otu))




#The above was just to determine which OTUs were at occupancy of 1 for a given week.
#now prepare the data for the actual occupancy/abundance plot.

#convert rarefied leaf OTU table to PA matrix 
leaf_otu <- leaf_otu[rowSums(leaf_otu)>0,]
leaf_otu_PA <- 1*((leaf_otu>0)==1)
Occ_leaf <- rowSums(leaf_otu_PA)/ncol(leaf_otu_PA)
leaf_otu_rel <- decostand(leaf_otu, method="total", MARGIN=2)
Mean_abund_leaf <- apply(leaf_otu_rel, 1, mean)

## get mean relative abundance of alternaria alternata, a pathogen
listing <-as.data.frame(leaf_otu['OTU_2',])
colnames(listing)<-c('rarefiedreads')
mean(listing$rarefiedreads/86710)
sd(listing$rarefiedreads/86710)

#create dataframe for occupancy data generated above.
leaf_df_occ <- data.frame(otu=names(Occ_leaf), occ=Occ_leaf) 
#create dataframe for abundance data generated above.
leaf_df_abun <- data.frame(otu=names(Mean_abund_leaf), abun=log10(Mean_abund_leaf))

#combine the two
leaf_occ_abun <- left_join(leaf_df_abun, leaf_df_occ, by='otu')

#label all otus as 'NotCore'
leaf_occ_abun$unique <- 'NotCore'
#change label  to 'core' for otus that were occ=1 at any time point.
leaf_occ_abun$unique[(leaf_occ_abun$otu %in% fungalcore_sw17)] <- 'Core (occ=1 at any time pt)'

#Figure 4A (occupancy abundance plot)
leaf_occ_abun_plot <- ggplot(data=leaf_occ_abun, aes(x=abun, y=occ, fill=unique)) +
  theme_bw()+
  geom_point(size=3, pch=21, alpha=.8) +
  scale_fill_manual(breaks=unique, values=c('green4','white')) +
  labs(x=paste('log(mean relative abundance per OTU)\n (n=',nrow(leaf_occ_abun),' OTUs)',sep=''), y=paste('Occupancy (n=',ncol(leaf_otu_PA),')', sep=''), fill=NULL) +
  theme(legend.position = 'none',
        legend.background = element_rect(fill=alpha(0.1)),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(breaks=seq(0,1,.2)) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))+
  theme(axis.title.y = element_text(size = 16))+
  theme(axis.title.x= element_text(size = 16))+
  labs(tag = "A")
leaf_occ_abun_plot


#############################
#Hierarchical clustering analysis
#############################


##All OTUs: relative abundance, metadata, and taxonomy table (one OTU per sample per row).
#subset for OTUs in the leaf core (328 OTUs, with 109 samples = 35752rows)
selected_otus_leaf <- selected_otus[selected_otus$otu %in% fungalcore_sw17,]


#remove records for the OTUs which were not phyllosphere
#(obviously sometimes core leaf OTUTs were found in soil samples too)
#Then group by samplingdate and taxonomy
#also remove rows with occ of zero (for that sampling date, that taxon was absent)

selected_otus_leaf %>%
  filter(source == 'phyllosphere') %>%
  group_by(sampling_date, final_names, Class, Order, Family, Genus) %>%
  dplyr::summarise(n=sum(abun>0)/length(abun),  
                   all=length(abun),
                   rep_ab=mean(abun),
                   sd_rep=sd(abun)
  ) %>%
  filter(n>0) -> temp

#df with stats for the whole dataset per OTU
#in other words, same thing as above, but pooling sampling dates for a given taxon.
selected_otus_leaf %>%
  filter(source == 'phyllosphere') %>%
  group_by(final_names,otu) %>%
  dplyr::summarise(
    all_ab=mean(abun),
    all_sd=sd(abun)
  ) -> temp2

#combining df and calculating the z-score for each taxon at each date
#(# of sd's that that taxon/date is above the mean for that taxon)
z_df <- left_join(temp, temp2)
z_df %>% arrange(otu) %>%
  mutate(
    z_score=(rep_ab-all_ab)/all_sd
  ) %>%
  arrange(Class, Family) -> tmp3

leaf_core <- tmp3[,c(1,14,11)]
leaf_core <- as.data.frame(leaf_core)
leaf_core_wide <- spread(leaf_core, key='sampling_date', value='z_score')


leaf_core_wide[is.na(leaf_core_wide)] <- 0
rownames(leaf_core_wide) <- leaf_core_wide$otu
leaf_core_wide$otu <-  NULL
set.seed(21)
clusters_leaf <- hclust(dist(leaf_core_wide),'complete')
memb_leaf <- cutree(clusters_leaf, k=7)
leaf_dend <- plot(clusters_leaf, main=NULL)
rect.hclust(clusters_leaf, k=7)

leaf_clusters <- data.frame(memb_leaf)
leaf_clusters$otu <-  rownames(leaf_clusters)
#assign the proper seasonality (determined after plotting leaf_core_abundance below)
leaf_clusters$stage <- 'late'
leaf_clusters$stage[leaf_clusters$memb_leaf==3] <- 'late'
leaf_clusters$stage[leaf_clusters$memb_leaf==2] <- 'mid'
leaf_clusters$stage[leaf_clusters$memb_leaf==4] <- 'late'
leaf_clusters$stage[leaf_clusters$memb_leaf==5] <- 'early'
leaf_clusters$stage[leaf_clusters$memb_leaf==6] <- 'early'
leaf_clusters$stage[leaf_clusters$memb_leaf==7] <- 'mid'


#Figure 4B
leaf_core_abundance <- data.frame(otu = as.factor(row.names(rel_otu_rare)), rel_otu_rare) %>% 
  gather(sequence_name, abun, -otu) %>%  
  left_join(map_ITS[, c('sequence_name','source', 'plant', 'sampling_date', 'year', 'sampling_week','rep')], by = 'sequence_name') %>%
  filter(source != 'soil') %>%
  filter(otu %in% fungalcore_sw17) %>%
  left_join(leaf_clusters, by='otu') %>%
  group_by(sampling_date) %>%
  mutate(sample_size = length(unique(sequence_name))) %>%
  group_by(sampling_week, sampling_date, memb_leaf, stage) %>%
  summarize(n_relabun=sum(abun)/unique(sample_size)) %>%
  filter(!is.na(memb_leaf)) %>%
  ggplot(aes(x=as.Date(sampling_date), y=n_relabun, color=as.factor(memb_leaf), group=memb_leaf))+
  geom_line(size=2,linetype = "dashed")+
  geom_point(size=4)+
  scale_color_manual(values = c('dodgerblue1','dodgerblue4','coral1','coral4','green2','green4',
                                'red'))+
  theme_classic() + theme(strip.background = element_blank(), 
                          legend.background = element_rect(fill = "transparent"),
                          legend.position = c(.4,.8), axis.text.x =element_blank(),
                          legend.key.size = unit(.4, "cm"), legend.title = element_text(size=9)) +
  labs(x=NULL, y=NULL, color='Core OTU Cluster')+
  theme(axis.title.y = element_text(size = 16))+
  ylab("Relative abundance")+
  ylim(0, 0.74)+
  labs(tag = "B")
leaf_core_abundance



#####
#Figure S2A: seasonal relative abundance at Class level
leaf_core_taxonomy_prep <- data.frame(otu = as.factor(row.names(rel_otu_rare)), rel_otu_rare) %>% 
  gather(sequence_name, abun, -otu) %>%  
  left_join(map_ITS[, c('sequence_name','source', 'plant', 'sampling_date', 'year', 'sampling_week','rep')], by = 'sequence_name') %>%
  filter(source != 'soil') %>%
  filter(otu %in% fungalcore_sw17) %>%
  left_join(leaf_clusters, by='otu') %>%
  left_join(keep_otus, by='otu') %>%
  group_by(Class, sampling_date) %>%
  summarise(n_count=sum(abun),
            n_sample=length(unique(sequence_name)),
            norm_ra=n_count/n_sample)

# Get levels and add "Unidentified"
levels <- levels(leaf_core_taxonomy_prep$Class)
levels[length(levels) + 1] <- "Unidentified"

# refactor Class to include "Unidentified" as a factor level
# 
leaf_core_taxonomy_prep$Class <- factor(leaf_core_taxonomy_prep$Class, levels = levels)
leaf_core_taxonomy_prep$Class[is.na(leaf_core_taxonomy_prep$Class)] <- "Unidentified"
  
leaf_core_taxonomy_prep$sampling_date <- gsub(pattern = " EDT", replacement = "", leaf_core_taxonomy_prep$sampling_date)
leaf_core_taxonomy_prep$sampling_Rdate <- as.Date(leaf_core_taxonomy_prep$sampling_date)

#subset dataset to remove taxa that where relative abundance not ever above 0.01.
#in other words, keep only taxa where at least one of the 7 dates has it above 0.01.
##aka, delete taxa where all 7 observations are below 0.01.

agg <- aggregate(leaf_core_taxonomy_prep$norm_ra, by = list(leaf_core_taxonomy_prep$Class), max)
agg <-  subset(agg, agg$x <0.01) 

leaf_core_taxonomy_prep <- leaf_core_taxonomy_prep[ ! leaf_core_taxonomy_prep$Class %in% agg$Group.1, ]


leaf_core_taxonomy<-  ggplot(data=leaf_core_taxonomy_prep,aes(x=sampling_Rdate, y=norm_ra, color=Class, group=Class)) +
  geom_line(size=2) +
  geom_point(size=4)+
  scale_color_manual(values = c('red','dodgerblue3','grey','green4','yellow4','coral3',
                                'greenyellow','indianred4','black'))+
  theme_classic() + theme(strip.background = element_blank(),
                          legend.position = 'bottom') +
  labs(x="Date", y="Relative abundance", color=NULL)+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"))+
  theme(axis.title.y = element_text(size = 16))+
  ylim(0, 0.99)+
  labs(tag = "A")
  
leaf_core_taxonomy
########


######
#Figure S2B: seasonal relative abundance at Phylum level

leaf_core_taxonomy_prep2 <- data.frame(otu = as.factor(row.names(rel_otu_rare)), rel_otu_rare) %>% 
  gather(sequence_name, abun, -otu) %>%  
  left_join(map_ITS[, c('sequence_name','source', 'plant', 'sampling_date', 'year', 'sampling_week','rep')], by = 'sequence_name') %>%
  filter(source != 'soil') %>%
  filter(otu %in% fungalcore_sw17) %>%
  left_join(leaf_clusters, by='otu') %>%
  left_join(keep_otus, by='otu') %>%
  group_by(Phylum, sampling_date) %>%
  summarise(n_count=sum(abun),
            n_sample=length(unique(sequence_name)),
            norm_ra=n_count/n_sample)

levels <- levels(leaf_core_taxonomy_prep2$Phylum)
levels[length(levels) + 1] <- "Unidentified"

# refactor to include "Unidentified" as a factor level
# 
leaf_core_taxonomy_prep2$Phylum <- factor(leaf_core_taxonomy_prep2$Phylum, levels = levels)
leaf_core_taxonomy_prep2$Phylum[is.na(leaf_core_taxonomy_prep2$Phylum)] <- "Unidentified"

leaf_core_taxonomy_prep2$sampling_date <- gsub(pattern = " EDT", replacement = "", leaf_core_taxonomy_prep2$sampling_date)
leaf_core_taxonomy_prep2$sampling_Rdate <- as.Date(leaf_core_taxonomy_prep2$sampling_date)

leaf_core_taxonomy2<-  ggplot(data=leaf_core_taxonomy_prep2,aes(x=sampling_Rdate, y=norm_ra, color=Phylum, group=Phylum)) +
  geom_line(size=2) +
  geom_point(size=4)+
  scale_color_manual(values = c('red','blue','grey','green4','chocolate4','magenta',
                                'greenyellow','indianred4','yellow4','violet',
                                'turquoise4','yellow','darksalmon',
                                'cadetblue','grey0'))+
  theme_classic() + theme(strip.background = element_blank(),
                          legend.position = 'bottom') +
  labs(x="Date", y="Relative abundance", color=NULL)+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"))+
  theme(axis.title.y = element_text(size = 16))+
  ylim(0, 0.99)+
  labs(tag = "B")
leaf_core_taxonomy2

#Figure S2
g <- arrangeGrob(leaf_core_taxonomy,
                 leaf_core_taxonomy2, 
                 nrow=2) #generates g
#ggsave(file="FigureS2.tiff", g, width=8, height=10, units=c("in"))





# #Contribution of the core to community dissimilarity
# leaffull <- otu_rare[,map_ITS$source == 'phyllosphere']
# leaffull <- leaffull[rowSums(leaffull)>0,]
# leafcore <- leaffull[rownames(leaffull) %in% fungalcore_sw17,]
# #unique(map_16S$sampling_week[map_16S$plant == 'switchgrass' & map_16S$year == 2017 & map_16S$source == 'phyllosphere'])
# 
# for(i in leaffull){
#   week = c(1:7)
#   leaf_BC_full <- NULL
#   for(i in week){
#     samples = map_ITS$sequence_name[map_ITS$sampling_week== i]    #selecting the unique dates for pairwise calcualtions
#     time_sub=leaffull[,colnames(leaffull) %in% samples]   #subsetting matrix to the selected dates
#     x <- apply(combn(ncol(time_sub), 2), 2, function(x) sum(abs(time_sub[,x[1]]- time_sub[,x[2]]))/2000)
#     x_names <- apply(combn(ncol(time_sub), 2), 2, function(x) paste(colnames(time_sub)[x], collapse=' - '))
#     df <- data.frame(x, x_names, i)
#     leaf_BC_full <- rbind(leaf_BC_full,df)
#   }
# }
# 
# for(i in switch17core){
#   week = c(1:7)
#   switch17_BC_core <- NULL
#   for(i in week){
#     samples = map_ITS$sequence_name[map_ITS$sampling_week== i]    #selecting the unique dates for pairwise calcualtions
#     time_sub=switch17core[,colnames(switch17core) %in% samples]   #subsetting matrix to the selected dates
#     x <- apply(combn(ncol(time_sub), 2), 2, function(x) sum(abs(time_sub[,x[1]]- time_sub[,x[2]]))/2000)
#     x_names <- apply(combn(ncol(time_sub), 2), 2, function(x) paste(colnames(time_sub)[x], collapse=' - '))
#     df <- data.frame(x, x_names, i)
#     switch17_BC_core <- rbind(switch17_BC_core,df)
#   }
# }
# 
# diffBC_switch17 <- left_join(switch17_BC_full,switch17_BC_core, by=c('x_names', 'i'))
# diffBC_switch17$diff_BC <- diffBC_switch17$x.y/diffBC_switch17$x.x
# mean(diffBC_switch17$diff_BC)
# sw17_BC_new <- switch17_BC_plot <- ggplot(diffBC_switch17, aes(x=factor(i), y=diff_BC)) +
#   geom_violin(trim=FALSE,cex=1, color='darkolivegreen3')+
#   geom_jitter(width=0.1, color='darkolivegreen3', pch=21)+
#   ylim(0,1)+
#   theme_bw()+
#   labs(x=NULL, y=NULL) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         axis.text.x = element_blank(), panel.border = element_blank())







######
#Figure 4C: seasonal relative abundance at Genus level

leaf_core_taxonomy_prep3 <- data.frame(otu = as.factor(row.names(rel_otu_rare)), rel_otu_rare) %>% 
  gather(sequence_name, abun, -otu) %>%  
  left_join(map_ITS[, c('sequence_name','source', 'plant', 'sampling_date', 'year', 'sampling_week','rep')], by = 'sequence_name') %>%
  filter(source != 'soil') %>%
  filter(otu %in% fungalcore_sw17) %>%
  left_join(leaf_clusters, by='otu') %>%
  left_join(keep_otus, by='otu') %>%
  group_by(Genus, sampling_date) %>%
  summarise(n_count=sum(abun),
            n_sample=length(unique(sequence_name)),
            norm_ra=n_count/n_sample)

levels <- levels(leaf_core_taxonomy_prep3$Genus)
levels[length(levels) + 1] <- "Unidentified"

# refactor Species to include "Unidentified" as a factor level
# 
leaf_core_taxonomy_prep3$Genus <- factor(leaf_core_taxonomy_prep3$Genus, levels = levels)
leaf_core_taxonomy_prep3$Genus[is.na(leaf_core_taxonomy_prep3$Genus)] <- "Unidentified"

leaf_core_taxonomy_prep3$sampling_date <- gsub(pattern = " EDT", replacement = "", leaf_core_taxonomy_prep3$sampling_date)
leaf_core_taxonomy_prep3$sampling_Rdate <- as.Date(leaf_core_taxonomy_prep3$sampling_date)



#subset dataset to remove taxa that where relative abundance not ever above 0.01.
#in other words, keep only taxa where at least one of the 7 dates has it above 0.01.
##aka, delete taxa where all 7 observations are below 0.01.

agg <- aggregate(leaf_core_taxonomy_prep3$norm_ra, by = list(leaf_core_taxonomy_prep3$Genus), max)
agg <-  subset(agg, agg$x <0.01) 

leaf_core_taxonomy_prep3 <- leaf_core_taxonomy_prep3[ ! leaf_core_taxonomy_prep3$Genus %in% agg$Group.1, ]



leaf_core_taxonomy3<-  ggplot(data=leaf_core_taxonomy_prep3,aes(x=sampling_Rdate, y=norm_ra, color=Genus, group=Genus)) +
  geom_line(size=2) +
  geom_point(size=4)+
  scale_color_manual(values = c('red','blue','grey','green4','chocolate4','magenta',
                                'greenyellow','indianred3','yellow4','violet',
                                'turquoise4','yellow','darksalmon',
                                'grey0'))+
  theme_classic() + theme(strip.background = element_blank(),
                          legend.position = 'bottom') +
  labs(x="Date", y="Relative abundance", color=NULL)+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"))+
  theme(axis.title.y = element_text(size = 16))+
  ylim(0, 0.74)+
  labs(tag = "C")
leaf_core_taxonomy3



#Figure 4
g <- arrangeGrob(leaf_occ_abun_plot,
                 leaf_core_abundance,
                 leaf_core_taxonomy3, 
                 nrow=3) 
#ggsave(file="Figure4.tiff", g, scale=2, width=5, height=10, units=c("in"))





#########################
#### Alpha Diversity ####
#########################

#Subset rarefied OTU table to only soil samples
soil_otu <- otu_rare[,map_ITS$source=="soil"]
map.soil<- map_ITS[map_ITS$source=="soil",]

#soil evenness and diversity
s <- specnumber(soil_otu,MARGIN=2)
h <- diversity(t(soil_otu), index = "shannon")
mean(h)
sd(h)
pielou<-h/log(s)
mean(pielou)
sd(pielou)


#leaf evenness and diversity
s <- specnumber(leaf_otu,MARGIN=2)
h <- diversity(t(leaf_otu), index = "shannon")
mean(h)
sd(h)
pielou<-h/log(s)
mean(pielou)
sd(pielou)

#for all samples (soil and phyllosphere)
s <- specnumber(otu_rare,MARGIN=2)
h <- diversity(t(otu_rare), index = "shannon")
pielou=h/log(s)
map_ITS$sampling_date <- gsub(pattern = " EDT", replacement = "", map_ITS$sampling_date)
map_ITS$sampling_Rdate <- as.Date(map_ITS$sampling_date)

### Setting up Contextual Data Maps
map.div <- map_ITS
map.div$Richness <- s
map.div$Shannon <- h
map.div$Pielou <- pielou 
map.div.plant <- map.div[map.div$source=="phyllosphere",]

cor.test(map.div.plant$time_numeric, map.div.plant$Richness)
cor.test(map.div.plant$time_numeric, map.div.plant$Shannon)
cor.test(map.div.plant$time_numeric, map.div.plant$Pielou)

### Phyllosphere switchgrass unique times
switch_unique_times <- unique(map.div.plant$time_numeric)[order(unique(map.div.plant$time_numeric))]

# Look at species accumulation in the phyllosphere
switch_accumulation <- rep(1, length(switch_unique_times))
z <- NULL
for( i in 1:length(unique(map.div.plant$time_numeric))){
  x <- leaf_otu[,map.div.plant$time_numeric==switch_unique_times[i]]
  z <- c(z, row.names(x[rowSums(x)!=0,]))
  switch_accumulation[i] <- length(unique(z))
}

Spec_accum <- NULL
Spec_accum$Species_Accumulation <- switch_accumulation
Spec_accum$Date <- unique(map.div.plant$sampling_Rdate)[order(unique(map.div.plant$sampling_Rdate))]

Spec_accum <- as.data.frame(Spec_accum)







#Figure 3A: Phyllosphere richness and accumulation through the season
Fig3A <- ggplot(Spec_accum, aes(x=Date, y=Species_Accumulation)) + 
  geom_point(size=3) +
  geom_line()+
  theme_bw()+
  labs(y="Total Observed Taxa") + 
  scale_y_continuous(breaks = seq(1000, 4000, by = 500))+
  #scale_color_manual(values=c("green4"))  + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"))+
  guides(shape=FALSE)+
  labs(tag = "A")
  
Fig3A

#Figure 3B: phyllosphere richness through growing season
Fig3B<-ggplot(map.div.plant, aes(x=sampling_Rdate, y=Richness,group=sampling_Rdate)) + 
  geom_boxplot() +
  theme_bw()+
  labs(y="Richness") + 
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
  guides(shape=FALSE)+
  theme(axis.title.y = element_text(size = 16))+
  labs(tag = "B")
Fig3B

#Figure 3C: phyllosphere Shannon diversity through growing season
Fig3C<-ggplot(map.div.plant, aes(x=sampling_Rdate, y=Shannon,group=sampling_Rdate)) + 
  geom_boxplot() +
  theme_bw()+
  labs(y="Shannon Diversity",x="Date") + 
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"))+
  #theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
  guides(shape=FALSE)+
  theme(axis.title.y = element_text(size = 16))+
  theme(axis.title.x = element_text(size = 16))+
  labs(tag = "C")
Fig3C

#Figure 3D: Pielou's evenness through growing season
Fig3D<-ggplot(map.div.plant, aes(x=sampling_Rdate, y=Pielou,group=sampling_Rdate)) + 
  geom_boxplot() +
  theme_bw()+
  labs(y="Pielou's Evenness",x="Date") + 
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"))+
  guides(shape=FALSE)+
  theme(axis.title.y = element_text(size = 16))+
  theme(axis.title.x = element_text(size = 16))+
  labs(tag = "D")
Fig3D

f <- arrangeGrob(Fig3A, Fig3B, Fig3C, Fig3D, nrow=2)

#ggsave("Figure3.tiff", f, width=9, height=8, units="in")


################################
#### Table 1 Analyses####
################################

#Table 1 data for soil: mean number of OTUs per sample
rare_OTU_pa <- as.data.frame(1*((soil_otu>0)==1))
mean(colSums(rare_OTU_pa))
sd(colSums(rare_OTU_pa))
ncol(soil_otu)
#list of otus which appear in soil
soilotus <- as.data.frame(row.names(soil_otu[which(rowSums(soil_otu) > 0),]))
colnames(soilotus)<- c("SoilOTUs")

#Table 1 data for phyllosphere: mean number of OTUs per sample
rare_OTU_pa <- as.data.frame(1*((leaf_otu>0)==1))
mean(colSums(rare_OTU_pa))
sd(colSums(rare_OTU_pa))
ncol(leaf_otu)
#list of otus which appear in phyllosphere
leafotus <- as.data.frame(row.names(leaf_otu[which(rowSums(leaf_otu) > 0),]))
colnames(leafotus)<- c("LeafOTUs")

#For Table 1: unique soil otus
uniquesoilotus <- as.data.frame(soilotus[!soilotus$SoilOTUs %in% leafotus$LeafOTUs,])

#For Table 1: unique leaf otus
uniqueleafotus <- as.data.frame(leafotus[!leafotus$LeafOTUs %in% soilotus$SoilOTUs,])


################################
#### Beta Diversity Analyses####
################################

#create a 'week number' column and add to metadata
Date_to_Week <- data.frame(sampling_date=unique(map_ITS$sampling_date)[order(unique(map_ITS$sampling_date))], Week=c(1:8))
map_ITS <- left_join(map_ITS, Date_to_Week, by="sampling_date")


#Generate Bray-Curtis Distance matrices for full dataset, leaf dataset, and soil dataset.
dist.otu <- vegdist(t(otu_rare), method="bray")
dist.leaf <- vegdist(t(leaf_otu),method="bray")
dist.soil <- vegdist(t(soil_otu), method="bray")

#test for community difference between phyllosphere and soil
set.seed(1)
adonis(dist.otu~map_ITS$source)
#F=90.06, R2=0.46, P=0.001.

#Leaf: test for community differences btween sampling dates, and between fertilization treatments.
set.seed(1)
adonis(dist.leaf~map_leaf$time_numeric*map_leaf$treatment)
set.seed(1)
adonis(dist.leaf~map_leaf$treatment*map_leaf$time_numeric)

#Soil: test for community differences btween sampling dates, and between fertilization treatments.
set.seed(1)
adonis(dist.soil~map.soil$time_numeric*map.soil$treatment)
set.seed(1)
adonis(dist.soil~map.soil$treatment*map.soil$time_numeric)


#######################################
### Figure 2B: Plant PCoA ###
#######################################

### Start by getting the dates for plant smaples
sg.dates <- unique(map_leaf$sampling_date)[order(unique(map_leaf$sampling_date))]

### Set up the points from the "entire" plant PCoA
pcoa.leaf <- cmdscale(dist.leaf, eig=TRUE)
leaf.points <- pcoa.leaf$points

### determine the average location on the PCoA at each timepoint
leaf.points.collapsed <- NULL
leaf.samples.per.date <- NULL
for (i in 1:length(sg.dates)){
  x <- leaf.points[map_leaf$sampling_date==sg.dates[i],]
  leaf.samples.per.date <- c(leaf.samples.per.date, nrow(x))
  leaf.points.collapsed <- rbind(leaf.points.collapsed, c(colSums(x), sd(x[,1]), sd(x[,2])))
}
leaf.points.collapsed <- as.data.frame(leaf.points.collapsed)
colnames(leaf.points.collapsed) <- c("Axis1", "Axis2", "sd_axis1", "sd_axis2")
leaf.points.collapsed$sampling_date <- sg.dates
leaf.points.collapsed[,1:2] <- leaf.points.collapsed[,1:2]/leaf.samples.per.date


leaf.points.collapsed$sampling_date <-as.factor(leaf.points.collapsed$sampling_date)
plant.points.collapsed <- left_join(leaf.points.collapsed, Date_to_Week, by="sampling_date")
plant.points.collapsed$WeekNumeric <- as.numeric(plant.points.collapsed$Week)

### Set up the weather maps for Weather Environmental fits
switch.weather <- map_leaf[,c("precipitation", "Air_temp_mean", "air_temp_max", "Air_Temp_Min", "Air_Pressure", "RH", "AH", "Wind_Speed_Mean", "Solar_Radiation", "PAR", "soil_temp_5_cm_bare_avg", "sampling_date", "time_numeric")]
switch.weather <- unique(switch.weather)

scale_arrow <- function(arrows, data, at = c(0, 0), fill = 0.75) {
  u <- c(range(data[,1], range(data[,2])))
  u <- u - rep(at, each = 2)
  r <- c(range(arrows[, 1], na.rm = TRUE), range(arrows[, 2], na.rm = TRUE))
  rev <- sign(diff(u))[-2]
  
  if (rev[1] < 0) {
    u[1:2] <- u[2:1]
  }
  if (rev[2] < 0) {
    u[3:4] <- u[4:3]
  }
  u <- u/r
  u <- u[is.finite(u) & u > 0]
  invisible(fill * min(u))
}


#arrow_scaling <- scale_arrow(data = scores(pcoa.leaf, display = "sites"), arrows = scores(envfit.phyllo.sc, display = "vectors"))

#arrow_scaling <- scale_arrow(data = scores(pcoa.leaf, display = "sites"), arrows = rbind(scores(envfit.phyllo.sc, display="vectors"), scores(envfit.phyllo.lc, display="vectors")))



### Set up EnvFit for ggplot (remembering to use vegan:::ordiArrowMul to adjust arrow sizes)

switch.weather <- switch.weather[order(switch.weather$sampling_date),]
switch.weather$sampling_date == plant.points.collapsed$sampling_date
switch.weather$Week <- plant.points.collapsed$WeekNumeric

envfit.phyllo.weather <- envfit(plant.points.collapsed[,1:2], switch.weather[,c(1:11,14)])
envfit.phyllo.weather.df<-as.data.frame(scores(envfit.phyllo.weather, display = "vectors"))
colnames(envfit.phyllo.weather.df) <- c("Dim1", "Dim2")

phyllo.leaf.chemistry <- map_leaf[,c("LDMC_mg_per_g", "nitrogen_percent", "carbon_percent", "carbon_per_nitrogen", "height_mean_cm")]
envfit.phyllo.lc <- envfit(pcoa.leaf$points, phyllo.leaf.chemistry)
envfit.phyllo.lc.df<-as.data.frame(scores(envfit.phyllo.lc, display = "vectors"))

phyllo.soil.chemistry <- map_leaf[,c("pH", "P_ppm", "K_ppm", "Ca_ppm", "Mg_ppm", "organic_matter", "NO3N_ppm", "NH4_ppm", "soil_moisture_percent", "soil_temp_10cm")]
envfit.phyllo.sc <- envfit(pcoa.leaf$points, phyllo.soil.chemistry)
envfit.phyllo.sc.df<-as.data.frame(scores(envfit.phyllo.sc, display = "vectors"))

envfit.phyllo.total <- rbind(envfit.phyllo.weather.df, envfit.phyllo.lc.df, envfit.phyllo.sc.df)

envfit.phyllo.total$r2 <- c(envfit.phyllo.weather$vectors$r, envfit.phyllo.lc$vectors$r, envfit.phyllo.sc$vectors$r)
envfit.phyllo.total$pval <- c(envfit.phyllo.weather$vectors$pvals, envfit.phyllo.lc$vectors$pvals, envfit.phyllo.sc$vectors$pvals)

arrow_scaling <- scale_arrow(arrows=envfit.phyllo.total[,1:2], data=plant.points.collapsed)

envfit.phyllo.total$Axis1 <- envfit.phyllo.total$Dim1 * arrow_scaling
envfit.phyllo.total$Axis2 <- envfit.phyllo.total$Dim2 * arrow_scaling
envfit.phyllo.total$Variable <- row.names(envfit.phyllo.total)

### Subset to p <0.05 & Rsquared >0.4
envfit.phyllo.sub <- subset(envfit.phyllo.total, envfit.phyllo.total$r2 >0.3999) 
envfit.phyllo.sub <- subset (envfit.phyllo.sub , envfit.phyllo.sub$pval<0.05)

#Write envfit results to Table S3
#write.table(file = "TableS3_EnvFit.txt", x=envfit.phyllo.total[,1:4], sep="\t", quote=FALSE)

### Set up some of the plotting specifics
Point_Sizes <- seq(from=2, to=6, length.out = 8)
Ax1.plant <- pcoa.leaf$eig[1]/sum(pcoa.leaf$eig)
Ax2.plant <- pcoa.leaf$eig[2]/sum(pcoa.leaf$eig)

#Abbreviate Variable names so easier to read in figure
envfit.phyllo.sub[1,7]<-c("MinTemp")
envfit.phyllo.sub[4,7]<-c("LDMC")
envfit.phyllo.sub[5,7]<-c("NLeaf")
envfit.phyllo.sub[6,7]<-c("C:N")
envfit.phyllo.sub[7,7]<-c("Hgt")

###Figure 2B: Plant PCoA
Fig2B <- ggplot(plant.points.collapsed, aes(x=Axis1, y=Axis2))+
  #coord_fixed()+
  coord_cartesian(xlim=c(-.5,0.5),ylim=c(-0.4,0.4))+
  scale_colour_identity() +
  geom_segment(data=plant.points.collapsed, aes(x=Axis1,xend=Axis1+sd_axis1,y=Axis2,yend=Axis2,color="grey"))+
  geom_segment(data=plant.points.collapsed, aes(x=Axis1,xend=Axis1-sd_axis1,y=Axis2,yend=Axis2,color="grey"))+
  geom_segment(data=plant.points.collapsed, aes(x=Axis1,xend=Axis1,y=Axis2,yend=Axis2+sd_axis2,color="grey"))+
  geom_segment(data=plant.points.collapsed, aes(x=Axis1,xend=Axis1,y=Axis2,yend=Axis2-sd_axis2,color="grey"))+
  xlab(label = paste(round(Ax1.plant,digits = 3)*100, "% Var. Explained", sep = ""))+
  ylab(label= paste(round(Ax2.plant,digits = 3)*100, "% Var. Explained", sep = ""))+
  geom_segment(data = envfit.phyllo.sub,
               aes(x = 0, xend = Axis1, y = 0, yend = Axis2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey0")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  geom_point(aes(size=Point_Sizes[2:8]),color='grey0',pch=22,fill="green4")+
  geom_text(data=envfit.phyllo.sub, aes(label=Variable),nudge_x=-.02)+
  labs(tags="B")

Fig2B


#######################################
### Figure 2A: Soil and Plant PCoA ###
#######################################

### Make Whole PCoA
pcoa.whole <- cmdscale(dist.otu, eig=TRUE)

### Set up dataframe for plotting
points.whole <- pcoa.whole$points
points.whole <- as.data.frame(points.whole)
colnames(points.whole) <- c("Axis1", "Axis2")
points.whole$source <- factor(map_ITS$source, levels=c("soil", "phyllosphere"))
points.whole$sampling_date <- map_ITS$sampling_date
points.whole <- left_join(points.whole, Date_to_Week, by="sampling_date")
points.whole$Week <- factor(points.whole$Week, levels=c(1:10))

points.whole$SourceFert <- map_ITS$treatment
points.whole$SourceFert <- factor(points.whole$SourceFert, levels = c("standard fertilization", "nitrogen free"))
points.whole$SampleType <- points.whole$SourceFert


### Determine % variation explained on each axis
Ax1.whole <- pcoa.whole$eig[1]/sum(pcoa.whole$eig)
Ax2.whole <- pcoa.whole$eig[2]/sum(pcoa.whole$eig)

###Figure 2A: plant and soil PCoA
Fig2A <- ggplot(points.whole, aes(x=Axis1, y=Axis2,color=source,fill=source))+
  #coord_fixed()+
  coord_cartesian(xlim=c(-.5,0.5),ylim=c(-0.5,0.5))+
  geom_point(aes(shape=SampleType, size=Week))+
  scale_shape_manual(values = c(21,25))+
  scale_size_manual(values=Point_Sizes)+
  scale_color_manual(values=c("grey0","grey0"))+
  scale_fill_manual(values=c("chocolate4","green4"))+
  xlab(label = paste(round(Ax1.whole, 3)*100, "% Var. Explained", sep=""))+
  ylab(label= paste(round(Ax2.whole, 3)*100, "% Var. Explained", sep=""))+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(tags="A")
Fig2A

#Produce legend for Figure 2
Fig2A_leg <- ggplot(points.whole, aes(x=Axis1, y=Axis2,color=source,fill=source))+
  coord_fixed()+
  geom_point(aes(shape=SampleType, size=sampling_date))+
  scale_shape_manual(values = c(21, 25))+
  scale_size_manual(values=Point_Sizes)+
  scale_color_manual(values=c("grey0","grey0"))+
  scale_fill_manual(values=c("chocolate4","green4"))+
  xlab(label = paste(round(Ax1.whole, 3)*100, "% Var. Explained", sep=""))+
  ylab(label= paste(round(Ax2.whole, 3)*100, "% Var. Explained", sep=""))+
  theme(legend.position = "right", legend.box = "horizontal")+
  guides(fill=guide_legend(override.aes=list(shape=22, color=c(soil="grey0",phyllosphere="grey0"),fill=c(soil="chocolate4",phyllosphere="green4"), size=5)), shape=guide_legend(override.aes = list(size=5)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
Fig2_Legend <- get_legend(Fig2A_leg)



#######################################
### Figure 2C: Soil PCoA ###
#######################################

### Start by getting the dates for soil smaples
sg.dates <- unique(map.soil$sampling_date)[order(unique(map.soil$sampling_date))]

### Set up the points from the "entire" plant PCoA
pcoa.soil <- cmdscale(dist.soil, eig=TRUE)
soil.points <- pcoa.soil$points

### determine the average location on the PCoA at each timepoint
soil.points.collapsed <- NULL
soil.samples.per.date <- NULL
for (i in 1:length(sg.dates)){
  x <- soil.points[map.soil$sampling_date==sg.dates[i],]
  soil.samples.per.date <- c(soil.samples.per.date, nrow(x))
  soil.points.collapsed <- rbind(soil.points.collapsed, c(colSums(x), sd(x[,1]), sd(x[,2])))
}
soil.points.collapsed <- as.data.frame(soil.points.collapsed)
colnames(soil.points.collapsed) <- c("Axis1", "Axis2", "sd_axis1", "sd_axis2")
soil.points.collapsed$sampling_date <- sg.dates
soil.points.collapsed[,1:2] <- soil.points.collapsed[,1:2]/soil.samples.per.date


soil.points.collapsed$sampling_date <-as.factor(soil.points.collapsed$sampling_date)
soil.points.collapsed <- left_join(soil.points.collapsed, Date_to_Week, by="sampling_date")
soil.points.collapsed$WeekNumeric <- as.numeric(soil.points.collapsed$Week)

### Set up the weather maps for Weather Environmental fits
switch.weather <- map.soil[,c("precipitation", "Air_temp_mean", "air_temp_max", "Air_Temp_Min", "Air_Pressure", "RH", "AH", "Wind_Speed_Mean", "Solar_Radiation", "PAR", "soil_temp_5_cm_bare_avg", "sampling_date", "time_numeric")]
switch.weather <- unique(switch.weather)

### Set up EnvFit for ggplot (remembering to use vegan:::ordiArrowMul to adjust arrow sizes)

switch.weather <- switch.weather[order(switch.weather$sampling_date),]
switch.weather$sampling_date == soil.points.collapsed$sampling_date
switch.weather$Week <- soil.points.collapsed$WeekNumeric

envfit.phyllo.weather <- envfit(soil.points.collapsed[,1:2], switch.weather[,c(1:11,14)])
envfit.phyllo.weather.df<-as.data.frame(scores(envfit.phyllo.weather, display = "vectors"))
colnames(envfit.phyllo.weather.df) <- c("Dim1", "Dim2")

phyllo.leaf.chemistry <- map.soil[,c("LDMC_mg_per_g", "nitrogen_percent", "carbon_percent", "carbon_per_nitrogen", "height_mean_cm")]
envfit.phyllo.lc <- envfit(pcoa.soil$points, phyllo.leaf.chemistry, na.rm=TRUE)
envfit.phyllo.lc.df<-as.data.frame(scores(envfit.phyllo.lc, display = "vectors"))

phyllo.soil.chemistry <- map.soil[,c("pH", "P_ppm", "K_ppm", "Ca_ppm", "Mg_ppm", "organic_matter", "NO3N_ppm", "NH4_ppm", "soil_moisture_percent", "soil_temp_10cm")]
envfit.phyllo.sc <- envfit(pcoa.soil$points, phyllo.soil.chemistry)
envfit.phyllo.sc.df<-as.data.frame(scores(envfit.phyllo.sc, display = "vectors"))

envfit.phyllo.total <- rbind(envfit.phyllo.weather.df, envfit.phyllo.lc.df, envfit.phyllo.sc.df)

envfit.phyllo.total$r2 <- c(envfit.phyllo.weather$vectors$r, envfit.phyllo.lc$vectors$r, envfit.phyllo.sc$vectors$r)
envfit.phyllo.total$pval <- c(envfit.phyllo.weather$vectors$pvals, envfit.phyllo.lc$vectors$pvals, envfit.phyllo.sc$vectors$pvals)

arrow_scaling <- scale_arrow(arrows=envfit.phyllo.total[,1:2], data=soil.points.collapsed)

envfit.phyllo.total$Axis1 <- envfit.phyllo.total$Dim1 * arrow_scaling
envfit.phyllo.total$Axis2 <- envfit.phyllo.total$Dim2 * arrow_scaling
envfit.phyllo.total$Variable <- row.names(envfit.phyllo.total)

### Subset to p <0.05 & Rsquared >0.4
envfit.phyllo.sub <- subset(envfit.phyllo.total, envfit.phyllo.total$r2 >0.3999) 
envfit.phyllo.sub <- subset (envfit.phyllo.sub , envfit.phyllo.sub$pval<0.05)

#Write table of envfit results (will manually add to Table S3 along with leaf results)
#write.table(file = "TableS4_EnvFit.txt", x=envfit.phyllo.total[,1:4], sep="\t", quote=FALSE)

### Set up some of the plotting specifics
Point_Sizes <- seq(from=2, to=6, length.out = 8)
Ax1.soil <- pcoa.soil$eig[1]/sum(pcoa.soil$eig)
Ax2.soil <- pcoa.soil$eig[2]/sum(pcoa.soil$eig)

#Abbreviate Variable names so easier to read in figure

envfit.phyllo.sub[1,7]<-c("Mean Wind Speed")



###Figure 2C: soil PCoA
Fig2C <- ggplot(soil.points.collapsed, aes(x=Axis1, y=Axis2))+
  #coord_fixed()+
  coord_cartesian(xlim=c(-.25,0.25),ylim=c(-0.2,0.2))+
  scale_colour_identity() +
  geom_segment(data=soil.points.collapsed, aes(x=Axis1,xend=Axis1+sd_axis1,y=Axis2,yend=Axis2,color="grey"))+
  geom_segment(data=soil.points.collapsed, aes(x=Axis1,xend=Axis1-sd_axis1,y=Axis2,yend=Axis2,color="grey"))+
  geom_segment(data=soil.points.collapsed, aes(x=Axis1,xend=Axis1,y=Axis2,yend=Axis2+sd_axis2,color="grey"))+
  geom_segment(data=soil.points.collapsed, aes(x=Axis1,xend=Axis1,y=Axis2,yend=Axis2-sd_axis2,color="grey"))+
  xlab(label = paste(round(Ax1.soil,digits = 3)*100, "% Var. Explained", sep = ""))+
  ylab(label= paste(round(Ax2.soil,digits = 3)*100, "% Var. Explained", sep = ""))+
  geom_segment(data = envfit.phyllo.sub,
               aes(x = 0, xend = Axis1, y = 0, yend = Axis2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey0")+
  geom_text(data=envfit.phyllo.sub, aes(label=Variable),nudge_x=-0.02)+
  geom_point(aes(size=Point_Sizes[1:8]),color='grey0',pch=22,fill="chocolate4") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  labs(tags="C")

Fig2C

Fig2 <- grid.arrange(Fig2A,Fig2B,Fig2_Legend,Fig2C,nrow=2)

#ggsave(file="Figure2.tiff", Fig2, width=11, height=9, units=c("in"))



####################################################
### Figure S1:  betadispersion through time ###
####################################################

#calculate betadispersion for leaf samples at each sampling date
Dispersion.switch <- betadisper(dist.leaf, map_leaf$sampling_date)

names(Dispersion.switch$distances)==map_leaf$sequence_name
Dispersion.switch.df <- data.frame(Distance_to_Median=Dispersion.switch$distances, Date=map_leaf$sampling_date)

Dispersion.switch.df$Date <- as.factor(Dispersion.switch.df$Date)

#transform for ANOVA
Dispersion.switch.df$Distance_to_Median_trans <- log10(Dispersion.switch.df$Distance_to_Median)
#Do ANOVA 
options(contrasts=c(unordered='contr.sum', ordered='contr.poly'))
betadisp_fit <- aov(Distance_to_Median_trans ~ Date, data=Dispersion.switch.df)
drop1(betadisp_fit,~.,test="F") # type III SS and F Tests

#TESTING ASSUMPTIONS
betadisp_resids <- residuals(betadisp_fit)
betadisp_preds <- predict(betadisp_fit)
plot(betadisp_resids ~ betadisp_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#test for homogeneity of variance
bartlett.test(Distance_to_Median_trans ~ Date, data=Dispersion.switch.df)
leveneTest(Distance_to_Median_trans ~ Date, data=Dispersion.switch.df)
#Test for normality
shapiro.test(betadisp_resids)
qqnorm(betadisp_resids)
qqline(betadisp_resids)
TukeyHSD(betadisp_fit) 
  
#Convert date format for plotting
Dispersion.switch.df$sampling_Rdate <- as.Date(Dispersion.switch.df$Date)

#Figure S1: beta dispersion through time
FigS1 <- ggplot(Dispersion.switch.df, aes(x=sampling_Rdate, y=Distance_to_Median,group=Date))+
  geom_boxplot(fill='green4')+
  geom_point()+
  ylab("Beta Dispersion")+
  xlab("Date")+
  theme(axis.text.x = element_text(angle = 60))+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"))+
  ylim(c(0,1))+
  theme_bw()

FigS1

#Figure S1
#ggsave("FigureS1.tiff", FigS1, width = 5, height = 4, units = "in")



######################
# Data preparation for MENA analysis (which will be used for Figure 5: network)
# Before continuing, you must run "RecreateGrady2019_OTUtable.R". It produces 3 files:
#1)otu_rare_16s.txt: the OTU table from Grady et al 2019 (rarefied to 1000 reads)
#2)core_16s.txt: the core phyllosphere OTUs from Grady et al. 2019
#3)tax_16s.txt: the taxonomy of the core phyllosphere OTUs from Grady et al. 2019
######################

#ITS:  subset the rarefied OTU table for core leaf otus
rare_core_fungal <- as.data.frame(subset(leaf_otu, rownames(leaf_otu) %in% fungalcore_sw17))

#16S: Upload the rarefied OTU table from Grady et al. 2019 code. 
otu_rare_16s <- read.table("otu_rare_16s.txt",header=TRUE, sep='\t',row.names=1)
#upload switchgrass 2017 core OTU names and format as vector
core_16s <-read.table("core_16s.txt",sep='\t')
core_16s <- as.vector(core_16s$V1)
#subset 16s OTU table to only the samples that are also in the fungal data.
otu_rare_16s<- otu_rare_16s[,colnames(otu_rare_16s)%in%colnames(rare_core_fungal)] 

#subset the rarefied 16s OTU table for 16s core leaf otus
rare_core_16s <- as.data.frame(subset(otu_rare_16s, rownames(otu_rare_16s) %in% core_16s))
#calculate minimum number of samples in which each OTU appears.
bactarch_otu_PA <- as.data.frame(1*((rare_core_16s>0)==1))
bactarch_otu_PA$rowsum <- rowSums(bactarch_otu_PA)
min(bactarch_otu_PA$rowsum)

#subset rarefied fungal phyllosphere table for only samples also in the bactarch data.
rare_core_fungal<- rare_core_fungal[,colnames(rare_core_fungal)%in%colnames(rare_core_16s)] 
#calculate number of samples in which each OTU appears.
fungal_otu_PA <- as.data.frame(1*((rare_core_fungal>0)==1))
fungal_otu_PA$rowsum <- rowSums(fungal_otu_PA)
min(fungal_otu_PA$rowsum)

# Order the samples in the fungal phyllosphere core otu table
rare_core_fungal <- rare_core_fungal[,order(colnames(rare_core_fungal))]
# Order the 16s phyllosphere core otu table the same way
rare_core_16s <- rare_core_16s[,order(colnames(rare_core_16s))]
# check that samples are in the same order in both tables
colnames(rare_core_fungal) == colnames(rare_core_16s)
#bind the otu tables
rare_core_all <- rbind(rare_core_fungal, rare_core_16s)

#replace zeros with an empty space (for input to MENA analysis)
#zeroes meaning not detectable.
rare_core_all[rare_core_all == 0] <- ''

#save as .txt file.
#write.table(rare_core_all, file="rare_core_all.txt", sep = '\t', quote=FALSE, row.names=TRUE, col.names=TRUE)
#manually add a tab as the first character of that file. then conduct
#Molecular Ecological Network Analysis as described in manuscript.

#loaded it on 20190905.

#The number of nodes: 239
#The number of links: 359
#The average path: 4.417
#R square of power-law: 0.876
#Average Connectivity : 3.00418410041841
#Average length : 4.417
#R square of power-law: 0.876
#Power-law: lamda = 1.410943, R2 of power fitting= 0.9762555



#######
## Calculate various results from MENA analysis
## Create data needed for Cytoscape to produce Figure 5 (network)
#######

#upload initial network nodes from MENA and add metadata.
initialnodes <- read.csv("initialuploadnodes.csv",header=TRUE, sep=',')

#get list of fungal OTUs, and include Fungi as kingdom name
fungal_tax <- subset(keep_otus, select=c(Kingdom, otu))
#get list of 16s OTUs, and include Kingdom name (all are bacteria, no archaea)
tax_16s <- read.table("tax_16s.txt", header=TRUE, sep='\t')
bactarch_tax <- subset(tax_16s, select=c(Kingdom, otu))




nodes_withtax <- merge(initialnodes,bactarch_tax, by.x='shared.name', by.y='otu',all.x=TRUE)
nodes_withtax <- merge(nodes_withtax, fungal_tax, by.x='shared.name', by.y='otu',all.x=TRUE)
#combine the two columns
nodes_withtax$Kingdom.x <-as.character(nodes_withtax$Kingdom.x)
nodes_withtax$Kingdom.x <- replace_na(nodes_withtax$Kingdom.x, rep('Fungi'))

#clean up the table
nodes_withtax <- subset(nodes_withtax, select=-c(Kingdom.y))
nodes_withtax$Kingdom.x <- gsub(pattern = "k:", replacement = "", x = nodes_withtax$Kingdom.x )

#find taxonomy of the most highly connected fungal OTUs in the microbial network
subset(keep_otus, keep_otus$otu=="OTU_4280")
subset(keep_otus, keep_otus$otu=="OTU_4321")
subset(keep_otus, keep_otus$otu=="OTU_3081")

#find taxonomy of the most highly connected bacterial OTU in the microbial network
subset(tax_16s, tax_16s$otu=="OTU632")


#upload the node attribute file (which was an output from MENA)
node_degree <- read.table("MENA_Node_attrib_file.txt", sep='\t',header=TRUE)
node_degree <- subset(node_degree, select=c(Name, node.degree))

nodes_withtax <- merge(nodes_withtax,node_degree, by.x='shared.name',by.y='Name',all.x=TRUE)

colnames(nodes_withtax)[colnames(nodes_withtax) == 'Kingdom.x'] <- 'Kingdom'


#Write table to text file: will be used by Cytoscape to prepare Figure 5.
#write.csv(nodes_withtax, file="nodes_withtax.csv", quote=FALSE, sep=',', row.names=FALSE)



#Read in the MENA edge file. Will inform on number and type of connections between taxa.
mena_int<- read.csv('MENA_Edge_attrib_file.txt', sep='\t', header=F)

#Get list of network nodes (OTUs)
coreOTUs <- nodes_withtax$name

#parse apart the Mena Int file into a dataframe.
interTAXA <- mena_int %>% separate(V1, c('start','connect', 'end','value'), sep='([\\(\\)\\=])', remove=T) %>%
  
  mutate(start=str_trim(start, side = "both"),
         end=str_trim(end, side = "both"),
         start.memb=if_else(start %in% coreOTUs, 1,0),
         end.memb=if_else(end %in% coreOTUs, 1, 0),
         core.int=start.memb+end.memb,
  
         domain.start= if_else(start %in% tax_16s$otu, '16S', 'ITS'),
         domain.end= if_else(end %in% tax_16s$otu, '16S', 'ITS'),
         fungus=if_else(start %in% fungal_tax$otu, 1,0),
         fungus=if_else(end %in% fungal_tax$otu,1,fungus),
         bacteria= if_else(start %in% tax_16s$otu, 1, 0),
         bacteria= if_else(end %in% tax_16s$otu, 1, bacteria),
         bact_fung=if_else(bacteria+fungus==2, 1, 0))

#number of fungal OTUs in the network
nrow(nodes_withtax[nodes_withtax$Kingdom =="Fungi",])
#number of bacterial OTUs in the network
nrow(nodes_withtax[nodes_withtax$Kingdom =="Bacteria",])

#total number of bacteria-fungi interactions
sum(interTAXA$bact_fung)
#or
nrow(interTAXA[interTAXA$bact_fung == "1",])
#positive bacteria-fungi interactions
nrow(interTAXA[interTAXA$bact_fung == "1" & interTAXA$connect=="pp",])
#negative bacteria-fungi interactions
nrow(interTAXA[interTAXA$bact_fung == "1" & interTAXA$connect=="np",])

#total number of fungi-fungi interactions
nrow(interTAXA[interTAXA$fungus == "1" & interTAXA$bacteria=="0",])
#positive fungi-fungi interactions
nrow(interTAXA[interTAXA$fungus == "1" & interTAXA$bacteria=="0" & interTAXA$connect=="pp",])
#negative fungi-fungi interactions
nrow(interTAXA[interTAXA$fungus == "1" & interTAXA$bacteria=="0" & interTAXA$connect=="np",])


#total number of bacteria-bacteria interactions
nrow(interTAXA[interTAXA$fungus == "0" & interTAXA$bacteria=="1",])
#positive bacteria-bacteria interactions
nrow(interTAXA[interTAXA$fungus == "0" & interTAXA$bacteria=="1" & interTAXA$connect=="pp",])
#negative bacteria-bacteria interactions
nrow(interTAXA[interTAXA$fungus == "0" & interTAXA$bacteria=="1" & interTAXA$connect=="np",])




##get list of bacterial OTUs that were associated with fungi
intbactOTUs <-interTAXA[interTAXA$bact_fung == "1",]
##subset bacterial taxonomy for those OTUs
intbactOTUtax <- tax_16s[tax_16s$otu %in% intbactOTUs$end,]
#inspect the above table to determine common bacterial taxa

##get list of fungal OTUs that were associated with bacteria
intfungiOTUs <-interTAXA[interTAXA$bact_fung == "1",]
##subset bacterial taxonomy for those OTUs
intfungiOTUtax <- keep_otus[keep_otus$otu %in% intfungiOTUs$start,]
#inspect the above table to determine common fungal taxa
