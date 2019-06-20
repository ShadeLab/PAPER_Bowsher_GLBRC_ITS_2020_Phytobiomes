# ******* DATA ANALYSIS ********************************************** -----
# Project name: swihcgrass leaf and soil microbiome - GLBRC
# Manuscript:   
# Authors:      
# Affiliation:  Michigan State University
# Journal:      Phytobiome special issue (Microbiome meeting at Berkeley) 
# Date:         April 29, 2019
# ******************************************************************** -----

# WORKING ENVIRONMENT SETUP ------------------------------------------------
options(scipen = 9999) #to use decimals
options(max.print=100000000) # to print more lines on the display

# loading required packages ------------------------------------------------
library(phyloseq); packageVersion("phyloseq")
library(Biostrings)
library(ggplot2)
library(indicspecies)
library(vegan)
library(ape)


# phyloseq object TRIMMED reads --------------------------------------------

otus_ITS_trim <- read.delim("Nico_ITS_trimmed_otus/otu_table_ITS_UPARSE_filtered_180.txt",row.names=1) 
otus_phy_ITS_trim <-otu_table(otus_ITS_trim, taxa_are_rows = TRUE)

metadata_ITS_trim <-read.csv("Nico_ITS_trimmed_otus/metadata.csv", row.names=1, header=TRUE, sep=",")
metadata_phy_ITS_trim <-sample_data(metadata_ITS_trim)

taxonomy_ITS_trim <-read.delim("Nico_ITS_trimmed_otus/constax_taxonomy/consensus_taxonomy.txt", header=TRUE, row.names=1)
taxonomy_ITS_trim <- subset(taxonomy_ITS_trim, select=c(Kingdom, Phylum, Class, Order, Family, Genus, Species))
taxonomy_phy_ITS_trim <- tax_table(as.matrix(taxonomy_ITS_trim))

otus_seq_ITS_trim <- readDNAStringSet("Nico_ITS_trimmed_otus/filtered_otus.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)

physeq_obj_ITS_trim <- phyloseq(otus_phy_ITS_trim, 
                                metadata_phy_ITS_trim,
                                taxonomy_phy_ITS_trim,
                                otus_seq_ITS_trim) 

physeq_obj_ITS_trim
str(physeq_obj_ITS_trim)

# correcting taxonomies --------------------------------------------
tax_table(physeq_obj_ITS_trim)[tax_table(physeq_obj_ITS_trim)==""]<- NA
tax_table(physeq_obj_ITS_trim)[is.na(tax_table(physeq_obj_ITS_trim))]<-"Unclassified"

# checking for non-target taxa ------------------ 
unique(as.data.frame(tax_table(physeq_obj_ITS_trim))$Kingdom) -> uniques_Kingdom
as.character(uniques_Kingdom)

unique(as.data.frame(tax_table(physeq_obj_ITS_trim))$Phylum) -> uniques_Phylum
as.character(uniques_Phylum)

physeq_obj_Cercozoa <- subset_taxa(physeq_obj_ITS_trim, Phylum == "Cercozoa")
tax_table(physeq_obj_Cercozoa)

physeq_obj_Chlorophyta <- subset_taxa(physeq_obj_ITS_trim, Phylum == "Chlorophyta")
tax_table(physeq_obj_Chlorophyta)


# if you want to re-check OTU taxonomy -------------------
physeq_obj_Unclassified <- subset_taxa(physeq_obj_ITS_trim, Phylum == c("Unclassified"))
tax_table(physeq_obj_Unclassified)

write.dna(refseq(physeq_obj_Unclassified), format="fasta", 
          colsep="", file="physeq_obj_unclassified.fasta")

# >>> Filtering out OTUs ----------------------------------
# Oliver et al. 2015, PCR errors and tag switching
# Lindhal et al. 2013, tag switching - that's a good  one!
# Barberan et al. 2012, removing OTUs that appear in less than x samples

physeq_obj_fungi_trim <- subset_taxa(physeq_obj_ITS_trim, Phylum != "Cercozoa")
physeq_obj_fungi_trim <- subset_taxa(physeq_obj_fungi_trim, Phylum != "Chlorophyta")
physeq_obj_fungi_trim


physeq_obj_fungi_trim -> physeq_obj_fungi_trim_filt
otu_table(physeq_obj_fungi_trim_filt) <- otu_table(physeq_obj_fungi_trim_filt)[which(rowSums(otu_table(physeq_obj_fungi_trim_filt)) >= 10),] ### PCR Errors 
otu_table(physeq_obj_fungi_trim_filt)[otu_table(physeq_obj_fungi_trim_filt) <= 4] <- 0 ### tag switching
otu_table(physeq_obj_fungi_trim_filt) <- otu_table(physeq_obj_fungi_trim_filt)[which(rowSums(otu_table(physeq_obj_fungi_trim_filt)) > 0),] 
physeq_obj_fungi_trim_filt


# rarefaction curves ---------------------------------------
library("iNEXT")
library("vegan")

otu_fungi_soil_trim <- as.data.frame(otu_table(physeq_obj_fungi_trim_filt))

rarecurve_soil_fungi <- rarecurve(t(otu_fungi_soil_trim), main = "soil 180bp OTUs", step = 1000, label = FALSE)
rarecurve_soil_fungi

out_fungi_trim <- iNEXT(otu_fungi_soil_trim, q=0, datatype="abundance")
out_fungi_trim

curves_trim <- ggiNEXT(out_fungi_trim, type = 1, se = FALSE) +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) + 
  theme(axis.title = element_text(size = 14)) + 
  labs(title="trimmed 180bp OTUs", x="Number or DNA reads", y="Number of OTUs") +
  theme(legend.position="bottom") +
  #scale_x_continuous(limits=c(0, 1000000),) +
  scale_x_continuous(breaks=c(10000, 50000, 80000,
                              200000, 500000, 1000000)) +
  theme(axis.text.x = element_text(angle=90, size=9, vjust = .5)) +
  guides(colour=FALSE, shape=FALSE) 

curves_trim$layers
curves_trim$layers <- curves_trim$layers[-1]
curves_trim

# rarefying data at min library size--------------------------------------------
physeq_obj_fungi_trim
min(sample_sums(physeq_obj_fungi_trim))
physeq_obj_fungi_trim_ev = rarefy_even_depth(physeq_obj_fungi_trim,
                                             sample.size =86710,rngseed = FALSE, 
                                             replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
physeq_obj_fungi_trim_ev
colSums(otu_table(physeq_obj_fungi_trim_ev))
any(taxa_sums(physeq_obj_fungi_trim_ev) == 0)


# alpha-diversity -----------------------------------

sample_data(physeq_obj_fungi_trim_ev)$Date <- as.factor(sample_data(physeq_obj_fungi_trim_ev)$Date)
sample_data(physeq_obj_fungi_trim_ev)$Date <- factor(sample_data(physeq_obj_fungi_trim_ev)$Date, 
                                                     levels=c("4/24/2017","5/15/2017","6/5/2017","6/26/2017","7/17/2017","8/7/2017","8/28/2017","9/18/2017"))

alpha_ITS_trim = plot_richness(physeq_obj_fungi_trim_ev, x="Date", color="Type", measures = c("Observed")) +
  labs(title="trimmed 180bp OTUs", x="Date", y = "Rarefied Richness") +
  geom_point(size=0.1, alpha = 0.9) + 
  geom_point(position = position_dodge(width = 0.75)) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black", outlier.shape = 8) +
  facet_grid(~Treatment~Type, scales = "free_x") + 
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) + 
  #expand_limits(x = 0, y = 0) +
  ylim(0, 2000) +
  theme(legend.position="bottom")

alpha_ITS_trim
alpha_ITS_trim$layers
alpha_ITS_trim$layers <- alpha_ITS_trim$layers[-1]
alpha_ITS_trim$layers <- alpha_ITS_trim$layers[-1]
alpha_ITS_trim$layers <- alpha_ITS_trim$layers[-1]


# beta diversity - PCoA ------------------------------------------

pcoa_physeq_obj_fungi_trim_ev <- ordinate(physeq_obj_fungi_trim_ev, method ="PCoA", distance="bray")
pcoa_physeq_obj_fungi_trim_ev


plot_pcoa_trim = plot_ordination(physeq_obj_fungi_trim_ev, pcoa_physeq_obj_fungi_trim_ev, 
                                 color="Date", shape="Treatment") + 
  geom_point(size=2.5, alpha=0.9, aes(shape=Treatment)) +
  labs(title="PCoA - trimmed 180bp OTUs") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  theme(legend.position="right")

plot_pcoa_trim

