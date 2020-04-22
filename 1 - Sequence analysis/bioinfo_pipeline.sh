# Analyis GLBRC soil/phyllosphere switcgrass fungal communities
# **************************************************** -----------------------

# recommended softwares:
USEARCH https://www.drive5.com/usearch/
VSEARCH https://github.com/torognes/vsearch
cutadapt https://cutadapt.readthedocs.io/en/stable/
seqtk https://github.com/lh3/seqtk
bowtie2 http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
seaview http://doua.prabi.fr/software/seaview
CONSTAX https://github.com/Gian77/CONSTAX/

#1) check quality of original raw sequence data---------------------------------
mkdir stats
cd raw_reads/

for i in *.fastq
do echo "$i : `grep -c "^+$" $i`"
done > ../stats/reads_raw.counts

fastqc ./*R1_001.fastq
fastqc ./*R2_001.fastq 

#2) mergin PE reads ------------------------------------------------------------
mkdir ../merged

./usearch10 -fastq_mergepairs *.fastq -relabel @ -fastqout merged.fastq -fastq_pctid 85 -interleaved -threads 20

mv merged.fastq ../merged/

grep -c "^+$" ../merged/merged.fastq  > ../stats/merged.counts

# 3) Removing Phix genome form reads -------------------------------------------
mkdir ../no_phix
cd ../merged/

./bowtie2 -t -S --un no_Phix.fastq /<path-to-Phix-index>/phix_index/ merged.fastq Phix_alignments.sam 2> no_Phix.log

grep -c "^+$" no_Phix.fastq  > ../stats/no_Phix.counts

# 4) Removing primers form reads
mkdir ../stripped
cd ../no_phix

#Fungal ITS2 primers used in the study -----------------------------------------
#ITS9F GAACGCAGCRAAIIGYGA
#ITS4R TCCTCCGCTTATTGATATGC
#ITS4R_rc (reverse complement):GCATATCAATAAGCGGAGGA

./cutadapt -g GAACGCAGCRAANNGYGA -a GCATATCAATAAGCGGAGGA -f fastq -n 2 --discard-untrimmed --match-read-wildcards -o stripped.fastq no_Phix.fastq > stripped.txt

mv stripped.fastq ../stripped/
mv stripped.txt ../stripped/

cd ../stripped

# subsample reads for checking and converting .fastq to .fasta
./seqtk sample -s100 stripped.fastq 500 > sub_stripped.fastq
./seqtk seq -aQ64 sub_stripped.fastq > sub_stripped.fasta

grep -c "^+$" stripped.fastq > ../stats/stripped.counts

# 5) trimming out conserved regions and generating statistics -----------------
# this step is based on the alignmne of the 500 randomply subsetted reads
# use seaview to visualize and align
mkdir ../filtered
cd ../filtered

./usearch10 -fastq_filter stripped.fastq -fastq_stripleft 95 -fastqout trimmed_left.fastq
./seqtk trimfq -e 16 trimmed_left.fastq > trimmed_both.fastq

# getting starts
./vsearch -fastq_stats trimmed_both.fastq -log ../stats/stats_trimmed.txt
./usearch10 -fastq_eestats2 trimmed_both.fastq -output ../stats/eestats2_trimmed.txt -length_cutoffs 100,500,1

# 6) Filtering out reads with errors -------------------------------------------
./usearch10 -fastq_filter trimmed_both.fastq -fastq_maxee 1.0 -fastq_trunclen 180 -fastq_maxns 0 -fastqout filtered_trimmed.fastq -fastaout filtered_trimmed.fasta

# check trimming reults on a subsample of the reads
seqtk sample -s100 filtered_trimmed.fasta 500 > sub_filtered_trimmed.fasta 

# 7) Clustering OTUs -----------------------------------------------------------
mkdir ../clustered
cd ../clustered

./usearch -fastx_uniques ../filtered/filtered_trimmed.fastq -fastqout uniques.fastq -sizeout
./usearch -cluster_otus uniques.fastq -minsize 2 -otus otus.fasta -uparseout otus.txt -relabel OTU_
./usearch -usearch_global ../filtered_trimmed.fasta -db otus.fasta -strand plus -id 0.97 -otutabout otu_table.txt

# 8) assigning taxonomy to OTU representative sequences ------------------------

For this step please refer to the official page of CONSTAX
https://github.com/Gian77/CONSTAX/

# ***************************** please cite ************************************
# Gdanetz K, Benucci GMN, Vande Pol N, Bonito G (2017) CONSTAX: a tool for improved 
# taxonomic resolution of environmental fungal ITS sequences. BMC Bioinformatics 
# 18:538 doi 10.1186/s12859-017-1952-x


