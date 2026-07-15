```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("microbiome")
install.packages("vegan")
install.packages("devtools")
BiocManager::install("microbiome")

install.packages(pkgs=c("BiodiversityR", "vegan","Rcmdr", "MASS", "mgcv","cluster", "RODBC", "rpart", "effects", "multcomp","ellipse", "maptree", "sp", "splancs", "spatial","akima", "nnet", "dismo", "raster", "rgdal", "bootstrap", "PresenceAbsence","maxlike", "gbm", "randomForest", "gam", "earth", "mda","kernlab", "e1071", "glmnet", "sem", "rgl", "relimp","lmtest", "leaps", "Hmisc", "colorspace", "aplpack","abind", "XLConnect", "car", "markdown", "knitr","geosphere", "maptools", "rgeos", "ENMeval", "red"),dependencies=c("Depends", "Imports"))
```

# Load Pacakges
```{r}
#All the packages you will need for entire code, some of these could be unnecessary pending what you want to graph
library(dada2)
library(BiocManager)
library(ggplot2); packageVersion("ggplot2")
library(devtools)
library(vegan)
library(dbplyr)
library(microbiome)
library(tidyverse)
library(conflicted)
library(dplyr)
library(cowplot)
library(ggpubr)


formatPvalues <- function(pvalue) {
  ra<- ""
  if(pvalue <= 0.1) ra<- "."
  if(pvalue <= 0.05) ra<- "*"
  if(pvalue <= 0.01) ra<- "**"
  if(pvalue <= 0.001) ra<- "***"
  return(ra)
}
```

```{r}
# Path needs to have fastq files unzipped!
path <- "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/fastq"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_S"), `[`, 1)
list(sample.names)

#plotQualityProfile(fnFs[1:2])
#plotQualityProfile(fnRs[1:2])

FWD <- "GTGYCAGCMGCCGCGGTAA"
REV <- "CCGYCAATTYMTTTRAGTTT"
trimLeft = c(FWD, REV)
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(280,200), maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE, trimLeft = c(19,20))
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#plotErrors(errF, nominalQ=TRUE)
#plotErrors(errR, nominalQ=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
#dadaFs[[1]]

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
#dadaRs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
#head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#Identified 26085 bimeras out of 26370 input sequences.
#[1]  37 285
#[1] 0.8278426

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/track_sequences.csv", sep=",", quote=F, col.names=NA)
```

```{r}
taxa <- assignTaxonomy (seqtab.nochim,"/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/fastq/silva_nr99_v138.2_toGenus_trainset.gz")
```

```{r}
taxa.print <- taxa
rownames(taxa.print) <- NULL
#head(taxa.print)

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
asv_headers[i] <- paste(">ASV", i, sep="_")
}
```

```{r}
asv_fasta <- c(rbind(asv_headers, asv_seqs))

asv_otu <- t(seqtab.nochim)
row.names(asv_otu) <- sub(">", "", asv_headers)

asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)

otu_tax_table <- merge(asv_otu, asv_tax, by=0)
```

```{r}
write.table(asv_fasta, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/asv_fasta_enrich.fa")
write.table(asv_otu, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/asv_otu_enrich.csv", sep=",", quote=F, col.names=NA)
write.table(asv_tax, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/asv_tax_enrich.csv", sep=",", quote=F, col.names=NA)
write.table(otu_tax_table, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/asv_otu_tax_enrich.csv", sep=",", quote=F, col.names=NA)
```

```{r}
metadata <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/enrichment_metadata.csv")

otu_counts <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/asv_otu_enrich.csv") %>%
  pivot_longer(-ASV, names_to="sample_id", values_to = "count")
otu_counts

taxonomy <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/asv_tax_enrich.csv")
taxonomy
```

```{r}
library(phyloseq)
library(Biostrings)
```

```{r}
otu_rel_abund <- inner_join(metadata, otu_counts, by="sample_id") %>%
  inner_join(., taxonomy, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Order", "Class", "Genus", "ASV"),
         names_to="level",
         values_to="taxon")
otu_rel_abund

write.table(otu_rel_abund, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/otu_rel_abund_enrich.csv", sep=",", quote=F, col.names=NA)

otu_rel_abund <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/otu_rel_abund_enrich.csv")
otu_rel_abund
```

# Relative Abundance Based on Groups
```{r}
#change LEVEL to what ever taxon level you want, group_by factors you want to compare, and taxon == to taxon of interest
#change GROUP_BY(XX) %>% to just look at taxon as a whole
otu_rel_abund %>%
  filter(level=="Family") %>%
  group_by(microbe_type, sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(microbe_type, taxon) %>%
    summarize( mean_rel_abund = 100 * mean(rel_abund), se_rel_abund   = 100 * sd(rel_abund) / sqrt(n()), n = n(), .groups = "drop") %>%
  filter(taxon == "Thiomicrospira")
```

# Enrich Barcharts: FeOB + SRB
```{r}
## Phylum
otu_rel_abund %>%
  filter(level=="Phylum") %>%
  group_by(sample_id, metal_type, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(metal_type, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=metal_type, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=metal_type, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, y="Mean Relative Abundance (%)") + theme_classic()
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/phylum_stacked_barchart_metal.tiff", width=10, height=10)

## Class
otu_rel_abund %>%
  filter(level=="Class") %>%
  group_by(sample_id, metal_type, microbe_type, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(microbe_type, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=microbe_type, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=microbe_type, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, y="Mean Relative Abundance (%)") + theme_classic()
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/class_stacked_barchart_enrich_metal.tiff", width=12, height=10)

## Order
otu_rel_abund %>%
  filter(level=="Order") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/order_stacked_barchart_enrich.tiff", width=15, height=10, limitsize = FALSE)

## Family
otu_rel_abund %>%
  filter(level=="Family") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/family_stacked_barchart_enrich.tiff", width=15, height=10, limitsize = FALSE)

## Genus
otu_rel_abund %>%
  filter(level=="Genus") %>%
  group_by(sample_id, metal_type, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(metal_type, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=metal_type, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=metal_type, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/genus_stacked_barchart_enrich.tiff", width=13, height=10, limitsize = FALSE)
```

# Phyloseq: Enrich
```{r}
library(phyloseq)
library(Biostrings)

# Read data into R
otu_tab <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/asv_otu_enrich.csv")
otu_tab

tax_tab <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/asv_tax_enrich.csv")
tax_tab

samples_df <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/enrichment_metadata.csv")
samples_df

# Phyloseq objects need to have row.names
otu_tab <- otu_tab %>%
  tibble::column_to_rownames("ASV")

tax_tab <- tax_tab %>%
  tibble::column_to_rownames("ASV")

samples_df <- samples_df %>%
  tibble::column_to_rownames("sample_id")

# Transform OTU & Tax table into matrices
otu_tab <- as.matrix(otu_tab)
tax_tab <- as.matrix(tax_tab)

#Transform into Phyloseq Objects
ASV = otu_table(otu_tab, taxa_are_rows = TRUE)
TAX = tax_table(tax_tab)
sample_id = sample_data(samples_df)
  
ps <- phyloseq(ASV, TAX, sample_id)
ps

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 285 taxa and 37 samples ]
#sample_data() Sample Data:       [ 37 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 285 taxa by 6 taxonomic ranks ]
```

```{r}
# Visualize Data
sample_names(ps)
rank_names(ps)
sample_variables(ps)

# Normalize number of reads in each sample using median sequencing depth
total = median(sample_sums(ps))
standf = function(x, t=total) round(t * (x / sum(x)))
ps = transform_sample_counts(ps, standf)

# Boxplot based on division
plot_taxa_prevalence(ps, fill = "Family") +
  geom_boxplot(aes(color=Family), stat="identity", position="stack")

# Bar graphs based on division
plot_bar(ps, fill = "Family") + 
  geom_bar(aes(fill=Family), stat="identity", position="stack") +
  facet_grid(~metal_type)
```

# Boxplots Raxa Rank
```{r}
require(dplyr) #required for data table filtering
require(ggplot2) #require for plotting
require(scales) #required for percentage conversion
require(ggthemes) #required for colorblind palette
require(ggpubr) #required for wilcoxon rank sum test statistical comparison between explanatory variables (categorical)

ps_family <- ps %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance (or use ps.ra)
  psmelt() %>%                                         # Melt to long format for easy ggploting
  filter(Abundance > 0.01)                             # Filter out low abundance taxa

write.table(ps_family, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/phyloseq_family_abund.csv", sep=",", quote=F, col.names=NA)
ps_family

box_plot <- ggplot(ps_family, aes(x=microbe_type, y=Abundance, fill=Family)) +
	geom_boxplot() +
  facet_grid(~metal_type + location) +
	scale_y_continuous(labels=percent) +
	xlab("Location") + 
	ylab("Average Relative Abundance\n Out of Total Sequences") +
	theme(axis.text.x=element_text(angle=60, hjust=1)) + #this will tilt the phyla text so that it looks cleaner +
  scale_colour_colorblind() +
	scale_colour_colorblind(name="Family") #+ #this will make the geom_points the same colors as the geom_boxplot
	#geom_point(position=position_dodge(width=0.75), aes(color=Family),show.legend=FALSE)
box_plot #prints box plot to standard out

#Include this to export this figure as a .tiff publication quality file
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/Taxa_Rank_Enrich_Family.tiff", plot=box_plot, device=NULL, path=NULL, scale=1, width=15, height=7, dpi=900, limitsize=TRUE)
```

# OTU Table Creating
```{r}
df_otu <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/df_asv_otu_enrich.csv")
df_otu
```

```{r}
otu_count_all <- df_otu %>%
  pivot_longer(-sample_id, names_to = "ASV", values_to = "count")
otu_count_all

write.table(otu_count_all, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/otu_count_all_enrich.csv", sep=",", quote=F, col.names=NA)

otu_count_all <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/otu_count_all_enrich.csv")
otu_count_all

otu_count_all %>%
  group_by(sample_id) %>%
  mutate(total = sum(count)) %>%
  filter(total > 5000) %>%
  group_by(ASV) %>%
  mutate(total=sum(count)) %>% 
  filter(total != 0) %>%
  as.data.frame()
#Going to set threshold at 5000
```

# NMDS
```{r}
nmds_asv_otu_all <- inner_join(metadata, df_otu, by="sample_id")
nmds_asv_otu_all

write.table(nmds_asv_otu_all, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/nmds_asv_otu_enrich.csv", sep=",", quote=F, col.names=NA)

pc <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/nmds_asv_otu_enrich.csv")
pc

#make community matrix: extract columns with ASV information
com <- pc[,6:ncol(pc)]
com

#turn ASV information into a matrix
m_com <- as.matrix(com)

#Run NMDS using Bray-Curtis distance
set.seed(123)
nmds <- metaMDS(m_com, distance="bray") #stress = nearly zero 
nmds
plot(nmds)

#access the specific points data of the NMDS plot & scores
str(nmds)
nmds$points
scores(nmds)

#extract NMDS scores
data.scores = as.data.frame(scores(nmds)$sites)

#add columns to data frame
data.scores$sample_id = pc$sample_id
data.scores$location = pc$location
data.scores$microbe_type = pc$microbe_type
data.scores$metal_type = pc$metal_type

head(data.scores)
```

```{r}
# By Location
xx = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = location, shape = metal_type))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "Wreck Location", y = "NMDS2")
xx
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/NMDS_enrich_metal_location.tiff", width = 30, height = 10)
```

```{r}
# By Metal Type
xx = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = metal_type))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "Metal", y = "NMDS2")
xx
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/NMDS_enrich_metal_type.tiff", width = 10, height = 10)
```

```{r}
# By Microbe Type
xx = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = microbe_type, shape = location))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "Enrichment Classification", y = "NMDS2")
xx
#ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/NMDS_enrich_microbe.tiff", width = 15, height = 10)
```

# By Metal Type & Location
```{r}
xx = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(shape = metal_type, colour = location))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "Location", y = "NMDS2")
xx
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/NMDS_all_samples_metal_location_type.tiff", width = 10, height = 10)
```

# Shannon, Richness, Simpson and Evenness: All Enrich
```{r}
richness <- function(x){
  sum(x > 0)}

shannon <- function(x){
  relabund <- x[x>0]/sum(x)
  -sum(relabund * log(relabund))
}

simpson <- function(x){
  n <- sum(x)
  sum(x * (x-1) /(n*n-1))
}
```

```{r}
otu_count_all <- otu_count_all %>%
  group_by(sample_id) %>%
  summarize(richness = richness(count),
            shannon = shannon(count), 
            evenness = shannon/log(richness),
            n=sum(count))

write.table(otu_count_all, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/otu_count_all_enrich_diversity_metrics.csv", sep=",", quote=F, col.names=NA)

diversity_metrics <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/otu_count_all_enrich_diversity_metrics.csv")
diversity_metrics

write.table(diversity_metrics, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/enrich_diversity_metrics.csv", sep=",", quote=F, col.names=NA)
```

# Boxplots: All Enrich
```{r}
diversity_metrics <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/enrich_diversity_metrics.csv")
diversity_metrics

diversity_metrics %>%
  group_by(location) %>%
  pivot_longer(cols=c(richness, shannon, evenness), 
               names_to="metric") %>%
ggplot(aes(x=n, y=value, fill= microbe_type)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 8, outlier.size = 2) +
  facet_wrap(~metric, nrow=1, scales="free_y")

ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/alpha_diversity_metrics_enrich_microbe_location.tiff", width = 10, height = 7)
```

```{r}
diversity_metrics %>%
  group_by(location) %>%
  pivot_longer(cols=c(richness, shannon, evenness), 
               names_to="metric") %>%
ggplot(aes(x=n, y=value, fill= location)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 8, outlier.size = 2) +
  facet_wrap(~metric, nrow=1, scales="free_y")

ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/alpha_diversity_metrics_enrich_sample_location.tiff", width = 10, height = 7)
```

```{r}
diversity_metrics %>%
  group_by(location) %>%
  pivot_longer(cols=c(richness, shannon, evenness), 
               names_to="metric") %>%
ggplot(aes(x=n, y=value, fill= metal_type)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 8, outlier.size = 2) +
  facet_wrap(~metric, nrow=1, scales="free_y")

ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/alpha_diversity_metrics_enrich_metal_type.tiff", width = 10, height = 5)

#Each point represents a sample, (n) Sum of Count, (X) Total number of sequences for each sample & (Y) Value of diversity metric
```

# ANOSIM Statistical Testing:
```{r}
#ANOSIM give you the P value (i.e. significance levels) and a R value (i.e. the strength of the factors on the samples)
# R value is supposed to vary between 0 and 1 (not between -1 and +1) but you can obtained negative values but they are always close to 0. 
# R value close to 1 indicates high separation between levels of your factor (e.g. control vs treatment samples), while R value close to 0 indicate no separation between levels of your factor.
```

```{r}
pc_ano <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/nmds_asv_otu_enrich.csv")
pc_ano
```

```{r}
# By Location
com_ano = pc_ano[,6:ncol(pc_ano)]
m_com_ano = as.matrix(com_ano)
ano_all = anosim(m_com_ano, pc_ano$location, distance = "bray", permutations = 9999)
ano_all

# ANOSIM statistic R: 0.02383
      # Significance: 0.2663
```

```{r}
# By Metal Type
com_ano = pc_ano[,6:ncol(pc_ano)]
m_com_ano = as.matrix(com_ano)
ano_all = anosim(m_com_ano, pc_ano$metal_type, distance = "bray", permutations = 9999)
ano_all

# ANOSIM statistic R: -0.02724 (No group separation)
      # Significance: 0.5849
```

```{r}
# By Microbe Classification
com_ano = pc_ano[,6:ncol(pc_ano)]
m_com_ano = as.matrix(com_ano)
ano_all = anosim(m_com_ano, pc_ano$microbe_type, distance = "bray", permutations = 9999)
ano_all

# ANOSIM statistic R: 0.9944
      # Significance: 1e-04
```

# SIMPER Analysis
```{r}
library(vegan)

mydata <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/simper_asv_otu.csv")
mydata

library(vegan)

# Grouping variable
group <- mydata$location
group2 <- mydata$metal_type
group3 <- mydata$microbe_type

# Community matrix (remove SampleID and Location)
community <- mydata[, !(names(mydata) %in% c("sample_id", "location", "metal_type"))]

# Check that all remaining columns are numeric
sapply(community, class)
```

```{r}
# Run SIMPER: Location
simper_result <- simper(community, group)
summary(simper_result)

#install.packages("openxlsx", dependencies = TRUE)
library(openxlsx)

# View the comparison names
names(simper_result)

# Extract the first comparison
simper_df <- as.data.frame(simper_result[[1]])

# Write to Excel
write.xlsx(simper_df, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/SIMPER_results_location.xlsx", rowNames = TRUE)
```

```{r}
# Run SIMPER: Metal Type
simper_result2 <- simper(community, group2)
summary(simper_result2)

#install.packages("openxlsx", dependencies = TRUE)
library(openxlsx)

# View the comparison names
names(simper_result2)

# Extract the first comparison
simper_df2 <- as.data.frame(simper_result2[[1]])

# Write to Excel
write.xlsx(simper_df2, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/SIMPER_results_metal.xlsx", rowNames = TRUE)
```

```{r}
# Run SIMPER: Microbe Type
simper_result3 <- simper(community, group3)
summary(simper_result3)

#install.packages("openxlsx", dependencies = TRUE)
library(openxlsx)

# View the comparison names
names(simper_result3)

# Extract the first comparison
simper_df3 <- as.data.frame(simper_result3[[1]])

# Write to Excel
write.xlsx(simper_df3, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/SIMPER_results_microbe.xlsx", rowNames = TRUE)
```

# FeOB Enrichment Only
```{r}
metadata_F <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/metadata_FeOB.csv")
metadata_F

otu_counts_F <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/asv_otu_FeOB.csv") %>%
  pivot_longer(-ASV, names_to="sample_id", values_to = "count")
otu_counts_F

taxonomy_F <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/asv_tax_FeOB.csv")
taxonomy_F

otu_rel_abund_F <- inner_join(metadata_F, otu_counts_F, by="sample_id") %>%
  inner_join(., taxonomy_F, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Family", "ASV"),
         names_to="level",
         values_to="taxon")
otu_rel_abund_F

write.table(otu_rel_abund_F, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/otu_rel_abund_FeOB.csv", sep=",", quote=F, col.names=NA)

otu_rel_abund_F <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/otu_rel_abund_FeOB.csv")
otu_rel_abund_F
```

# Feob Barcharts
```{r}
## Phylum
otu_rel_abund_F %>%
  filter(level=="Phylum") %>%
  group_by(sample_id, metal_type, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(metal_type, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=metal_type, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=metal_type, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, y="Mean Relative Abundance (%)") + theme_classic()
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/FeOB/phylum_stacked_barchart_FeOB_metal.tiff", width=15, height=10)

## Class
otu_rel_abund_F %>%
  filter(level=="Class") %>%
  group_by(sample_id, metal_type, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(metal_type, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=metal_type, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=metal_type, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/FeOB/class_stacked_barchart_FeOB_metal.tiff", width=15, height=10)

## Order
otu_rel_abund_F %>%
  filter(level=="Order") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/FeOB/order_stacked_barchart_FeOB.tiff", width=15, height=10, limitsize = FALSE)
```

```{r}
## Family
otu_rel_abund_F %>%
  filter(level=="Family") %>%
  group_by(sample_id, metal_type, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(metal_type, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=metal_type, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=metal_type, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/FeOB/family_stacked_barchart_FeOB_metal_type.tiff", width=15, height=10, limitsize = FALSE)
```

```{r}
## Genus
otu_rel_abund_F %>%
  filter(level=="Genus") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/FeOB/genus_stacked_barchart_FeOB.tiff", width=15, height=10, limitsize = FALSE)
```

# Phyloseq: FeOB ONLY
```{r}
library(phyloseq)
library(Biostrings)

# Read data into R
otu_tabF <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/asv_otu_FeOB.csv")
otu_tabF

tax_tabF <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/asv_tax_FeOB.csv")
tax_tabF

samples_dfF <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/metadata_FeOB.csv")
samples_dfF

# Phyloseq objects need to have row.names
otu_tabF <- otu_tabF %>%
  tibble::column_to_rownames("ASV")

tax_tabF <- tax_tabF %>%
  tibble::column_to_rownames("ASV")

samples_dfF <- samples_dfF %>%
  tibble::column_to_rownames("sample_id")

# Transform OTU & Tax table into matrices
otu_tabF <- as.matrix(otu_tabF)
tax_tabF <- as.matrix(tax_tabF)

#Transform into Phyloseq Objects
ASV = otu_table(otu_tabF, taxa_are_rows = TRUE)
TAX = tax_table(tax_tabF)
sample_id = sample_data(samples_dfF)
  
psF <- phyloseq(ASV, TAX, sample_id)
psF

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 285 taxa and 4 samples ]
#sample_data() Sample Data:       [ 4 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 285 taxa by 6 taxonomic ranks ]
```

```{r}
# Visualize Data
sample_names(psF)
rank_names(psF)
sample_variables(psF)

# Normalize number of reads in each sample using median sequencing depth
total = median(sample_sums(psF))
standf = function(x, t=total) round(t * (x / sum(x)))
psF = transform_sample_counts(psF, standf)

# Bar graphs based on division
plot_bar(psF, fill = "Family") + 
  geom_bar(aes(fill=Family), stat="identity", position="stack") +
  facet_grid(~metal_type + location)
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/FeOB/metal_location_bar.tiff", width=15, height=10, limitsize = FALSE)
```

# Boxplots Raxa Rank
```{r}
require(dplyr) #required for data table filtering
require(ggplot2) #require for plotting
require(scales) #required for percentage conversion
require(ggthemes) #required for colorblind palette
require(ggpubr) #required for wilcoxon rank sum test statistical comparison between explanatory variables (categorical)

psF_family <- psF %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance (or use ps.ra)
  psmelt() %>%                                         # Melt to long format for easy ggploting
  filter(Abundance > 0.01)                             # Filter out low abundance taxa

write.table(psF_family, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/phyloseq_family_abund_FEOB.csv", sep=",", quote=F, col.names=NA)
psF_family

box_plot <- ggplot(psF_family, aes(x=microbe_type, y=Abundance, fill=Family)) +
	geom_boxplot() +
  facet_grid(~metal_type + location) +
	scale_y_continuous(labels=percent) +
	xlab("Location") + 
	ylab("Average Relative Abundance\n Out of Total Sequences") +
	theme(axis.text.x=element_text(angle=60, hjust=1)) + #this will tilt the phyla text so that it looks cleaner +
  scale_colour_colorblind() +
	scale_colour_colorblind(name="Family") #+ #this will make the geom_points the same colors as the geom_boxplot
	#geom_point(position=position_dodge(width=0.75), aes(color=Family),show.legend=FALSE)
box_plot #prints box plot to standard out

#Include this to export this figure as a .tiff publication quality file
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/Taxa_Rank_FEOB_Family.tiff", plot=box_plot, device=NULL, path=NULL, scale=1, width=15, height=7, dpi=900, limitsize=TRUE)
```

# FeOB Relative Abundance Based on Groups
```{r}
#change LEVEL to what ever taxon level you want, group_by factors you want to compare, and taxon == to taxon of interest
otu_rel_abund_F %>%
  filter(level=="Family") %>%
  group_by(metal_type, sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(metal_type,taxon) %>%
    summarize(
    mean_rel_abund = 100 * mean(rel_abund),
    se_rel_abund   = 100 * sd(rel_abund) / sqrt(n()),
    n              = n(),
    .groups = "drop"
  ) %>%
  filter(taxon == "Desulfovibrionaceae")
```

```{r}
df_otu_F <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/df_otu_FeOB.csv")
df_otu_F

otu_count_all_F <- df_otu_F %>%
  pivot_longer(-sample_id, names_to = "ASV", values_to = "count")
otu_count_all_F

write.table(otu_count_all_F, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/otu_count_all_FeOB.csv", sep=",", quote=F, col.names=NA)

otu_count_F <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/otu_count_all_FeOB.csv")
otu_count_F

otu_count_all_F %>%
  group_by(sample_id) %>%
  mutate(total = sum(count)) %>%
  filter(total > 5000) %>%
  group_by(ASV) %>%
  mutate(total=sum(count)) %>% 
  filter(total != 0) %>%
  as.data.frame()
#Going to set threshold at 5000

nmds_asv_otu_F <- inner_join(metadata_F, df_otu_F, by="sample_id")
nmds_asv_otu_F

write.table(nmds_asv_otu_FeOB, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/nmds_asv_otu_FeOB.csv", sep=",", quote=F, col.names=NA)

pc <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/nmds_asv_otu_FeOB.csv")
pc

#make community matrix: extract columns with ASV information
com <- pc[,6:ncol(pc)]
com

#turn ASV information into a matrix
m_com <- as.matrix(com)

#Run NMDS using Bray-Curtis distance
set.seed(123)
nmdsF <- metaMDS(m_com, distance="bray") #stress = nearly zero 
nmdsF
plot(nmdsF)

#access the specific points data of the NMDS plot & scores
str(nmdsF)
nmds$points
scores(nmdsF)

#extract NMDS scores
data.scores = as.data.frame(scores(nmdsF))

#add columns to data frame
data.scores$sample_id = pc$sample_id
data.scores$location = pc$location
data.scores$microbe_type = pc$microbe_type
data.scores$metal_type = pc$metal_type

head(data.scores)
```

## NMDS: FeOB
```{r}
# By Location
xxF = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = location))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "Wreck Location", y = "NMDS2")
xxF
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/FeOB/NMDS_FeOB_samples_location.tiff", width = 10, height = 10)

# By Metal Type
xxF = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = metal_type, shape = location))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "Metal", y = "NMDS2")
xxF
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/FeOB/NMDS_FeOB_metal_location.tiff", width = 10, height = 10)

# By Microbe Type
xx = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = microbe_type))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "Microbe Classification", y = "NMDS2")
xx
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/FeOB/NMDS_FeOB_samples_microbe.tiff", width = 10, height = 10)

# By Metal Type & Location
xx = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(shape = metal_type, colour = location))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "Location", y = "NMDS2")
xx
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/FeOB/NMDS_FeOB_samples_metal_location_type.tiff", width = 10, height = 10)
```

# Shannon, Richness, Eveness: FeOB
```{r}
otu_count_F <- otu_count_F %>%
  group_by(sample_id) %>%
  summarize(richness = richness(count),
            shannon = shannon(count), 
            evenness = shannon/log(richness),
            n=sum(count))

write.table(otu_count_F, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/otu_count_FEOB_diversity_metrics.csv", sep=",", quote=F, col.names=NA)

diversity_metrics_F <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/otu_count_FEOB_diversity_metrics.csv")
diversity_metrics_F

write.table(diversity_metrics_F, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/enrich_diversity_metrics_FEOB.csv", sep=",", quote=F, col.names=NA)
```

# Boxplots: FEOB
```{r}
diversity_metrics_F <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/enrich_diversity_metrics_FEOB.csv")
diversity_metrics_F

diversity_metrics_F %>%
  group_by(location) %>%
  pivot_longer(cols=c(richness, shannon, evenness), 
               names_to="metric") %>%
ggplot(aes(x=n, y=value, fill= metal_type)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 8, outlier.size = 2) +
  facet_wrap(~metric, nrow=1, scales="free_y")

anovaS4 <- aov(shannon ~ metal_type, data = diversity_metrics_F)
anovaS5 <- aov(richness ~ metal_type, data = diversity_metrics_F)
anovaS6 <- aov(evenness ~ metal_type, data = diversity_metrics_F)
summary(anovaS4)
summary(anovaS5)
summary(anovaS6)


ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/FEOB/alpha_diversity_metrics_FEOB_metal_location.tiff", width = 10, height = 7)
```

## ANOSIM: FeOB
```{r}
pc_ano <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/nmds_asv_otu_FeOB.csv")
pc_ano

# All Samples: Location
com_ano <- pc_ano[,6:ncol(pc_ano)]
com_ano

m_com_ano <- as.matrix(com_ano)
m_com_ano
```

###: ANOSIM FEOB: Location & Metal Type
```{r}
ano_all <- anosim(m_com_ano, pc_ano$location, distance = "bray", permutations = 9999)
ano_all

#ANOSIM statistic R:  1
#      Significance:  0.16666


# All Samples: Metal Type
com_ano = pc_ano[,6:ncol(pc_ano)]
m_com_ano = as.matrix(com_ano)

ano_all = anosim(m_com_ano, pc_ano$metal_type, distance = "bray", permutations = 9999)
ano_all

#ANOSIM statistic R:  1
#      Significance: 0.1666
```

# SRB Enrichment Only
```{r}
metadata_S <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/metadata_SRB.csv")

otu_counts_S <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/asv_otu_SRB.csv") %>%
  pivot_longer(-ASV, names_to="sample_id", values_to = "count")
otu_counts_S

taxonomy_S <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/asv_tax_SRB.csv")
taxonomy_S

otu_rel_abund_S <- inner_join(metadata_S, otu_counts_S, by="sample_id") %>%
  inner_join(., taxonomy_S, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom","Phylum", "Class", "Order", "Family", "ASV"),
         names_to="level",
         values_to="taxon")
otu_rel_abund_S

write.table(otu_rel_abund_S, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/otu_rel_abund_SRB.csv", sep=",", quote=F, col.names=NA)

otu_rel_abund_S <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/otu_rel_abund_SRB.csv")
otu_rel_abund_S
```

```{r}
## Phylum
otu_rel_abund_S %>%
  filter(level=="Phylum") %>%
  group_by(sample_id, metal_type, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(metal_type, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=metal_type, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=metal_type, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, y="Mean Relative Abundance (%)") + theme_classic()
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/SRB/phylum_stacked_barchart_SRB_metal.tiff", width=15, height=10)

## Class
otu_rel_abund_S %>%
  filter(level=="Class") %>%
  group_by(sample_id, metal_type, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(metal_type, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=metal_type, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=metal_type, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/SRB/class_stacked_barchart_SRB_metal.tiff", width=15, height=10)

## Order
otu_rel_abund_S %>%
  filter(level=="Order") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/SRB/order_stacked_barchart_SRB.tiff", width=15, height=10, limitsize = FALSE)

## Family
otu_rel_abund_S %>%
  filter(level=="Family") %>%
  group_by(sample_id, metal_type, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(metal_type, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=metal_type, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=metal_type, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/SRB/family_stacked_barchart_SRB_metal_type.tiff", width=15, height=10, limitsize = FALSE)

## Genus
otu_rel_abund_S %>%
  filter(level=="Genus") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/SRB/genus_stacked_barchart_SRB.tiff", width=15, height=10, limitsize = FALSE)
```

# Phyloseq: SRB ONLY
```{r}
library(phyloseq)
library(Biostrings)

# Read data into R
otu_tabS <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/asv_otu_SRB.csv")
otu_tabS

tax_tabS <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/asv_tax_SRB.csv")
tax_tabS

samples_dfS <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/metadata_SRB.csv")
samples_dfS

# Phyloseq objects need to have row.names
otu_tabS <- otu_tabS %>%
  tibble::column_to_rownames("ASV")

tax_tabS <- tax_tabS %>%
  tibble::column_to_rownames("ASV")

samples_dfS <- samples_dfS %>%
  tibble::column_to_rownames("sample_id")

# Transform OTU & Tax table into matrices
otu_tabS <- as.matrix(otu_tabS)
tax_tabS <- as.matrix(tax_tabS)

#Transform into Phyloseq Objects
ASV = otu_table(otu_tabS, taxa_are_rows = TRUE)
TAX = tax_table(tax_tabS)
sample_id = sample_data(samples_dfS)
  
psS <- phyloseq(ASV, TAX, sample_id)
psS

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 285 taxa and 33 samples ]
#sample_data() Sample Data:       [ 33 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 285 taxa by 6 taxonomic ranks ]
```

```{r}
# Visualize Data
sample_names(psS)
rank_names(psS)
sample_variables(psS)

# Normalize number of reads in each sample using median sequencing depth
total = median(sample_sums(psS))
standS = function(x, t=total) round(t * (x / sum(x)))
psS = transform_sample_counts(psS, standS)

# Bar graphs based on division
plot_bar(psS, fill = "Family") + 
  geom_bar(aes(fill=Family), stat="identity", position="stack") +
  facet_grid(~metal_type + location)
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/FeOB/metal_location_bar.tiff", width=15, height=10, limitsize = FALSE)
```

# Boxplots Taxa Rank
```{r}
require(dplyr) #required for data table filtering
require(ggplot2) #require for plotting
require(scales) #required for percentage conversion
require(ggthemes) #required for colorblind palette
require(ggpubr) #required for wilcoxon rank sum test statistical comparison between explanatory variables (categorical)

psS_family <- psS %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance (or use ps.ra)
  psmelt() %>%                                         # Melt to long format for easy ggploting
  filter(Abundance > 0.01)                             # Filter out low abundance taxa

write.table(psS_family, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/phyloseq_family_abund_SRB.csv", sep=",", quote=F, col.names=NA)
psS_family

box_plot <- ggplot(psS_family, aes(x=microbe_type, y=Abundance, fill=Family)) +
	geom_boxplot() +
  facet_grid(~metal_type + location) +
	scale_y_continuous(labels=percent) +
	xlab("Location") + 
	ylab("Average Relative Abundance\n Out of Total Sequences") +
	theme(axis.text.x=element_text(angle=60, hjust=1)) + #this will tilt the phyla text so that it looks cleaner +
  scale_colour_colorblind() +
	scale_colour_colorblind(name="Family") #+ #this will make the geom_points the same colors as the geom_boxplot
	#geom_point(position=position_dodge(width=0.75), aes(color=Family),show.legend=FALSE)
box_plot #prints box plot to standard out

#Include this to export this figure as a .tiff publication quality file
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/Taxa_Rank_SRB_Family.tiff", plot=box_plot, device=NULL, path=NULL, scale=1, width=15, height=7, dpi=900, limitsize=TRUE)
```

# SRB Relative Abundance Based on Groups
```{r}
#change LEVEL to what ever taxon level you want, group_by factors you want to compare, and taxon == to taxon of interest
otu_rel_abund_S %>%
  filter(level == "Phylum") %>%
  group_by(metal_type, sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups = "drop") %>%
  group_by(metal_type, taxon) %>%
  summarize(
    mean_rel_abund = 100 * mean(rel_abund),
    se_rel_abund   = 100 * sd(rel_abund) / sqrt(n()),
    n              = n(),
    .groups = "drop"
  ) %>%
  filter(taxon == "Thermodesulfobacteriota")
```

```{r}
df_otu_S <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/df_otu_SRB.csv")
df_otu_S

otu_count_all_S <- df_otu_S %>%
  pivot_longer(-sample_id, names_to = "ASV", values_to = "count")
otu_count_all_S

write.table(otu_count_all_S, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/otu_count_all_SRB.csv", sep=",", quote=F, col.names=NA)

otu_count_all_S <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/otu_count_all_SRB.csv")
otu_count_all_S

otu_count_all_S %>%
  group_by(sample_id) %>%
  mutate(total = sum(count)) %>%
  filter(total > 5000) %>%
  group_by(ASV) %>%
  mutate(total=sum(count)) %>% 
  filter(total != 0) %>%
  as.data.frame()
#Going to set threshold at 5000

nmds_asv_otu_S <- inner_join(metadata_S, df_otu_S, by="sample_id")
nmds_asv_otu_S

write.table(nmds_asv_otu_S, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/nmds_asv_otu_S.csv", sep=",", quote=F, col.names=NA)

pc <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/nmds_asv_otu_S.csv")
pc

#make community matrix: extract columns with ASV information
com <- pc[,6:ncol(pc)]
com

#turn ASV information into a matrix
m_com <- as.matrix(com)

#Run NMDS using Bray-Curtis distance
set.seed(123)
nmdsS <- metaMDS(m_com, distance="bray") #stress = nearly zero 
nmds
plot(nmdsS)

#access the specific points data of the NMDS plot & scores
str(nmdsS)
nmds$points
scores(nmdsS)

#extract NMDS scores
data.scores = as.data.frame(scores(nmdsS))

#add columns to data frame
data.scores$sample_id = pc$sample_id
data.scores$location = pc$location
data.scores$microbe_type = pc$microbe_type
data.scores$metal_type = pc$metal_type

head(data.scores)
```

## NMDS: SRB
```{r}
# By Location
xxS = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = location))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "Wreck Location", y = "NMDS2")
xxS
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/SRB/NMDS_SRB_samples_location.tiff", width = 10, height = 10)

# By Metal Type
xxS = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = metal_type, shape = location))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "Metal", y = "NMDS2")
xxS
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/SRB/NMDS_SRB_samples_metal_location.tiff", width = 10, height = 10)

# By Microbe Type
xxS = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = microbe_type))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "Microbe Classification", y = "NMDS2")
xxS
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/SRB/NMDS_SRB_samples_microbe.tiff", width = 10, height = 10)

# By Metal Type & Location
xx = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(shape = metal_type, colour = location))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "Location", y = "NMDS2")
xx
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/SRB/NMDS_SRB_samples_metal_location_type.tiff", width = 10, height = 10)
```

```{r}
# Shannon, Richness, Eveness: FeOB
otu_count_S <- otu_count_all_S %>%
  group_by(sample_id) %>%
  summarize(richness = richness(count),
            shannon = shannon(count), 
            evenness = shannon/log(richness),
            n=sum(count))

write.table(otu_count_S, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/otu_count_SRB_diversity_metrics.csv", sep=",", quote=F, col.names=NA)

diversity_metrics_S <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/otu_count_SRB_diversity_metrics.csv")
diversity_metrics_S

write.table(diversity_metrics_S, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/enrich_diversity_metrics_SRB.csv", sep=",", quote=F, col.names=NA)
```

# Boxplots: SRB
```{r}
diversity_metrics_S <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/enrich_diversity_metrics_SRB.csv")
diversity_metrics_S

diversity_metrics_S %>%
  group_by(location) %>%
  pivot_longer(cols=c(richness, shannon, evenness), 
               names_to="metric") %>%
ggplot(aes(x=n, y=value, fill= metal_type)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 8, outlier.size = 2) +
  facet_wrap(~metric, nrow=1, scales="free_y") +
  facet_wrap(~metric) +
  stat_compare_means(method = "anova",
                     label = "p.format")

anovaS1 <- aov(shannon ~ metal_type, data = diversity_metrics_S)
anovaS2 <- aov(richness ~ metal_type, data = diversity_metrics_S)
anovaS3 <- aov(evenness ~ metal_type, data = diversity_metrics_S)
summary(anovaS1)
summary(anovaS2)
summary(anovaS3)

ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/FEOB/alpha_diversity_metrics_SRB_metal_location.tiff", width = 7, height = 5)
```

# ANOSIM: SRB
```{r}
pc_ano <- read.csv("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/data/nmds_asv_otu_S.csv")
pc_ano

# All Samples: Location
com_ano <- pc_ano[,6:ncol(pc_ano)]
com_ano

m_com_ano <- as.matrix(com_ano)
m_com_ano

ano_all <- anosim(m_com_ano, pc_ano$location, distance = "bray", permutations = 9999)
ano_all

#ANOSIM statistic R:  0.04279
#      Significance:  0.1732


# All Samples: Metal Type
com_ano = pc_ano[,6:ncol(pc_ano)]
m_com_ano = as.matrix(com_ano)

ano_all = anosim(m_com_ano, pc_ano$metal_type, distance = "bray", permutations = 9999)
ano_all

#ANOSIM statistic R:  -0.07131
#      Significance: 0.7774
```

# Combining Figures
```{r}
install.packages("patchwork")
library(patchwork)
```

# Barcharts
```{r}
# Phylum
p1_Fbar_phylum <- otu_rel_abund_F %>%
  filter(level=="Phylum") %>%
  group_by(sample_id, metal_type, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(metal_type, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=metal_type, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=metal_type, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, y="Mean Relative Abundance (%)") + theme_classic()

p2_Sbar_phylum <- otu_rel_abund_S %>%
  filter(level=="Phylum") %>%
  group_by(sample_id, metal_type, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(metal_type, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=metal_type, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=metal_type, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, y="Mean Relative Abundance (%)") + theme_classic()

(p1_Fbar_phylum / p2_Sbar_phylum) +
   plot_annotation(tag_levels = "A")
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/Combined_Figures/phyum_metal_SRB_FeOB.tiff", width = 10, height = 10)

# Class
p1_Fbar_class <- otu_rel_abund_F %>%
  filter(level=="Class") %>%
  group_by(sample_id, metal_type, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(metal_type, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=metal_type, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=metal_type, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, y="Mean Relative Abundance (%)") + theme_classic()

p2_Sbar_class <- otu_rel_abund_S %>%
  filter(level=="Class") %>%
  group_by(sample_id, metal_type, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(metal_type, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=metal_type, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=metal_type, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, y="Mean Relative Abundance (%)") + theme_classic()

(p1_Fbar_class / p2_Sbar_class) +
   plot_annotation(tag_levels = "A")
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/Combined_Figures/class_metal_SRB_FeOB.tiff", width = 10, height = 10)

# Family
p1_Fbar_family <- otu_rel_abund_F %>%
  filter(level=="Family") %>%
  group_by(sample_id, metal_type, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(metal_type, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=metal_type, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=metal_type, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, y="Mean Relative Abundance (%)") + theme_classic()

p2_Sbar_family <- otu_rel_abund_S %>%
  filter(level=="Family") %>%
  group_by(sample_id, metal_type, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(metal_type, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=metal_type, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=metal_type, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, y="Mean Relative Abundance (%)") + theme_classic()

(p1_Fbar_family / p2_Sbar_family) +
   plot_annotation(tag_levels = "A")
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/Combined_Figures/family_metal_SRB_FeOB.tiff", width = 10, height = 10)
```

# Alpha Diveristy Boxplots
```{r}
p1_FeOB_box <- diversity_metrics_F %>%
  group_by(metal_type) %>%
  pivot_longer(cols=c(richness, shannon, evenness), 
               names_to="metric") %>%
ggplot(aes(x=n, y=value, fill= metal_type)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 8, outlier.size = 2) +
  facet_wrap(~metric, nrow=1, scales="free_y")

p2_SRB_box <- diversity_metrics_S %>%
  group_by(metal_type) %>%
  pivot_longer(cols=c(richness, shannon, evenness), 
               names_to="metric") %>%
ggplot(aes(x=n, y=value, fill= metal_type)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 8, outlier.size = 2) +
  facet_wrap(~metric, nrow=1, scales="free_y")

(p1_FeOB_box / p2_SRB_box) +
   plot_annotation(tag_levels = "A")
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/Combined_Figures/alpha_diversity_metrics_SRB_FeOB_metal.tiff", width = 10, height = 7)
```

# Taxa Rank
```{r}
p1_FeOB_rank <- box_plot <- ggplot(psF_family, aes(x=microbe_type, y=Abundance, fill=Family)) +
	geom_boxplot() +
  facet_grid(~metal_type) +
	scale_y_continuous(labels=percent) +
	xlab("Location") + 
	ylab("Average Relative Abundance\n Out of Total Sequences") +
	theme(axis.text.x=element_text(angle=60, hjust=1)) + #this will tilt the phyla text so that it looks cleaner +
  scale_colour_colorblind() +
	scale_colour_colorblind(name="Family") #+ #this will make the geom_points the same colors as the geom_boxplot
	#geom_point(position=position_dodge(width=0.75), aes(color=Family),show.legend=FALSE)
box_plot #prints box plot to standard out

p2_SRB_rank <- box_plot <- ggplot(psS_family, aes(x=microbe_type, y=Abundance, fill=Family)) +
	geom_boxplot() +
  facet_grid(~metal_type) +
	scale_y_continuous(labels=percent) +
	xlab("Location") + 
	ylab("Average Relative Abundance\n Out of Total Sequences") +
	theme(axis.text.x=element_text(angle=60, hjust=1)) + #this will tilt the phyla text so that it looks cleaner +
  scale_colour_colorblind() +
	scale_colour_colorblind(name="Family") #+ #this will make the geom_points the same colors as the geom_boxplot
	#geom_point(position=position_dodge(width=0.75), aes(color=Family),show.legend=FALSE)
box_plot #prints box plot to standard out
  
(p1_FeOB_rank / p2_SRB_rank) +
   plot_annotation(tag_levels = "A")
ggsave("/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/Combined_Figures/SRB_FeOB_metal_rank.tiff", width = 10, height = 10)
```

# Taxa Rank Year
```{r}
install.packages("magick")
library(magick)

img23 <- image_read("/Users/maggieshostak/Desktop/Saipan_2023_R_Studio_V6V8/Saipan_2023/post_rarefaction/MIC_ONLY/Taxa_Rank_Family_23.tiff")

img25 <- image_read("/Users/maggieshostak/Desktop/Saipan_2025_R_Studio_V6V8/data/MIC_ONLY/Taxa_Rank_Family_25.tiff")

combined <- image_append(c(img23, img25))

image_write(combined, "/Users/maggieshostak/Desktop/Enrichment_Sequencing_2026/results/Combined_Figures/Taxa_Rank_23_25.tiff")
```
