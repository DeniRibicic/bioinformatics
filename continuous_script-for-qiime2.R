library(ggplot2)
library(phyloseq)
library(data.table)
library(ggrepel)
library(viridis)
library(mixOmics)
library (plyr)

science_theme = theme(panel.grid.major = element_line(size = 0.5, color = "grey"), 
                      text = element_text(size = 14)) + 
  theme(axis.title.y = element_text(vjust=2)) +
  theme(panel.background = element_rect(fill = NA))


setwd("/home/denir/Projects/MonMic/Round_Updated/qiime2/exported/")
#otufile="table_w_tax.from_biom.txt"
biomfile="table-with-taxonomy.biom"
mapfile="/home/denir/Projects/MonMic/Round_Updated/qiime2/Mapping-file.txt"
rs_file="dna-sequences.fasta"
trefile="tree.nwk"
#import data into phyloseq as an object (qiimedata) by "import_qiime" command
qiimedata=import_biom(biomfile, trefile, rs_file)
map <- import_qiime_sample_data(mapfile)
#merge the files 
Q <- merge_phyloseq(qiimedata, map, trefile)
#change the taxonomy labels if using other database than Greengenes
colnames(tax_table(Q)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus", "Species")
Q

Q1 <- prune_samples(sample_sums(Q) > 5000, Q)
Q1

#Since we are interested in alpha and  beta diversity, it is probably not a bad idea to prune OTUs that are not present in any of the samples
Q1 <- prune_taxa(taxa_sums(Q1) > 0, Q1)
Q1

Q2 <- transform_sample_counts(Q1, function(x)  x/sum(x)*100)

#only if you want to merge samples-replicates
sample_data(Q2)$unique_sample <- paste0(sample_data(Q2)$Location, sep="_", 
                                        sample_data(Q2)$SamplingPosition, sep="_", 
                                        sample_data(Q2)$TimePoint)
Q3 = merge_samples(Q2, "unique_sample", fun=mean)
Q3 <- transform_sample_counts(Q3, function(x)  x/sum(x)*100)
#retrive original labels for metadata after merging
sample_data(Q3)$Location <- levels(sample_data(Q2)$Location)[get_variable(Q3, "Location")]
sample_data(Q3)$SamplingPosition <- levels(sample_data(Q2)$SamplingPosition)[get_variable(Q3, "SamplingPosition")]
sample_data(Q3)$TimePoint <- levels(sample_data(Q2)$TimePoint)[get_variable(Q3, "TimePoint")]
head(sample_data(Q3))

df <- psmelt(Q3)

#################################################################################

df$Genus <- as.character(df$Genus)
df$Genus[is.na(df$Genus)] <- "Unassigned"
df$Genus <- gsub("D_5__", "", df$Genus)
df$Family <- as.character(df$Family)
df$Family[is.na(df$Family)] <- "Unassigned"
df$Family <- gsub("D_4__", "", df$Family)
df$Order <- as.character(df$Order)
df$Order[is.na(df$Order)] <- "Unassigned"
df$Order <- gsub("D_3__", "", df$Order)

df1 <- df
#df1 <- aggregate(Abundance~Order + SamplingPosition + SampleNumbers + replicate_alt + Location + Date + TimePoint, df1, sum)
df1 <- aggregate(Abundance~Order + Sample + SamplingPosition + Location + TimePoint, df1, sum)

df1 <- as.data.table(df1)
# convert Genus to a character vector from a factor because R
df1$Order <- as.character(df1$Order)
# group dataframe by Genus, calculate median rel. abundance
df1[, max := max(Abundance, na.rm = TRUE), 
    by = "Order"]
# Change name to "Other" of Family less than 1%
df1[(max <= 5), Order := "Other"]
unique(df1$Order)

#df2 <- aggregate(Abundance~Order + SamplingPosition + SampleNumbers + replicate_alt + Location + Date + TimePoint, df1, FUN=sum)
df2 <- aggregate(Abundance~Order + Sample + SamplingPosition + Location + TimePoint, df1, sum)
#write.table(df2, file = "replicates_Genus-natural-10pct-long", row.names = FALSE, sep = "\t", quote = FALSE, dec = ",")
unique(df2$Order)

#df2$taxonomy <- gsub('Unassigned/Unassigned/Unassigned', 'Unassigned', df2$taxonomy)

df2 <- df2[order(df2$Order),]
unique(df2$Order)

df2 <- df2[order(df2$Abundance),]

cols <- c(Chlamydiales="red",
          Betaproteobacteriales="#745700",  
          Aeromonadales="#a77d00", 
          Alteromonadales="#daa300",
          Caldilineales="#c13100",
          Caulobacterales="#007558",
          Corynebacteriales="#290a00",
          Microtrichales="#a5db00",
          Flavobacteriales="#314100",
          Cytophagales="#577400",
          Nitrospirales="#00c131",
          Pseudomonadales="#0090c1",
          Rhodobacterales="#8e006a",
          Francisellales="hotpink",
          Myxococcales="#c10090",
          Micrococcales="#f400b6",
          Oceanospirillales="blue",
          Planctomycetales="green",
          Deinococcales="darkolivegreen2",
          Rhizobiales="pink",
          Pirellulales="darkred",
          Sphingobacteriales="brown2",
          Sphingomonadales="#ffc510",
          Thiotrichales="#dca700",
          Chitinophagales="darkorchid4",
          Verrucomicrobiales="yellow",
          Gemmatales="gray100",
          Gemmatimonadales="gray95",
          Chthoniobacterales="gold3",
          Xanthomonadales="gold",
          Other="gray41",
          Unassigned="gray90")

df2$Order = factor(df2$Order, levels=c("Chlamydiales",
                                       "Betaproteobacteriales",  
                                       "Aeromonadales", 
                                       "Alteromonadales",
                                       "Caldilineales",
                                       "Caulobacterales",
                                       "Corynebacteriales",
                                       "Microtrichales",
                                       "Flavobacteriales",
                                       "Cytophagales",
                                       "Nitrospirales",
                                       "Pseudomonadales",
                                       "Rhodobacterales",
                                       "Francisellales",
                                       "Myxococcales",
                                       "Micrococcales",
                                       "Oceanospirillales",
                                       "Planctomycetales",
                                       "Deinococcales",
                                       "Rhizobiales",
                                       "Pirellulales",
                                       "Sphingobacteriales",
                                       "Sphingomonadales",
                                       "Thiotrichales",
                                       "Chitinophagales",
                                       "Verrucomicrobiales",
                                       "Gemmatales",
                                       "Gemmatimonadales",
                                       "Xanthomonadales",
                                       "Chthoniobacterales",
                                       "Other",
                                       "Unassigned"))

df2$TimePoint = factor(df2$TimePoint, levels=c("t0", "t1", "t2", "t3", "t4",
                                               "t5", "t6", "t7", "t8", "t9",
                                               "t10", "t11", "t12", "t13", "t14"))

#df2$pos_rep <- paste0(df2$SamplingPosition, sep="-", df2$replicate_alt)
p5 <- ggplot(subset(df2,Location %in% "E"), aes(x=SamplingPosition, y=Abundance, fill=Order, )) + #remember to change location
  geom_bar(stat = "identity", color="black") +
  scale_fill_manual(values = cols, name="Taxonomical composition\nby Order level") +
  theme(legend.position="right") +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("Relative abundance (% of total sequences)") +
  xlab("Sample") +
  facet_grid(~TimePoint) +
  ggtitle("Farm E") + #remember to change the title for each plot 
  science_theme
p5

#save 2200x450
#name Farm-C_microbiota-t0-t9

######Alternatively plotting after Gaute Nedberg suggestion##############

p5 <- ggplot(subset(df2,Location %in% "E"), aes(x=TimePoint, y=Abundance, fill=Order, )) + #remember to change location
  geom_bar(stat = "identity", color="black") +
  scale_fill_manual(values = cols, name="Taxonomical composition\nby Order level") +
  theme(legend.position="right") +
  ylab("Relative abundance (% of total sequences)") +
  xlab("Sample") +
  facet_wrap(~SamplingPosition, nrow=2) +
  ggtitle("Farm E") + #remember to change the title for each plot 
  science_theme
p5

#save 1400x613
#name:Farm_A-microbiota2-t0-t9
#########################################################################################

#quick overview without legend and specific colors
p6 <- ggplot(df2, aes(x=pos_rep, y=Abundance, fill=Order, )) + 
  geom_bar(stat = "identity", color="black") +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("Relative abundance (% of total sequences)") +
  xlab("Replicate") +
  facet_grid(Location~TimePoint, scale="free_x") +
  ggtitle("") + 
  science_theme
p6

##All with legend
p7 <- ggplot(df2, aes(x=pos_rep, y=Abundance, fill=Order, )) + #remember to change location
  geom_bar(stat = "identity", color="black") +
  scale_fill_manual(values = cols, name="Taxonomical composition\nby Order level") +
  theme(legend.position="right") +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("Relative abundance (% of total sequences)") +
  xlab("Replicate") +
  facet_grid(Location~TimePoint, scale="free_x", space="free") +
  ggtitle("All Farms") + #remember to change the title for each plot 
  science_theme
p7

##################################################################
##HEATMAPS

ht <- df
ht$taxonomy <- paste(ht$Order, ht$Family, ht$Genus, sep = "/")
ht1 <- aggregate(Abundance~taxonomy + SamplingPosition + SampleNumbers + replicate_alt + Location + Date + TimePoint, ht, sum)

ht2 <- as.data.table(ht1)
# convert Genus to a character vector from a factor because R
ht2$taxonomy <- as.character(ht2$taxonomy)
# group dataframe by Genus, calculate median rel. abundance
ht2[, max := max(Abundance, na.rm = TRUE), 
    by = "taxonomy"]
# Change name to "Other" of Family less than 1%
ht2[(max <= 1), taxonomy := "Other"]
unique(ht2$taxonomy)

ht3 <- aggregate(Abundance~taxonomy + SamplingPosition + SampleNumbers + replicate_alt + Location + Date + TimePoint, ht2, FUN=sum)
#write.table(df2, file = "replicates_Genus-natural-10pct-long", row.names = FALSE, sep = "\t", quote = FALSE, dec = ",")
unique(ht3$taxonomy)

ht3$logAbundance <- log10(ht3$Abundance+1)

ht4 <- ht3
ht5 <- subset(ht4, ! taxonomy %in% c("Other", "Unassigned/Unassigned/Unassigned"))

ht<- ggplot(subset(ht5, Location %in% "A"), aes(x = SamplingPosition, y = taxonomy)) + 
  geom_tile(aes(fill = logAbundance)) + 
  scale_fill_gradient(name = 'Log scale\nrelative abundance', low = 'darkblue', high = 'yellow') + 
  ylab("Taxonomy (Order/Family/Genus)") +
  xlab("Sampling position ") +
  facet_grid(~TimePoint) +
  ggtitle("Farm A")
ht
#save 2000x2500
#Farm-C_heatmap_t0-t9