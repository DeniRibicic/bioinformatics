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

#otufile="table_w_tax.from_biom.txt"
biomfile="new_otu_table_json.biom"
mapfile="Mapping-file.txt"
rs_file="new_refseqs.fna"
trefile="rep_set.tre"

setwd("/media/Storage/denir/Projects/MonMic/Round_B/only_qc/7_otus/")
#import data into phyloseq as an object (qiimedata) by "import_qiime" command
qiimedata=import_biom(biomfile, trefile, rs_file)
map <- import_qiime_sample_data(mapfile)
#merge the files 
Q <- merge_phyloseq(qiimedata, map, trefile)
#change the taxonomy labels if using other database than Greengenes
colnames(tax_table(Q)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus", "Species")
Q

setwd("/media/Storage/denir/Projects/MonMic/Round_A/7_otus/")
#import data into phyloseq as an object (qiimedata) by "import_qiime" command
qiimedata=import_biom(biomfile, trefile, rs_file)
map <- import_qiime_sample_data(mapfile)
#merge the files 
P <- merge_phyloseq(qiimedata, map, trefile)
#change the taxonomy labels if using other database than Greengenes
colnames(tax_table(P)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus", "Species")
P

setwd("/media/Storage/denir/Projects/MonMic/Round_C/7_otus/")
#import data into phyloseq as an object (qiimedata) by "import_qiime" command
qiimedata=import_biom(biomfile, trefile, rs_file)
map <- import_qiime_sample_data(mapfile)
#merge the files 
C <- merge_phyloseq(qiimedata, map, trefile)
#change the taxonomy labels if using other database than Greengenes
colnames(tax_table(C)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus", "Species")
C

setwd("/media/Storage/denir/Projects/MonMic/Round_D/7_otus/")
#import data into phyloseq as an object (qiimedata) by "import_qiime" command
qiimedata=import_biom(biomfile, trefile, rs_file)
map <- import_qiime_sample_data(mapfile)
#merge the files 
D <- merge_phyloseq(qiimedata, map, trefile)
#change the taxonomy labels if using other database than Greengenes
colnames(tax_table(D)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus", "Species")
D

Q2 = transform_sample_counts(Q, function(x)  x/sum(x)*100)
P2 = transform_sample_counts(P, function(x)  x/sum(x)*100)
C2 = transform_sample_counts(C, function(x)  x/sum(x)*100)
D2 = transform_sample_counts(D, function(x)  x/sum(x)*100)

dfA <- psmelt(Q2)
dfB <- psmelt(P2)
dfC <- psmelt(C2)
dfD <- psmelt(D2)

df <- rbind(dfA, dfB, dfC, dfD)

dfm <- df

#Now we have total dataset already analyzed in QIIME, so previous inputs can be avoided

setwd("/media/Storage/denir/Projects/MonMic/Round_Updated/7_otus/")
#import data into phyloseq as an object (qiimedata) by "import_qiime" command
qiimedata=import_biom(biomfile, trefile, rs_file)
map <- import_qiime_sample_data(mapfile)
#merge the files 
Q <- merge_phyloseq(qiimedata, map, trefile)
#change the taxonomy labels if using other database than Greengenes
colnames(tax_table(Q)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus", "Species")
Q

Q1 <- transform_sample_counts(Q, function(x)  x/sum(x)*100)
df <- psmelt(Q1)

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
df1 <- aggregate(Abundance~Order + SamplingPosition + SampleNumbers + replicate_alt + Location + Date + TimePoint, df1, sum)

df1 <- as.data.table(df1)
# convert Genus to a character vector from a factor because R
df1$Order <- as.character(df1$Order)
# group dataframe by Genus, calculate median rel. abundance
df1[, max := max(Abundance, na.rm = TRUE), 
  by = "Order"]
# Change name to "Other" of Family less than 1%
df1[(max <= 5), Order := "Other"]
unique(df1$Order)

df2 <- aggregate(Abundance~Order + SamplingPosition + SampleNumbers + replicate_alt + Location + Date + TimePoint, df1, FUN=sum)
#write.table(df2, file = "replicates_Genus-natural-10pct-long", row.names = FALSE, sep = "\t", quote = FALSE, dec = ",")
unique(df2$Order)

#df2$taxonomy <- gsub('Unassigned/Unassigned/Unassigned', 'Unassigned', df2$taxonomy)

df2 <- df2[order(df2$Order),]
unique(df2$Order)

df2 <- df2[order(df2$Abundance),]

cols <- c(Acidimicrobiales="red",
          Burkholderiales="#745700",
          Aeromonadales="#a77d00",
          Alteromonadales="#daa300",
          Caldilineales="#c13100",
          Caulobacterales="#007558",
          Corynebacteriales="#290a00",
          Deinococcales="#a5db00",
          Flavobacteriales="#314100",
          Cytophagales="#577400",
          Nitrospirales="#00c131",
          Nitrosomonadales="#f43e00",
          Pseudomonadales="#0090c1",
          Rhodobacterales="#8e006a",
          Frankiales="#c10090",
          Micrococcales="#f400b6",
          Oceanospirillales="blue",
          Planctomycetales="green",
          Rhizobiales="pink",
          Rhodospirillales="darkred",
          Rickettsiales="cornsilk1",
          Sphingobacteriales="brown2",
          Sphingomonadales="#ffc510",
          Thiotrichales="#dca700",
          `Subgroup 6`="darkorchid4",
          `uncultured bacterium`="pink1",
          Verrucomicrobiales="yellow",
          Vibrionales="gray100",
          Xanthomonadales="gold",
          Other="gray41",
          Unassigned="gray90")

df2$Order = factor(df2$Order, levels=c("Acidimicrobiales",
                      "Burkholderiales",
                      "Aeromonadales",
                      "Alteromonadales",
                      "Caldilineales",
                      "Caulobacterales",
                      "Corynebacteriales",
                      "Deinococcales",
                      "Flavobacteriales",
                      "Cytophagales",
                      "Nitrospirales",
                      "Nitrosomonadales",
                      "Pseudomonadales",
                      "Rhodobacterales",
                      "Frankiales",
                      "Micrococcales",
                      "Oceanospirillales",
                      "Planctomycetales",
                      "Rhizobiales",
                      "Rhodospirillales",
                      "Rickettsiales",
                      "Sphingobacteriales",
                      "Sphingomonadales",
                      "Thiotrichales",
                      "Subgroup 6",
                      "uncultured bacterium",
                      "Verrucomicrobiales",
                      "Vibrionales",
                      "Xanthomonadales",
                      "Other",
                      "Unassigned"))

df2$TimePoint = factor(df2$TimePoint, levels=c("t0", "t1", "t2", "t3", "t4", "t5", "t6", "t7"))
  
df2$pos_rep <- paste0(df2$SamplingPosition, sep="-", df2$replicate_alt)
p5 <- ggplot(subset(df2,Location %in% "C"), aes(x=pos_rep, y=Abundance, fill=Order, )) + #remember to change location
  geom_bar(stat = "identity", color="black") +
  scale_fill_manual(values = cols, name="Taxonomical composition\nby Order level") +
  theme(legend.position="right") +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("Relative abundance (% of total sequences)") +
  xlab("Replicate") +
  facet_grid(~TimePoint) +
  ggtitle("Farm C") + #remember to change the title for each plot 
  science_theme
p5

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