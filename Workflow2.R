setwd("/home/denir/Projects/MonMic/Round_G/exported/")

library(ggplot2)
library(phyloseq)
library(data.table)
library(ggrepel)
library(dplyr)
library(viridis)
library(vegan)
library(metagenomeSeq)


science_theme = theme(panel.grid.major = element_line(size = 0.5, color = "grey"), 
                      text = element_text(size = 14)) + 
  theme(axis.title.y = element_text(vjust=2), axis.title.x =element_blank()) + 
  theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
  theme(panel.background = element_rect(fill = NA))

#set alternate theme for scietnific graphics
science_theme2 = theme(panel.grid.major = element_line(size = 0.5, color = "grey"), 
                       text = element_text(size = 14)) + 
  theme(axis.title.y = element_text(vjust=2)) +
  theme(panel.background = element_rect(fill = NA))

science_theme3 = theme(panel.grid.major = element_blank(), 
                       text = element_text(size = 14)) + 
  theme(axis.title.y = element_text(vjust=2), axis.title.x =element_blank()) + 
  theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
  theme(panel.background = element_rect(fill = NA))

#otufile="table_w_tax.from_biom.txt"
biomfile="new_otu_table_json.biom"
mapfile="Mapping-file.txt"
rs_file="new_refseqs.fna"
trefile="rep_set.tre"
#or this if qiime2 converted files used
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

#Let's explore sequencing depth of our data (all reads that passed QC checks in QIIME):
#explore read counts
sdt = data.table(as(sample_data(Q), "data.frame"),
TotalReads = sample_sums(Q), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth

#Majority of the samples consist of about 50k reads. Two samples have more than 400K reads
#prefiltering; remove samples with low number of counts (reads) based on sequencing depth histogram
Q1 <- prune_samples(sample_sums(Q) > 5000, Q)
Q1

#Since we are interested in alpha and  beta diversity, it is probably not a bad idea to prune OTUs that are not present in any of the samples
Q1 <- prune_taxa(taxa_sums(Q1) > 0, Q1)
Q1
#No sample was filtered out based on low number of reads.

#Explore alpha-diveristy (indicies based)
#plot alpha diveristy Q1:

sample_data(Q1)$SampleNumbers <- as.factor(sample_data(Q1)$SampleNumbers)
                                      
p <- plot_richness(Q1, x = "SamplingPosition", color = "SamplingPosition", shape = "Location", measures = c("Shannon", "Observed")) + 
  geom_point(size=4) + 
  theme_bw(base_size = 12, base_family = "Helvetica") +
  facet_wrap(~TimePoint)
p <- p + science_theme2
p

#####For each Farm separately
A <- subset_samples(Q1, Location == "A")
B <- subset_samples(Q1, Location == "B")
C <- subset_samples(Q1, Location == "C" & TimePoint %in% c("t12","t13"))
#seems that Farm C has additional sampling date, that why to select specific ones that you want 
D <- subset_samples(Q1, Location == "D")
E <- subset_samples(Q1, Location == "E")

a <- plot_richness(A, x = "SamplingPosition", color = "SamplingPosition", shape="TimePoint", measures = c("Shannon", "Observed")) + 
  geom_point(size=4) + 
  theme_bw(base_size = 12, base_family = "Helvetica") +
  ggtitle("Farm A")
a <- a + science_theme2
a

b <- plot_richness(B, x = "SamplingPosition", color = "SamplingPosition", shape="TimePoint", measures = c("Shannon", "Observed")) + 
  geom_point(size=4) + 
  theme_bw(base_size = 12, base_family = "Helvetica") +
  ggtitle("Farm B")
b <- b + science_theme2
b

c <- plot_richness(C, x = "SamplingPosition", color = "SamplingPosition", shape="TimePoint", measures = c("Shannon", "Observed")) + 
  geom_point(size=4) + 
  theme_bw(base_size = 12, base_family = "Helvetica") +
  ggtitle("Farm C")
c <- c + science_theme2
c

d <- plot_richness(D, x = "SamplingPosition", color = "SamplingPosition", shape="TimePoint", measures = c("Shannon", "Observed")) + 
  geom_point(size=4) + 
  theme_bw(base_size = 12, base_family = "Helvetica") +
  ggtitle("Farm D")
d <- d + science_theme2
d

e <- plot_richness(E, x = "SamplingPosition", color = "SamplingPosition", shape="TimePoint", measures = c("Shannon", "Observed")) + 
  geom_point(size=4) + 
  theme_bw(base_size = 12, base_family = "Helvetica") +
  ggtitle("Farm E")
e <- e + science_theme2
e
######
#Let's have a look a those rarefaction curves:
#Make rarefaction curves
psdata <- Q1
psdata
sample_sums(psdata)

#create function
calculate_rarefaction_curves <- function(psdata, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}
####

####
rarefaction_curve_data <- calculate_rarefaction_curves(psdata, c('Observed', 'Shannon'), rep(c(1, 10, 100, 1000, 1:100 * 10000), each = 10))
summary(rarefaction_curve_data)

#summarize alpha_diversity
rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))
#additional step when you have sample names starting with number. For some reason previous steps add an prefix X. 
#We need to get rid of this
rarefaction_curve_data_summary$Sample <- gsub("X", "", paste(rarefaction_curve_data_summary$Sample))

#add sample data
rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(psdata)), by.x = 'Sample', by.y = 'row.names')
#####

#####
#plot rarefactions
ggplot(
  data = rarefaction_curve_data_summary_verbose,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    colour = SamplingPosition,
    group = Sample
  )
) + geom_line(
) + geom_pointrange(
) + facet_wrap(
  facets = ~ Measure,
  scales = 'free_y'
)

#just to run farms A, B and C on timepoints t8 and t9
abc <- subset_samples(Q1, Location %in% c("A","B","C") & TimePoint %in% c("t12","t13"))
abc

sample_data(abc)$sampl_rep <- paste0(sample_data(abc)$SamplingPosition, sep="-", sample_data(abc)$SampleNumbers)
ordu2 = ordinate(abc, "PCoA", "unifrac", weighted = TRUE)
p2 <- plot_ordination(abc, ordu2, color = "Location", shape = "SamplingPosition") + 
  geom_point(size=4) +
  geom_text_repel(size = 6, aes(label=TimePoint), show.legend = FALSE) + 
  geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0)) +
  theme_bw(base_size = 12, base_family = "Helvetica")
  p2

p2 <- p2 + 
  stat_ellipse(type = "norm", level = .95, aes(group=Location)) +
  #stat_ellipse(type = "t") +
  theme_bw()
p2

#Run normal data without subsetting
#What about beta-diversity:
#Make PCoA on Weighted Unifrac
sample_data(Q1)$sampl_rep <- paste0(sample_data(Q1)$SamplingPosition, sep="-", sample_data(Q1)$SampleNumbers)
ordu2 = ordinate(Q1, "PCoA", "unifrac", weighted = TRUE)
p2 <- plot_ordination(Q1, ordu2, color = "Location", shape = "SamplingPosition") + 
geom_point(size=4) +
#geom_text_repel(size = 6, aes(label=TimePoint), show.legend = FALSE) + 
geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0)) +
theme_bw(base_size = 12, base_family = "Helvetica")
p2

p2 <- p2 + 
  stat_ellipse(type = "norm", level = .95, aes(group=Location)) +
  #stat_ellipse(type = "t") +
  theme_bw()
p2

###########################################################################
###########################################################################
#new function instead of rarefying
mseq <- phyloseq_to_metagenomeSeq(Q1) #transform physeq to metagenomeseq

p = cumNormStatFast(mseq)
mseq = cumNorm(mseq, p = p)

otu_table(mseq) <- otu_table(mseq, taxa_are_rows = TRUE)


myround <- function(x) { trunc(x + 0.5) }

scale_reads <- function(physeq, n = min(sample_sums(physeq)), round = "floor") {
  
  # transform counts to n
  physeq.scale <- transform_sample_counts(physeq, 
                                          function(x) {(n * x/sum(x))}
  )
  
  # Pick the rounding functions
  if (round == "floor"){
    otu_table(physeq.scale) <- floor(otu_table(physeq.scale))
  } else if (round == "round"){
    otu_table(physeq.scale) <- myround(otu_table(physeq.scale))
  }
  
  # Prune taxa and return new phyloseq object
  physeq.scale <- prune_taxa(taxa_sums(physeq.scale) > 0, physeq.scale)
  return(physeq.scale)
}

#scale
Q2 <- Q1 %>% scale_reads(round = "round") 

#stat test
ordu3 <- distance(Q2, "wunifrac")
dfSD <- as(sample_data(Q2), "data.frame")

groups <- dfSD[["TimePoint"]]
groups

temp_adonis_45 <- adonis(ordu3 ~ SamplingPosition*TimePoint, data = dfSD)
temp_adonis_45

mod <- betadisper(ordu3, groups)
permutest(mod)

#This output tells us that our adonis test is significant so we can reject the null hypothesis that our 
#TimePoints and SamplingPositions have the same centroid.

#Additionally, our betadisper results are significant, meaning we cannot reject the null hypothesis that our 
#groups have the same dispersions. This means we can be more confident that our adonis result is a 
#real result, and not due to differences in group dispersions.

#the dispersion is different between groups, then examine
plot(mod)
boxplot(mod)
mod.HSD <- TukeyHSD(mod)
mod.HSD
plot(mod.HSD)

####################################################################################
###################################################################################
# C (Bremnes) merging Mapping file with water quality (env. parameteres)
setwd("/home/denir/Projects/MonMic/water_quality/")
df <- read.table('dfm_bremnes.csv',header=T,sep=',',dec='.',fileEncoding="UTF-8",na.strings = '')
df = na.omit(df)
df1 <- subset(df, Dato %in% c("2017-11-06", "2017-11-20", 
                              "2017-12-04", "2017-12-18", 
                              "2018-01-02", "2018-01-15"))

df1$Dato <- gsub(df1$Dato, pattern = "2017-11-06", replacement = "t0")
df1$Dato <- gsub(df1$Dato, pattern = "2017-11-20", replacement = "t1")
df1$Dato <- gsub(df1$Dato, pattern = "2017-12-04", replacement = "t2")
df1$Dato <- gsub(df1$Dato, pattern = "2017-12-18", replacement = "t3")
df1$Dato <- gsub(df1$Dato, pattern = "2018-01-02", replacement = "t4")
df1$Dato <- gsub(df1$Dato, pattern = "2018-01-15", replacement = "t5")

df2 <- subset(df1, variable %in% c("Biomasse_tonn", "Dødelighet_‰", "O3_Styrke",
                                   "Lys regime", "Alkalitet_CaCO3_mg/L", "NO3_mg/L",
                                   "NO2_mg/L", "NH3_mg/L", "CO2 mg/L", "O2 mg/L",
                                   "Temp_Filter_ ⁰C", "Temp_Spedevann_ ⁰C"))

df3 <- subset(df2, select = c("Dato", "variable", "value"))

#an alternative to dcast# df4 <- reshape(df3, idvar = "Dato", timevar = "variable", direction = "wide", )
df4 <- dcast(df3, Dato ~ variable)

setwd("/media/Storage/denir/Projects/MonMic/Round_Updated/7_otus/")
map <- read.delim2("Mapping-file.txt")
map2 <- subset(map, Location %in% "C")

map_new <- merge(x = map2, y = df4, by.x = "TimePoint", by.y = "Dato")
#write.table(map_new, file="Mapping-file_C.txt", quote = F, sep = "\t", row.names = F)

###########
q2 <- phyloseq::distance(physeq = Q2, method = "wunifrac")

cap_ord <- ordinate(physeq = Q2, method = "CAP", distance = q2,
                    formula = ~ Biomasse_fish + Temp_Filter + O2 + O3 + Alkalitet_CaCO3)

cap_plot <- plot_ordination(Q2, cap_ord, color = "TimePoint", shape = "SamplingPosition", axes = c(1,2)) + 
  geom_point(size=4) +
  geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0)) +
  theme_bw(base_size = 12, base_family = "Helvetica")

# Now add the environmental variables as arrows
arrowmat = vegan::scores(cap_ord, display = "bp")
# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
# Define the arrow aesthetic mapping
arrow_map = aes(xend = CAP1, yend = CAP2, x = 0, y = 0, shape = NULL, color = NULL, 
                label = labels)
label_map = aes(x = 1.2 * CAP1, y = 1.2 * CAP2, shape = NULL, color = NULL, 
                label = labels)

arrowhead = arrow(length = unit(0.05, "npc"))

# Make a new graphic
cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  ) +
  ggtitle("Canonical Analysis of Principal Coordinates; Farm C")

anova(cap_ord)

#####################################################
#########################################################################################
# A (Laksefjor) merging Mapping file with water quality (env. parameteres)
setwd("/media/Storage/denir/Projects/MonMic/water_quality/")
df <- read.table('dfm_laksefjord.csv',header=T,sep=',',dec='.',fileEncoding="UTF-8",na.strings = '')
df = na.omit(df)
df1 <- subset(df, Dato %in% c("2017-11-07", "2017-11-20", 
                              "2017-12-05", "2017-12-18", 
                              "2018-01-02", "2018-01-17"))

df1$Dato <- gsub(df1$Dato, pattern = "2017-11-07", replacement = "t0")
df1$Dato <- gsub(df1$Dato, pattern = "2017-11-20", replacement = "t1")
df1$Dato <- gsub(df1$Dato, pattern = "2017-12-05", replacement = "t2")
df1$Dato <- gsub(df1$Dato, pattern = "2017-12-18", replacement = "t3")
df1$Dato <- gsub(df1$Dato, pattern = "2018-01-02", replacement = "t4")
df1$Dato <- gsub(df1$Dato, pattern = "2018-01-17", replacement = "t5")

df2 <- subset(df1, variable %in% c("Biomasse_tonn", "Dødelighet_‰", "O3_Styrke",
                                   "Lys regime", "Alkalitet_CaCO3_mg/L", "NO3_mg/L",
                                   "NO2_mg/L", "NH3_mg/L", "CO2 mg/L", "O2 mg/L",
                                   "Temp_Filter_ ⁰C", "Temp_Spedevann_ ⁰C"))

df3 <- subset(df2, select = c("Dato", "variable", "value"))

#an alternative to dcast# df4 <- reshape(df3, idvar = "Dato", timevar = "variable", direction = "wide", )
df4 <- dcast(df3, Dato ~ variable)

setwd("/media/Storage/denir/Projects/MonMic/Round_Updated/7_otus/")
map <- read.delim2("Mapping-file.txt")
map2 <- subset(map, Location %in% "C")

map_new <- merge(x = map2, y = df4, by.x = "TimePoint", by.y = "Dato")
#write.table(map_new, file="Mapping-file_C.txt", quote = F, sep = "\t", row.names = F)

######
q2 <- phyloseq::distance(physeq = Q2, method = "wunifrac")

cap_ord <- ordinate(physeq = Q2, method = "CAP", distance = q2,
                    formula = ~ Biomasse_fish + Temp_Filter + O2 + O3 + Alkalitet_CaCO3)

cap_plot <- plot_ordination(Q2, cap_ord, color = "TimePoint", shape = "SamplingPosition", axes = c(1,2)) + 
  geom_point(size=4) +
  geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0)) +
  theme_bw(base_size = 12, base_family = "Helvetica")

# Now add the environmental variables as arrows
arrowmat = vegan::scores(cap_ord, display = "bp")
# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
# Define the arrow aesthetic mapping
arrow_map = aes(xend = CAP1, yend = CAP2, x = 0, y = 0, shape = NULL, color = NULL, 
                label = labels)
label_map = aes(x = 1.2 * CAP1, y = 1.2 * CAP2, shape = NULL, color = NULL, 
                label = labels)

arrowhead = arrow(length = unit(0.05, "npc"))

# Make a new graphic
cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  ) +
  ggtitle("Canonical Analysis of Principal Coordinates; Farm C")

anova(cap_ord)

####################################################################################
####################################################################################
# B (Lerøy Belsvik) merging Mapping file with water quality (env. parameteres)
setwd("/media/Storage/denir/Projects/MonMic/water_quality/")
df <- read.table('dfm_belsvik_ST21.csv',header=T,sep=',',dec='.',fileEncoding="UTF-8",na.strings = '')
df = na.omit(df)
df1 <- subset(df, Dato %in% c("2017-11-06", "2017-11-20", 
                              "2017-12-04", "2017-12-18", 
                              "2018-01-02", "2018-01-15"))

df1$Dato <- gsub(df1$Dato, pattern = "2017-11-06", replacement = "t0")
df1$Dato <- gsub(df1$Dato, pattern = "2017-11-20", replacement = "t1")
df1$Dato <- gsub(df1$Dato, pattern = "2017-12-04", replacement = "t2")
df1$Dato <- gsub(df1$Dato, pattern = "2017-12-18", replacement = "t3")
df1$Dato <- gsub(df1$Dato, pattern = "2018-01-02", replacement = "t4")
df1$Dato <- gsub(df1$Dato, pattern = "2018-01-15", replacement = "t5")

df2 <- subset(df1, variable %in% c("Antall døde", "Temperatur Snitt [°C]", "Ammonium Snitt [NH4-N mg/l]",
                                   "Nitrat Snitt [NO3- mg/l]", "Nitrit Snitt [NO2- mg/l]",
                                   "Saltinnhold, promille Snitt [‰]", "Vannforbruk (m3/døgn)"))

df3 <- subset(df2, select = c("Dato", "variable", "value"))

#an alternative to dcast# df4 <- reshape(df3, idvar = "Dato", timevar = "variable", direction = "wide", )
df4 <- dcast(df3, Dato ~ variable)

setwd("/media/Storage/denir/Projects/MonMic/Round_Updated/7_otus/")
map <- read.delim2("Mapping-file.txt")
map2 <- subset(map, Location %in% "B")

map_new <- merge(x = map2, y = df4, by.x = "TimePoint", by.y = "Dato")
#write.table(map_new, file="Mapping-file_Farm-B.txt", quote = F, sep = "\t", row.names = F)

####
q2 <- phyloseq::distance(physeq = Q2, method = "wunifrac")

cap_ord <- ordinate(physeq = Q2, method = "CAP", distance = q2,
                    formula = ~ Ammonium + Mortality + Nitrat + Nitrit + Temperatur)

cap_plot <- plot_ordination(Q2, cap_ord, color = "TimePoint", shape = "SamplingPosition", axes = c(1,2)) + 
  geom_point(size=4) +
  geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0)) +
  theme_bw(base_size = 12, base_family = "Helvetica")

# Now add the environmental variables as arrows
arrowmat = vegan::scores(cap_ord, display = "bp")
# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
# Define the arrow aesthetic mapping
arrow_map = aes(xend = CAP1, yend = CAP2, x = 0, y = 0, shape = NULL, color = NULL, 
                label = labels)
label_map = aes(x = 1.2 * CAP1, y = 1.2 * CAP2, shape = NULL, color = NULL, 
                label = labels)

arrowhead = arrow(length = unit(0.05, "npc"))

# Make a new graphic
cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  ) +
  ggtitle("Canonical Analysis of Principal Coordinates; Farm B")

anova(cap_ord)

#######################################################################
#######################################################################
#Normalize abundances within samples
Q2 = transform_sample_counts(Q1, function(x)  x/sum(x)*100)

#for Krona
tax <- tax_glom(Q2, "Species")
tax2 <- psmelt(Q2)
tax3 <- subset(tax2, select=c("Abundance" ,"SamplingPosition", "SampleNumbers", "replicate_alt", "Location", "TimePoint", "Date", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
#write.table(tax3, file = "total_dataset-D", row.names = FALSE, sep = "\t", quote = FALSE, dec = ",") #for Trond
#For Roman
krona <- tax2
krona$key <- paste(krona$Location ,krona$SamplingPosition, krona$TimePoint, krona$replicate_alt, sep = "_")

##VERY DANGEROUS TO AGGREGATE IF NAs in the dataset- risk of loosing info!!!
#it's okay to aggregate if previously combined columns, given no NAs alone!!

#krona <- aggregate(Abundance~ + SamplingPosition + SampleNumbers + replicate_alt + Location + 
 #                    Date + TimePoint + Kingdom + Phylum + Class + Order + Family + 
  #                   Genus + Species + key, krona, sum)
#krona <- krona[order(krona$Abundance),]

krona2 <- subset(krona, select = c("Abundance", "key", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
krona2 <- krona2[!krona2$Abundance == "0", ] #exclude all raws containing abundance 0

krona2 <- krona2[order(krona2$Abundance),]

#create a list containing dataframes separated by condition
setwd("/media/Storage/denir/Projects/MonMic/analysis_output/Round_D/Krona/")

sample_seq <- krona2$key
for(i in sample_seq)
{
  print (i)
  df <- subset(krona2, krona2$key==i)
  filename <- i
  df$key <- NULL
  write.table(df,filename,sep='\t',row.names = F, quote = FALSE)
}

#batch renaming in linux to resamble Roman's output
#for file in A_*; do  mv "$file" "${file/A_/}";  done #remove index for fishfarm

#ktImportText -o A/A-krona.html A/*

##################################################################
##################################################################old batch sample renaming
sample_seq <- c(1:108)
for(i in 1:length(sample_seq))
{
  # sample[[i]]<- subset(krona, krona$SampleNumbers==sample_seq)
  print (i)
  df <- subset(krona, krona$SampleNumbersII==i)
  filename <- paste('sample_',i,".txt", sep="")
  df$SampleNumbers <- NULL
  write.table(df,filename,sep='\t',row.names = F, quote = FALSE)
}
###############################################################################################
##############################################################################################

#agregate to Genus/Family/Order level
df <- psmelt(Q2)
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
df1$taxonomy <- paste(df1$Order, df1$Family, df1$Genus, sep = "/")
df1 <- aggregate(Abundance~taxonomy + SamplingPosition + SampleNumbers + replicate_alt + Location + Date + TimePoint, df1, sum)

A1= subset(df1, Location %in% c("A"), drop = TRUE)
B1= subset(df1, Location %in% c("B"), drop = TRUE)
C1= subset(df1, Location %in% c("C"), drop = TRUE)

A1 <- as.data.table(A1)
# convert Genus to a character vector from a factor because R
A1$taxonomy <- as.character(A1$taxonomy)
# group dataframe by Genus, calculate median rel. abundance
A1[, max := max(Abundance, na.rm = TRUE), 
    by = "taxonomy"]
# Change name to "Other" of Family less than 1%
A1[(max <= 5), taxonomy := "Other"]
unique(A1$taxonomy)

B1 <- as.data.table(B1)
# convert Genus to a character vector from a factor because R
B1$taxonomy <- as.character(B1$taxonomy)
# group dataframe by Genus, calculate median rel. abundance
B1[, max := max(Abundance, na.rm = TRUE), 
  by = "taxonomy"]
# Change name to "Other" of Family less than 1%
B1[(max <= 5), taxonomy := "Other"]
unique(B1$taxonomy)

C1 <- as.data.table(C1)
# convert Genus to a character vector from a factor because R
C1$taxonomy <- as.character(C1$taxonomy)
# group dataframe by Genus, calculate median rel. abundance
C1[, max := max(Abundance, na.rm = TRUE), 
  by = "taxonomy"]
# Change name to "Other" of Family less than 1%
C1[(max <= 5), taxonomy := "Other"]
unique(C1$taxonomy)

df2 <- rbind(A1, B1, C1)

#creating taxonomy table with replicates for excel-Genus level
tax <- df2
tax2 <- aggregate(Abundance~taxonomy + SamplingPosition + SampleNumbers + replicate_alt + Location + Date + TimePoint, tax, sum)
dw2 <- dcast(tax2, taxonomy ~ SampleNumbers + SamplingPosition + replicate_alt + Location + Date + TimePoint, value.var="Abundance")
#write.table(dw2, file = "replicates_Genus-natural-5pct", row.names = FALSE, sep = "\t", quote = FALSE, dec = ",")


plot_theme = theme(panel.grid.major = element_line(size = 0.5, color = "grey"), 
                   text = element_text(size = 18)) + 
  theme(axis.title.y = element_text(vjust=2)) +
  theme(panel.background = element_rect(fill = NA))

####################################################################################
####################################################################################
###plot


df2 <- aggregate(Abundance~taxonomy + SamplingPosition + SampleNumbers + replicate_alt + Location + Date + TimePoint, df2, FUN=sum)
#write.table(df2, file = "replicates_Genus-natural-10pct-long", row.names = FALSE, sep = "\t", quote = FALSE, dec = ",")
unique(df2$taxonomy)

#df2$TimePoint <- gsub('T1', 'T3', df2$TimePoint)
#df2$TimePoint <- gsub('T2', 'T4', df2$TimePoint)
df2$taxonomy <- gsub('Unassigned/Unassigned/Unassigned', 'Unassigned', df2$taxonomy)

df2 <- df2[order(df2$taxonomy),]
unique(df2$taxonomy)

cols <- c(`Acidimicrobiales/Acidimicrobiales Incertae Sedis/Candidatus Microthrix`="red",
          `Burkholderiales/Comamonadaceae/Sphaerotilus`="#745700",
          `Burkholderiales/Comamonadaceae/Unassigned`="#a77d00",
          `Burkholderiales/Oxalobacteraceae/Massilia`="#daa300",
          `Caldilineales/Caldilineaceae/uncultured`="#c13100",
          `Caulobacterales/Hyphomonadaceae/Woodsholea`="#007558",
          `Corynebacteriales/Unassigned/Unassigned`="#290a00",
          `Flavobacteriales/Flavobacteriaceae/Flavobacterium`="#314100",
          `Flavobacteriales/Flavobacteriaceae/Maribacter`="#577400",
          `Nitrospirales/Nitrospiraceae/Nitrospira`="#00c131",
          `Pseudomonadales/Moraxellaceae/Psychrobacter`="#0090c1",
          `Rhodobacterales/Rhodobacteraceae/Loktanella`="#8e006a",
          `Rhodobacterales/Rhodobacteraceae/Pseudorhodobacter`="#c10090",
          `Rhodobacterales/Rhodobacteraceae/Unassigned`="#f400b6",
          `Sphingomonadales/Sphingomonadaceae/Novosphingobium`="#ffc510",
          `Thiotrichales/Thiotrichaceae/Thiothrix`="#dca700",
          Other="gray41",
          Unassigned="gray90")

df2$taxonomy = factor(df2$taxonomy, levels=c("Acidimicrobiales/Acidimicrobiales Incertae Sedis/Candidatus Microthrix",
                                             "Burkholderiales/Comamonadaceae/Sphaerotilus",                           
                                             "Burkholderiales/Comamonadaceae/Unassigned",                             
                                             "Burkholderiales/Oxalobacteraceae/Massilia",                             
                                             "Caldilineales/Caldilineaceae/uncultured",                               
                                             "Caulobacterales/Hyphomonadaceae/Woodsholea",                            
                                             "Corynebacteriales/Unassigned/Unassigned",                               
                                             "Flavobacteriales/Flavobacteriaceae/Flavobacterium",                     
                                             "Flavobacteriales/Flavobacteriaceae/Maribacter",                         
                                             "Nitrospirales/Nitrospiraceae/Nitrospira",
                                             "Pseudomonadales/Moraxellaceae/Psychrobacter",                           
                                             "Rhodobacterales/Rhodobacteraceae/Loktanella",                           
                                             "Rhodobacterales/Rhodobacteraceae/Pseudorhodobacter",                    
                                             "Rhodobacterales/Rhodobacteraceae/Unassigned",                           
                                             "Sphingomonadales/Sphingomonadaceae/Novosphingobium",                    
                                             "Thiotrichales/Thiotrichaceae/Thiothrix",
                                             "Other",
                                             "Unassigned"))

df2$pos_rep <- paste0(df2$SamplingPosition, sep="-", df2$Replicate)
p5 <- ggplot(df2, aes(x=pos_rep, y=Abundance, fill=taxonomy, )) + 
  geom_bar(stat = "identity", color="black") +
  scale_fill_manual(values = cols) +
  theme(legend.position="right") +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("Relative abundance (% of total sequences)") +
  xlab("Replicate") +
  facet_grid(Location~TimePoint, scale="free_x") +
  ggtitle("") + 
  science_theme2
p5