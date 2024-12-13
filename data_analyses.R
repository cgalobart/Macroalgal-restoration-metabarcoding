####READ ME####
#FILE NAME: data_analyses.R     
#DATE: last update 25/11/2024
#TITLE: Metabarcoding identifies macroalgal composition as a driver of benthic invertebrate assemblages in restored habitats
#AUTHORS: Cristina Galobart, Jesús Zarcero, Adrià Antich, Xavier Turon, Emma Cebrian
#SCRIPT: C Galobart (cgalobart@ceab.csic.es; cgalobart@gmail.com) with support from J Zarcero and A Antich
#JOURNAL: Scientific Reports - under revision

# DISCLAIMER: This script has been developed by an ecologist, not a professional programmer, 
# so please take into account that there may be room for code optimization. 
# Positive feedback is always welcome and appreciated!

# Script Content
# 0. Plot theme and colours
# 1. Figure 1 - Biomass main macrophytes
# 2. Load data and data arrangements
# 3. Figure 2 - macroinvertebrate species diversity and evenness + ANOVAs
# 4. Figure 3 - non-metric multidimensional scaling (nMDS) + PERMANOVA
# 5. Load data and data arrangements
# 6. Figure 4 upper - Frequency of read abundance
# 7. Figure 4 lower - Frequency of MOTU abundance
# 8. Load data and data arrangements
# 9. Table 1 - MVABUND analyses
# 10. Figure 5 - IndVal analyses (indicator species)
# 11. Supplementary information, Figure S2 - Rarefaction and accumulation curves

#Load libraries
library(xlsx)
library(car)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(vegan)
library(mvabund)
library(indicspecies)
library(labdsv)
library(reshape2)
library(BiodiversityR)

#Set working directory 
setwd("C:/")
options(scipen=999)

####0. Plot theme and colours####
#Set the plot theme and colours
theme_set(theme_light() + 
            theme(plot.title = element_blank(),
                  strip.text = element_text(color = "black", size = 18),
                  axis.text.y = element_text(color = "black", size = 12),
                  axis.text.x = element_text(color = "black", size = 12, angle = 90)))

pal <- c("#DD542C", "#C08E70", "#98C190", "#777539")

#####1. Figure 1 - Biomass main macrophytes####
biomass <- read.csv('biomass_macrophytes.csv', sep = ';', stringsAsFactors = T)

#Summary of biomass data (mean and sd)
biomass_summary <- biomass %>% 
  group_by(assemblage) %>% 
  summarise(Cystoseira.sl_mean = mean(Cystoseira.sp),
            Cystoseira.sl_sd = sd(Cystoseira.sp, na.rm = TRUE),
            C.nodosa_mean = mean(Cymodocea.nodosa),
            C.nodosa_sd = sd(Cymodocea.nodosa, na.rm = TRUE),
            P.pavonica_mean = mean(Padina.pavonica),
            P.pavonica_sd = sd(Padina.pavonica, na.rm = TRUE),
            H.scoparia_mean = mean(Halopteris.scoparia),
            H.scoparia_sd = sd(Halopteris.scoparia, na.rm = TRUE))

#Prepare the dataframe for the plot
b_sum_long <- biomass_summary %>%
  pivot_longer(
    cols = -assemblage,         # Columns to pivot (exclude "assemblage")
    names_to = c("species", "type"),  # Split into "species" and "type"
    names_sep = "_"             # Separator to differentiate "species" and "type"
  ) %>%
  pivot_wider(
    names_from = type,          # Spread "type" into separate columns
    values_from = value         # Values come from the "value" column
  )

#Reorder the levels
b_sum_long$assemblage <- factor(b_sum_long$assemblage, levels = c("no_restored", "restored_forest", "cala_rotja_forest", "miami_forest"))
b_sum_long$species <- factor(b_sum_long$species, levels = c("Cystoseira.sl", "C.nodosa", "P.pavonica", "H.scoparia"))

#Create the plots
ggplot(b_sum_long, aes(x = species, y = mean, fill = assemblage)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1, position = position_dodge(0.9)) + # Add error bars
  labs(y = "biomass (g of wet weight)") +
  scale_fill_manual(values = pal) +
  ylim (0, 100) + 
  theme(legend.position = "none") +
  facet_wrap(~assemblage, nrow = 1)

ggsave("Figure_1.pdf", width=40, height=15, units= "cm")

####2. Load data and data arrangements####
metadata <- read.csv('metadata.csv', sep = ';', stringsAsFactors = T)
abundance <- read.csv('MOTUs_data.csv', sep = ';', stringsAsFactors = T)
str(abundance)

#Remove the "non identified MOTUs"
abundance <- abundance %>% 
  filter(!phylum_name == "Noid")

#Match metadata samples with the abundance table, in case there are more samples in the metadata
metadata <- metadata[metadata$codi %in% colnames(abundance),]
metadata$codi

#Split the metadata columns in "sample columns-> abundance" and "other information/taxonomic -> taxo"
abund <- abundance[,colnames(abundance) %in% metadata$codi]
#MOTU id as rownames
row.names(abund) <- abundance$id
taxo <- abundance[,!colnames(abundance) %in% metadata$codi]

#Read abundance transformed to relative read abundance (per sample)
#Prepare function
abund_rel <- function(x){
  if (sum(x)>0) {
    x <- x/sum(x)
  } else{
    x<-x
  }
  return(x)
}

std <- function(x) sd(x)/sqrt(length(x))

#Absolute reads to relative
abund_relative <- apply(abund, 2, abund_rel)
#To check
colSums(abund_relative)

####3. Figure 2 - macroinvertebrate species diversity and evenness + ANOVAs####

#Shannon index diversity
shan <- diversity(t(abund_relative))

#Pielou index evenness
#We first calculate species richness
sp_rich <- specnumber(t(abund_relative))
#Now the eveness index
even <- shan/log(sp_rich)

#Create a dataframe to store the outcomes (shan, even, and others later)
indices <- as.data.frame(add_rownames(as.data.frame(shan)))
colnames(indices)[colnames(indices) == 'rowname'] <- "sample_name"
#Add the outcomes
indices <- cbind(indices, even)

#Prepare data to create the plots
metadata_with_indices <- metadata
colnames(metadata_with_indices)[colnames(metadata_with_indices) %in% 'codi'] <- 'sample_name'
met_index <- left_join(metadata_with_indices, indices, by = 'sample_name')

#Reorder levels of "assemblage" factor
met_index$assemblage <- factor(met_index$assemblage, levels = c("Teulera_no_restaurat", "Teulera_restaurat", "Fornells_Cala_Roja", "Fornells_Miami"))

#Shannon index plot
ggplot(met_index, aes(x = assemblage, y = shan, fill = assemblage)) +
  geom_boxplot(size = 0.5) +
  facet_wrap(~fraction) + 
  scale_fill_manual(values = pal) + 
  theme (axis.text.x = element_text(angle = 45, vjust = 0.5),
         legend.position = "none")

ggsave("Figure_2_upper.pdf", width=20, height=15, units= "cm")

#One-way ANOVA for Shannon index
mod <- aov(shan ~ assemblage * fraction, data = met_index)
summary(mod)

#Anova assumptions
#Normality and Homocedasticity
qqPlot(residuals(mod))
plot(fitted(mod), residuals(mod), xlab="Fitted Values", ylab="Residuals")
abline(h=0, lty=2)
lines(lowess(fitted(mod), residuals(mod)))

shapiro.test(mod$residuals) # formal normality statistical test
leveneTest(shan ~ assemblage, data = met_index) #formal statistical test for homogeneity of variance
leveneTest(shan ~ fraction, data = met_index) #formal statistical test for homogeneity of variance

#Tukey test for pair-wise comparison
TukeyHSD(mod)

#Pielou index plot
ggplot(met_index, aes(x = assemblage, y = even, fill = assemblage)) +
  geom_boxplot(size = 0.5) +
  facet_wrap(~fraction) + 
  scale_fill_manual(values = pal) + 
  theme (axis.text.x = element_text(angle = 45, vjust = 0.5),
         legend.position = "none")

ggsave("Figure_2_lower.pdf", width=20, height=15, units= "cm")

#One-way ANOVA for Pielou index
mod <- aov(even ~ assemblage * fraction, data = met_index)
summary(mod)

#Anova assumptions
#Normality and Homocedasticity
qqPlot(residuals(mod))
plot(fitted(mod), residuals(mod), xlab="Fitted Values", ylab="Residuals")
abline(h=0, lty=2)
lines(lowess(fitted(mod), residuals(mod)))

shapiro.test(mod$residuals) # formal normality statistical test
leveneTest(even ~ assemblage, data = met_index) #formal statistical test for homogeneity of variance
leveneTest(even ~ fraction, data = met_index)#formal statistical test for homogeneity of variance

#Tukey test for pair-wise comparison
TukeyHSD(mod)

####4. Figure 3 - non-metric multidimensional scaling (nMDS) + PERMANOVA#### 

#fourth root of abundance
#prepare the function
arrel_quarta <- function(x){
  x <- sqrt(sqrt(x))
  return(x)
}

abund_relative_q <- arrel_quarta(t(abund_relative))
abund_relative_q <- as.data.frame(abund_relative_q)

#Calculate the distance matrix, Bray-Curtis dissimilarity
dist_matrix <- vegdist(abund_relative_q, method = "bray")

#Perform the nMDS, 2 dimensions (k=2)
mds_data <- metaMDS(dist_matrix, trymax = 600, k = 2)
#Check stress value
mds_data$stress
#Store the NMDS point coordenates
mds_points <- as.data.frame(add_rownames(as.data.frame(mds_data$points)))
colnames(mds_points)[colnames(mds_points) == 'rowname'] <- "sample_name"

#Prepare the dataframe to do the plot
metadata_nmds <- metadata
metadata_nmds$assemblage <- factor(metadata_nmds$assemblage, levels = c("Teulera_no_restaurat", "Teulera_restaurat", "Fornells_Cala_Roja", "Fornells_Miami"))
colnames(metadata_nmds)[colnames(metadata_nmds) %in% 'codi'] <- 'sample_name'
#Join both dataframes with necessary information
metadata_nmds <- left_join(metadata_nmds, mds_points, by = 'sample_name')

#NMDS plot
#First calculate the centoids and segments 
#AB fraction
metadata_nmds_AB <- metadata_nmds %>% 
  filter(fraction == "AB")
cent_AB <- aggregate(cbind(MDS1, MDS2) ~ assemblage, data = metadata_nmds_AB, FUN = mean)
segs_AB <- merge(metadata_nmds_AB, setNames(cent_AB, c('assemblage','MDS1_c','MDS2_c')),
                 by = 'assemblage', sort = FALSE)
#C fraction
metadata_nmds_C <- metadata_nmds %>% 
  filter(fraction == "C")
cent_C <- aggregate(cbind(MDS1, MDS2) ~ assemblage, data = metadata_nmds_C, FUN = mean)
segs_C <- merge(metadata_nmds_C, setNames(cent_C, c('assemblage','MDS1_c','MDS2_c')),
                by = 'assemblage', sort = FALSE)

#Do the plot
ggplot(metadata_nmds, aes(x = MDS1, y = MDS2, colour = assemblage, shape = fraction)) +
  geom_segment(data = segs_AB, mapping = aes(xend = MDS1_c, yend = MDS2_c)) + 
  geom_segment(data = segs_C, mapping = aes(xend = MDS1_c, yend = MDS2_c)) +  
  geom_point(data = cent_AB, size = 12, shape = 16) +   
  geom_point(data = cent_C, size = 12, shape = 17) +   
  geom_point(size = 3) +                               
  scale_color_manual(values = pal) +
  scale_shape_manual(values = c(16,17)) +
  coord_fixed() +
  theme_light() 

ggsave("Figure_3.pdf", width=25, height=15, units= "cm")

#Test differences with PERMANOVA
permanova <- adonis2(abund_relative_q ~ assemblage * fraction, data=metadata_nmds, permutations=999, method ="bray")
permanova

#Pair-wise test (assemblage)
#Activate the function - https://github.com/pmartinezarbizu/pairwiseAdonis
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni',reduce=NULL,perm=999)
{
  
  co <- combn(unique(as.character(factors)),2)
  pairs <- c()
  Df <- c()
  SumsOfSqs <- c()
  F.Model <- c()
  R2 <- c()
  p.value <- c()
  
  
  for(elem in 1:ncol(co)){
    if(inherits(x, 'dist')){
      x1=as.matrix(x)[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),
                      factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))]
    }
    
    else  (
      if (sim.function == 'daisy'){
        x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
      }
      else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    )
    
    x2 = data.frame(Fac = factors[factors %in% c(co[1,elem],co[2,elem])])
    
    ad <- adonis2(x1 ~ Fac, data = x2,
                  permutations = perm);
    pairs <- c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    Df <- c(Df,ad$Df[1])
    SumsOfSqs <- c(SumsOfSqs,ad$SumOfSqs[1])
    F.Model <- c(F.Model,ad$F[1]);
    R2 <- c(R2,ad$R2[1]);
    p.value <- c(p.value,ad$`Pr(>F)`[1])
  }
  p.adjusted <- p.adjust(p.value,method=p.adjust.m)
  
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  pairw.res <- data.frame(pairs,Df,SumsOfSqs,F.Model,R2,p.value,p.adjusted,sig)
  
  if(!is.null(reduce)){
    pairw.res <- subset (pairw.res, grepl(reduce,pairs))
    pairw.res$p.adjusted <- p.adjust(pairw.res$p.value,method=p.adjust.m)
    
    sig = c(rep('',length(pairw.res$p.adjusted)))
    sig[pairw.res$p.adjusted <= 0.1] <-'.'
    sig[pairw.res$p.adjusted <= 0.05] <-'*'
    sig[pairw.res$p.adjusted <= 0.01] <-'**'
    sig[pairw.res$p.adjusted <= 0.001] <-'***'
    pairw.res <- data.frame(pairw.res[,1:7],sig)
  }
  class(pairw.res) <- c("pwadonis", "data.frame")
  return(pairw.res)
}
#Method summary
summary.pwadstrata = function(object, ...) {
  cat("Result of pairwise.adonis2:\n")
  cat("\n")
  print(object[1], ...)
  cat("\n")
  
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
}

#Perform the pair-wise test
assemblage <- pairwise.adonis(x = abund_relative_q, factors = metadata_nmds$assemblage, 
                sim.function = "vegdist", sim.method = "bray", p.adjust.m = "bonferroni")
assemblage


#Remove all objects in the Environment tab
####5. Load data and data arrangements####
metadata <- read.csv('metadata.csv', sep = ";", stringsAsFactors = T)
metadata$assemblage <- factor(metadata$assemblage, levels = c("Teulera_no_restaurat", "Teulera_restaurat", "Fornells_Cala_Roja", "Fornells_Miami"))

abundance <- read.csv('MOTUs_data.csv', sep = ';', stringsAsFactors = T)

#Delete the "non identified MOTUs"
abundance <- abundance %>% 
  filter(!phylum_name == "Noid")

#Match metadata samples with the abundance table, in case there are more samples in the metadata
metadata <- metadata[metadata$codi %in% colnames(abundance),]
metadata$codi

#Split the metadata columns in "sample columns-> abundance" and "other information/taxonomic -> taxo"
abund <- abundance[,colnames(abundance) %in% metadata$codi]
#MOTU id as rownames
row.names(abund) <- abundance$id
taxo <- abundance[,!colnames(abundance) %in% metadata$codi]

####6. Figure 4 upper - Frequency of read abundance####
#Identify MOTUS with higher number of reads
reads_metazoa <- abundance %>% 
  group_by(phylum_name) %>% 
  summarise(reads = sum(total_reads))
reads_metazoa #Arthropoda, Mollusca, and Annelida

#Total number of reads
total_sum_all <- sum(reads_metazoa$reads) 

#Calculate percentage of reads for each phylum
reads_metazoa <- reads_metazoa %>%
  mutate(percentage_reads = (reads / total_sum_all) * 100)
reads_metazoa

#Keep only the information we want to calculate read proportions
metazoa <- abundance %>% 
  select(8, 14:37) %>%  #select phylum_name and all samples
  group_by(phylum_name) %>% 
  summarise(across(everything(), sum))

#Total number of reads per sample
total_reads <- metazoa %>%
  summarise(across(-phylum_name, sum))
total_reads

#Calculate percentage of reads for each phylum in each sample
percentage_reads <- metazoa %>%
  mutate(across(-phylum_name, ~ . / total_reads[[cur_column()]] * 100))

#Pivot table
pivot_table <- percentage_reads %>%
  pivot_longer(cols = -phylum_name, names_to = "codi", values_to = "percentage") %>%
  pivot_wider(names_from = phylum_name, values_from = percentage)

#Add assemblage & fraction information per sample
percent <- merge(pivot_table, metadata, by = "codi") %>% 
  select(2:14, assemblage, fraction) 

#Mean percentage of reads for phylum of the three replicates
percent_plot <- percent %>% 
  group_by(assemblage, fraction) %>% 
  summarise(across(everything(), mean))

#Get the data ready for the plot
all_info_reads <- gather(percent_plot, key = "phylum_name", value = "reads", -assemblage, -fraction)

#Select colors for each phylum (those phylums with less than 5% of reads are coloured all white)
palette <- c("#f17890", "#aa3355", "#722062", "#3366bb", "#0099cc", "#44dd88", "#e3dbdb",
             "#e3dbdb", "#e3dbdb", "#e3dbdb", "#e3dbdb", "#e3dbdb", "#e3dbdb")

all_info_reads$phylum_name <- factor(all_info_reads$phylum_name, 
                                     levels = c("Arthropoda", "Mollusca", "Annelida", "Porifera",  
                                                "Chordata", "Cnidaria", "Echinodermata", "Nemertea","Nematoda", "Bryozoa", "Rotifera", 
                                                "Platyhelminthes", "Kinorhyncha"))

#Do the plot
ggplot(all_info_reads, aes(x = assemblage, y = reads, fill = phylum_name)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = palette) +
  facet_wrap(~fraction)

ggsave("Figure_4_upper.pdf", width=22, height=17, units= "cm")


####7. Figure 4 lower - Frequency of MOTU abundance####
#Keep only the information we want to calculate MOTU proportions
metazoa <- abundance %>% 
  select(8, 14:37) 

#Calculate number of MOTUs for each phylum in each sample
met_summ <- metazoa %>% 
  group_by(phylum_name) %>% 
  summarise(across(everything(), ~ sum(. > 0)))
met_summ

#Total MOTUs per sample
total_motus <- met_summ %>%
  summarise(across(-phylum_name, sum))
total_motus

#Calculate percentage of MOTUs for each phylum in each sample
percentage_motus <- met_summ %>%
  mutate(across(-phylum_name, ~ . / total_reads[[cur_column()]] * 100))

#Pivot table
pivot_motus <- percentage_motus %>%
  pivot_longer(cols = -phylum_name, names_to = "codi", values_to = "percentage") %>%
  pivot_wider(names_from = phylum_name, values_from = percentage)

#Add assemblage & fraction information per sample
percent <- merge(pivot_motus, metadata, by = "codi") %>% 
  select(2:14, assemblage, fraction)

#Mean percentage of MOTUs for phylum of the three replicates
percent_plot <- percent %>% 
  group_by(assemblage, fraction) %>% 
  summarise(across(everything(), mean))

#Get the data ready for the plot
all_info_reads <- gather(percent_plot, key = "phylum_name", value = "reads", -assemblage, -fraction)

#Select colors for each phylum (those phylums with less than 5% of reads are coloured all white)
palette <- c("#f17890", "#aa3355", "#722062", "#3366bb", "#0099cc", "#44dd88", "#e3dbdb",
             "#e3dbdb", "#e3dbdb", "#e3dbdb", "#e3dbdb", "#e3dbdb", "#e3dbdb")

all_info_reads$phylum_name <- factor(all_info_reads$phylum_name, 
                                     levels = c("Arthropoda", "Mollusca", "Annelida", "Porifera",  
                                                "Chordata", "Cnidaria", "Echinodermata", "Nemertea","Nematoda", "Bryozoa", "Rotifera", 
                                                "Platyhelminthes", "Kinorhyncha"))

#Do the plot
ggplot(all_info_reads, aes(x = assemblage, y = reads, fill = phylum_name)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = palette) +
  facet_wrap(~fraction)

ggsave("Figure_4_lower.pdf", width=22, height=17, units= "cm")


#Remove all objects in the Environment tab
####8. Load data and data arrangements####
metadata <- read.csv('metadata.csv', sep = ";", stringsAsFactors = T)
metadata$assemblage <- factor(metadata$assemblage, levels = c("Teulera_no_restaurat", "Teulera_restaurat", "Fornells_Cala_Roja", "Fornells_Miami"))

#For mvabund and indval analyses (see methods) MOTUs with relative abundances >0.05% for each taxonomic group were kept
abundance_or <- read.csv('MOTUs_data.csv', sep = ';', stringsAsFactors = T)
levels(abundance_or$phylum_name)

#Remove those MOTUs that were not identified at phylum level
abundance_or <- abundance_or %>% 
  filter(!phylum_name == "Noid")

#Obtain the MOTUs with >0.05% abundance for each phylum
#ALL MOTUs TOGETHER
abundance <- abundance_or
metadata <- metadata[metadata$codi %in% colnames(abundance),]
abund <- abundance[,colnames(abundance) %in% metadata$codi]
row.names(abund) <- abundance$id
taxo <- abundance[,!colnames(abundance) %in% metadata$codi]
#Step 1: Calculate total reads of each MOTU across all samples
MOTUs_ab <- rowSums(abund)
#Step 2: Calculate the total reads across all MOTUs and samples
reads_all <- sum(abund)
#Step 3: 0.05% of total reads number
reads_005 <- reads_all*0.05/100
#Step 4: Keep those MOTUs that total reads is higher than the 0.05% value
abun_filt <- abund %>%
  filter(rowSums(abund) >= reads_005)

#ARTHROPODA
abundance <- abundance_or %>% 
  filter(phylum_name == "Arthropoda")
metadata <- metadata[metadata$codi %in% colnames(abundance),]
abund <- abundance[,colnames(abundance) %in% metadata$codi]
row.names(abund) <- abundance$id
taxo <- abundance[,!colnames(abundance) %in% metadata$codi]
#Step 1: Calculate total reads of each MOTU across all samples
MOTUs_ab <- rowSums(abund)
#Step 2: Calculate the total reads across all MOTUs and samples
reads_all <- sum(abund)
#Step 3: 0.05% of total reads number
reads_005 <- reads_all*0.05/100
#Step 4: Keep those MOTUs that total reads is higher than the 0.05% value
abun_filt_art <- abund %>%
  filter(rowSums(abund) >= reads_005)

#MOLLUSCA
abundance <- abundance_or %>% 
  filter(phylum_name == "Mollusca")
metadata <- metadata[metadata$codi %in% colnames(abundance),]
abund <- abundance[,colnames(abundance) %in% metadata$codi]
row.names(abund) <- abundance$id
taxo <- abundance[,!colnames(abundance) %in% metadata$codi]
#Step 1: Calculate total reads of each MOTU across all samples
MOTUs_ab <- rowSums(abund)
#Step 2: Calculate the total reads across all MOTUs and samples
reads_all <- sum(abund)
#Step 3: 0.05% of total reads number
reads_005 <- reads_all*0.05/100
#Step 4: Keep those MOTUs that total reads is higher than the 0.05% value
abun_filt_moll <- abund %>%
  filter(rowSums(abund) >= reads_005)

#ANNELIDA
abundance <- abundance_or %>% 
  filter(phylum_name == "Annelida")
metadata <- metadata[metadata$codi %in% colnames(abundance),]
abund <- abundance[,colnames(abundance) %in% metadata$codi]
row.names(abund) <- abundance$id
taxo <- abundance[,!colnames(abundance) %in% metadata$codi]
#Step 1: Calculate total reads of each MOTU across all samples
MOTUs_ab <- rowSums(abund)
#Step 2: Calculate the total reads across all MOTUs and samples
reads_all <- sum(abund)
#Step 3: 0.05% of total reads number
reads_005 <- reads_all*0.05/100
#Step 4: Keep those MOTUs that total reads is higher than the 0.05% value
abun_filt_ann <- abund %>%
  filter(rowSums(abund) >= reads_005)

#PORIFERA
abundance <- abundance_or %>% 
  filter(phylum_name == "Porifera")
metadata <- metadata[metadata$codi %in% colnames(abundance),]
abund <- abundance[,colnames(abundance) %in% metadata$codi]
row.names(abund) <- abundance$id
taxo <- abundance[,!colnames(abundance) %in% metadata$codi]
#Step 1: Calculate total reads of each MOTU across all samples
MOTUs_ab <- rowSums(abund)
#Step 2: Calculate the total reads across all MOTUs and samples
reads_all <- sum(abund)
#Step 3: 0.05% of total reads number
reads_005 <- reads_all*0.05/100
#Step 4: Keep those MOTUs that total reads is higher than the 0.05% value
abun_filt_por <- abund %>%
  filter(rowSums(abund) >= reads_005)

#CHORDATA
abundance <- abundance_or %>% 
  filter(phylum_name == "Chordata")
metadata <- metadata[metadata$codi %in% colnames(abundance),]
abund <- abundance[,colnames(abundance) %in% metadata$codi]
row.names(abund) <- abundance$id
taxo <- abundance[,!colnames(abundance) %in% metadata$codi]
#Step 1: Calculate total reads of each MOTU across all samples
MOTUs_ab <- rowSums(abund)
#Step 2: Calculate the total reads across all MOTUs and samples
reads_all <- sum(abund)
#Step 3: 0.05% of total reads number
reads_005 <- reads_all*0.05/100
#Step 4: Keep those MOTUs that total reads is higher than the 0.05% value
abun_filt_cho <- abund %>%
  filter(rowSums(abund) >= reads_005)

#CNIDARIA
abundance <- abundance_or %>% 
  filter(phylum_name == "Cnidaria")
metadata <- metadata[metadata$codi %in% colnames(abundance),]
abund <- abundance[,colnames(abundance) %in% metadata$codi]
row.names(abund) <- abundance$id
taxo <- abundance[,!colnames(abundance) %in% metadata$codi]
#Step 1: Calculate total reads of each MOTU across all samples
MOTUs_ab <- rowSums(abund)
#Step 2: Calculate the total reads across all MOTUs and samples
reads_all <- sum(abund)
#Step 3: 0.05% of total reads number
reads_005 <- reads_all*0.05/100
#Step 4: Keep those MOTUs that total reads is higher than the 0.05% value
abun_filt_cni <- abund %>%
  filter(rowSums(abund) >= reads_005)

#BRYOZOA
abundance <- abundance_or %>% 
  filter(phylum_name == "Bryozoa")
metadata <- metadata[metadata$codi %in% colnames(abundance),]
abund <- abundance[,colnames(abundance) %in% metadata$codi]
row.names(abund) <- abundance$id
taxo <- abundance[,!colnames(abundance) %in% metadata$codi]
#Step 1: Calculate total reads of each MOTU across all samples
MOTUs_ab <- rowSums(abund)
#Step 2: Calculate the total reads across all MOTUs and samples
reads_all <- sum(abund)
#Step 3: 0.05% of total reads number
reads_005 <- reads_all*0.05/100
#Step 4: Keep those MOTUs that total reads is higher than the 0.05% value
abun_filt_bry <- abund %>%
  filter(rowSums(abund) >= reads_005)

#ECHINODERMATA
abundance <- abundance_or %>% 
  filter(phylum_name == "Echinodermata")
metadata <- metadata[metadata$codi %in% colnames(abundance),]
abund <- abundance[,colnames(abundance) %in% metadata$codi]
row.names(abund) <- abundance$id
taxo <- abundance[,!colnames(abundance) %in% metadata$codi]
#Step 1: Calculate total reads of each MOTU across all samples
MOTUs_ab <- rowSums(abund)
#Step 2: Calculate the total reads across all MOTUs and samples
reads_all <- sum(abund)
#Step 3: 0.05% of total reads number
reads_005 <- reads_all*0.05/100
#Step 4: Keep those MOTUs that total reads is higher than the 0.05% value
abun_filt_ech <- abund %>%
  filter(rowSums(abund) >= reads_005)

#KINORHYNCHA
abundance <- abundance_or %>% 
  filter(phylum_name == "Kinorhyncha")
metadata <- metadata[metadata$codi %in% colnames(abundance),]
abund <- abundance[,colnames(abundance) %in% metadata$codi]
row.names(abund) <- abundance$id
taxo <- abundance[,!colnames(abundance) %in% metadata$codi]
#Step 1: Calculate total reads of each MOTU across all samples
MOTUs_ab <- rowSums(abund)
#Step 2: Calculate the total reads across all MOTUs and samples
reads_all <- sum(abund)
#Step 3: 0.05% of total reads number
reads_005 <- reads_all*0.05/100
#Step 4: Keep those MOTUs that total reads is higher than the 0.05% value
abun_filt_kin <- abund %>%
  filter(rowSums(abund) >= reads_005)

#NEMATODA
abundance <- abundance_or %>% 
  filter(phylum_name == "Nematoda")
metadata <- metadata[metadata$codi %in% colnames(abundance),]
abund <- abundance[,colnames(abundance) %in% metadata$codi]
row.names(abund) <- abundance$id
taxo <- abundance[,!colnames(abundance) %in% metadata$codi]
#Step 1: Calculate total reads of each MOTU across all samples
MOTUs_ab <- rowSums(abund)
#Step 2: Calculate the total reads across all MOTUs and samples
reads_all <- sum(abund)
#Step 3: 0.05% of total reads number
reads_005 <- reads_all*0.05/100
#Step 4: Keep those MOTUs that total reads is higher than the 0.05% value
abun_filt_nem <- abund %>%
  filter(rowSums(abund) >= reads_005)

#NEMERTEA
abundance <- abundance_or %>% 
  filter(phylum_name == "Nemertea")
metadata <- metadata[metadata$codi %in% colnames(abundance),]
abund <- abundance[,colnames(abundance) %in% metadata$codi]
row.names(abund) <- abundance$id
taxo <- abundance[,!colnames(abundance) %in% metadata$codi]
#Step 1: Calculate total reads of each MOTU across all samples
MOTUs_ab <- rowSums(abund)
#Step 2: Calculate the total reads across all MOTUs and samples
reads_all <- sum(abund)
#Step 3: 0.05% of total reads number
reads_005 <- reads_all*0.05/100
#Step 4: Keep those MOTUs that total reads is higher than the 0.05% value
abun_filt_neme <- abund %>%
  filter(rowSums(abund) >= reads_005)

#PLATHYELMINTHES
abundance <- abundance_or %>% 
  filter(phylum_name == "Plathyelminthes")
metadata <- metadata[metadata$codi %in% colnames(abundance),]
abund <- abundance[,colnames(abundance) %in% metadata$codi]
row.names(abund) <- abundance$id
taxo <- abundance[,!colnames(abundance) %in% metadata$codi]
#Step 1: Calculate total reads of each MOTU across all samples
MOTUs_ab <- rowSums(abund)
#Step 2: Calculate the total reads across all MOTUs and samples
reads_all <- sum(abund)
#Step 3: 0.05% of total reads number
reads_005 <- reads_all*0.05/100
#Step 4: Keep those MOTUs that total reads is higher than the 0.05% value
abun_filt_pla <- abund %>%
  filter(rowSums(abund) >= reads_005)

#ROTIFERA
abundance <- abundance_or %>% 
  filter(phylum_name == "Rotifera")
metadata <- metadata[metadata$codi %in% colnames(abundance),]
abund <- abundance[,colnames(abundance) %in% metadata$codi]
row.names(abund) <- abundance$id
taxo <- abundance[,!colnames(abundance) %in% metadata$codi]
#Step 1: Calculate total reads of each MOTU across all samples
MOTUs_ab <- rowSums(abund)
#Step 2: Calculate the total reads across all MOTUs and samples
reads_all <- sum(abund)
#Step 3: 0.05% of total reads number
reads_005 <- reads_all*0.05/100
#Step 4: Keep those MOTUs that total reads is higher than the 0.05% value
abun_filt_rot <- abund %>%
  filter(rowSums(abund) >= reads_005)

####9. Table 1 - MVABUND analyses####
#obtained from https://environmentalcomputing.net/statistics/mvabund/ and references therein

#Load biomass data, as categories (0, 1-25, 26-50, 51-75, >75 gr of wet weight per quadrat)
biomass <- read.csv('biomass_macrophytes_categorical.csv', sep = ';', stringsAsFactors = T)
rownames(biomass) <- biomass$X 
biomass$X <- NULL

#Relative function
abund_rel <- function(x){
  if (sum(x)>0) {
    x <- x*100/sum(x)
  } else{
    x<-x
  }
  return(x)
}

#Prepare database for each mvabund analyses, one for all MOTUs together and one for each phylum
#ALL MOTUs TOGETHER
#Compute relative read abundance
motus_abund_rel <- as.data.frame(t(apply(abun_filt, 2, abund_rel)))
rowSums(motus_abund_rel)
#Round relative abundances to "integer numbers" (has to be this way for the mvabund analyses)
abund_rel_c <- ceiling(motus_abund_rel)
#Convert dataframe to an mvabund object format
abund_mvabund <- mvabund(abund_rel_c)
#Perform the model
glm.data <- manyglm(abund_mvabund ~ Cystoseira.sp + Cymodocea.nodosa + Padina.pavonica + Halopteris.scoparia, family = "negative.binomial", data = biomass)
#Assumptions to check
#Plot of the residuals (we should see a random scatter of points, not linear or curvilinear relationship, or a fan shape)
plot(glm.data)
#Statistical test of the fitted model and each variable
glm_result <- anova(glm.data)
glm_result

#ARTHROPODA
#Compute relative read abundance
motus_abund_rel <- as.data.frame(t(apply(abun_filt_art, 2, abund_rel)))
rowSums(motus_abund_rel)
#Round relative abundances to "integer numbers" (has to be this way for the mvabund analyses)
abund_rel_c <- ceiling(motus_abund_rel)
#Convert dataframe to an mvabund object format
abund_mvabund <- mvabund(abund_rel_c)
#Perform the model
glm.data <- manyglm(abund_mvabund ~ Cystoseira.sp + Cymodocea.nodosa + Padina.pavonica + Halopteris.scoparia, family = "negative.binomial", data = biomass)
#Assumptions to check
#Plot of the residuals (we should see a random scatter of points, not linear or curvilinear relationship, or a fan shape)
plot(glm.data)
#Statistical test of the fitted model and each variable
glm_result <- anova(glm.data)
glm_result

#MOLLUSCA
#Compute relative read abundance
motus_abund_rel <- as.data.frame(t(apply(abun_filt_moll, 2, abund_rel)))
rowSums(motus_abund_rel)
#Round relative abundances to "integer numbers" (has to be this way for the mvabund analyses)
abund_rel_c <- ceiling(motus_abund_rel)
#Convert dataframe to an mvabund object format
abund_mvabund <- mvabund(abund_rel_c)
#Perform the model
glm.data <- manyglm(abund_mvabund ~ Cystoseira.sp + Cymodocea.nodosa + Padina.pavonica + Halopteris.scoparia, family = "negative.binomial", data = biomass)
#Assumptions to check
#Plot of the residuals (we should see a random scatter of points, not linear or curvilinear relationship, or a fan shape)
plot(glm.data)
#Statistical test of the fitted model and each variable
glm_result <- anova(glm.data)
glm_result

#ANNELIDA
#Compute relative read abundance
motus_abund_rel <- as.data.frame(t(apply(abun_filt_ann, 2, abund_rel)))
rowSums(motus_abund_rel)
#Round relative abundances to "integer numbers" (has to be this way for the mvabund analyses)
abund_rel_c <- ceiling(motus_abund_rel)
#Convert dataframe to an mvabund object format
abund_mvabund <- mvabund(abund_rel_c)
#Perform the model
glm.data <- manyglm(abund_mvabund ~ Cystoseira.sp + Cymodocea.nodosa + Padina.pavonica + Halopteris.scoparia, family = "negative.binomial", data = biomass)
#Assumptions to check
#Plot of the residuals (we should see a random scatter of points, not linear or curvilinear relationship, or a fan shape)
plot(glm.data)
#Statistical test of the fitted model and each variable
glm_result <- anova(glm.data)
glm_result

#PORIFERA
#Compute relative read abundance
motus_abund_rel <- as.data.frame(t(apply(abun_filt_por, 2, abund_rel)))
rowSums(motus_abund_rel)
#Round relative abundances to "integer numbers" (has to be this way for the mvabund analyses)
abund_rel_c <- ceiling(motus_abund_rel)
#Convert dataframe to an mvabund object format
abund_mvabund <- mvabund(abund_rel_c)
#Perform the model
glm.data <- manyglm(abund_mvabund ~ Cystoseira.sp + Cymodocea.nodosa + Padina.pavonica + Halopteris.scoparia, family = "negative.binomial", data = biomass)
#Assumptions to check
#Plot of the residuals (we should see a random scatter of points, not linear or curvilinear relationship, or a fan shape)
plot(glm.data)
#Statistical test of the fitted model and each variable
glm_result <- anova(glm.data)
glm_result

#CHORDATA
#Compute relative read abundance
motus_abund_rel <- as.data.frame(t(apply(abun_filt_cho, 2, abund_rel)))
rowSums(motus_abund_rel)
#Round relative abundances to "integer numbers" (has to be this way for the mvabund analyses)
abund_rel_c <- ceiling(motus_abund_rel)
#Convert dataframe to an mvabund object format
abund_mvabund <- mvabund(abund_rel_c)
#Perform the model
glm.data <- manyglm(abund_mvabund ~ Cystoseira.sp + Cymodocea.nodosa + Padina.pavonica + Halopteris.scoparia, family = "negative.binomial", data = biomass)
#Assumptions to check
#Plot of the residuals (we should see a random scatter of points, not linear or curvilinear relationship, or a fan shape)
plot(glm.data)
#Statistical test of the fitted model and each variable
glm_result <- anova(glm.data)
glm_result

#CNIDARIA
#Compute relative read abundance
motus_abund_rel <- as.data.frame(t(apply(abun_filt_cni, 2, abund_rel)))
rowSums(motus_abund_rel)
#Round reltive abundances to "integer numbers" (has to be this way for the mvabund analyses)
abund_rel_c <- ceiling(motus_abund_rel)
#Convert dataframe to an mvabund object format
abund_mvabund <- mvabund(abund_rel_c)
#Perform the model
glm.data <- manyglm(abund_mvabund ~ Cystoseira.sp + Cymodocea.nodosa + Padina.pavonica + Halopteris.scoparia, family = "negative.binomial", data = biomass)
#Assumptions to check
#Plot of the residuals (we should see a random scatter of points, not linear or curvilinear relationship, or a fan shape)
plot(glm.data)
#Statistical test of the fitted model and each variable
glm_result <- anova(glm.data)
glm_result

####10. Figure 5 - IndVal analyses (indicator species)####
#Join all >0.05% filtered databases in one dataset
MOTUs <- rbind(abun_filt_art, abun_filt_moll, abun_filt_ann, abun_filt_por, abun_filt_cho, 
               abun_filt_cni, abun_filt_bry, abun_filt_ech, abun_filt_kin, abun_filt_nem, 
               abun_filt_neme, abun_filt_pla, abun_filt_rot)

#Arrange database
MOTUs_t <- as.data.frame(t(MOTUs))
rownames(MOTUs_t) <- NULL

#Create a vector with the order of the samples and different assemblages (each sample belongs to one group)
#1 = No Restored
#2 = Restored Forest
#3 = Cala Rotja Forest
#4 = Miami Forest
groups <- c("1", "1", "1", "1", "1", "1",
            "2", "2", "2", "2", "2", "2",
            "3", "3", "3", "3", "3", "3",
            "4", "4", "4", "4", "4", "4")

#Run the IndVal function with the previous set grouping
iva <- indval(MOTUs_t, groups)

#Specify the significance level for different parameters
#We set p-value < 0.02
gr <- iva$maxcls[iva$pval<=0.02] 
iv <- iva$indcls[iva$pval<=0.02] 
pv <- iva$pval[iva$pval<=0.02]
fr <- apply(MOTUs_t[,-1]>0, 2, sum)[iva$pval<=0.02]
#Summarise results in a dataframe
indvalsummary <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)

#Store the different information of the indval output (for each assemblage 1 to 4)
q <- cbind(iva$relabu, iva$relfrq, iva$indval)
#Columns 1-4 are relative abundance
#Columns 5-8 are relative frequency
#Columns 9-12 are the indval value

#Merge the IndVal summary with detailed indicator species metrics
all <- merge(indvalsummary, q, by="row.names")

#Keep those with an IndVal > 0.4
all.sub <- subset(all, indval > 0.4)
all.sub <- all.sub %>% 
  rename(id = Row.names)

#Reload metadata file containing MOTUs information
metadata<-read.csv2("MOTUs_data.csv", sep = ";", dec = ".", stringsAsFactors = T)
#Remove non-identified MOTUs (NOIDs)
metadata_1 <- metadata %>% 
  filter(!phylum_name == "Noid")
#Extract the "id" and "name_2" columns for merging (name_2 has a combination of the phylum, lower identified name and MOTUs id)
id_nam <- metadata_1 %>% 
  select(id, name_2)

#Merge the indval output with the id and name_2 information
all2 <- merge(all.sub, id_nam, by = "id")

#Keep specific columns: id, group, relative abundance (columns 6 to 9), and name_2
#The relative frequency could also be used
all.sub <- all2[,c(1,2,6:9, 18)]

#Convert the resulting dataset into a molten dataframe for visualization
all.sub.gg <- melt(all.sub, id.vars = c("id", "name_2", "group"))

#Change the grouping numbers to the correct assemblage name and store it as factor
all.sub.gg$variable<-as.character(all.sub.gg$variable)
all.sub.gg$variable[all.sub.gg$variable == "1"] <- "No Rest"
all.sub.gg$variable[all.sub.gg$variable == "2"] <- "Rest"
all.sub.gg$variable[all.sub.gg$variable == "3"] <- "Rotja"
all.sub.gg$variable[all.sub.gg$variable == "4"] <- "Miami"
all.sub.gg$variable<-as.factor(all.sub.gg$variable)
#Do the same for group
all.sub.gg$group<-as.character(all.sub.gg$group)
all.sub.gg$group[all.sub.gg$group == 1] <- "No Rest"
all.sub.gg$group[all.sub.gg$group == 2] <- "Rest"
all.sub.gg$group[all.sub.gg$group == 3] <- "Rotja"
all.sub.gg$group[all.sub.gg$group == 4] <- "Miami"
all.sub.gg$group<-as.factor(all.sub.gg$group)

#Set the rest of columns/variables as factors
all.sub.gg$id <- as.factor(all.sub.gg$id)
all.sub.gg$group <- as.factor(all.sub.gg$group)

#Set the order of factors for group and variable columns
all.sub.gg$group <- factor(all.sub.gg$group, levels = c("No Rest", "Rest", "Rotja", "Miami"))
all.sub.gg$variable <- factor(all.sub.gg$variable, levels = c("No Rest", "Rest", "Rotja", "Miami"))

#Set the order of factors for name_2 column in reverse sorted order
all.sub.gg$name_2 <- factor(all.sub.gg$name_2, levels = rev(sort(unique(all.sub.gg$name_2))))

#Create a plot to visualise the IndVal output
ggplot(data = all.sub.gg, aes(variable, name_2))+
  geom_point(aes(size=value, colour=group))+
  scale_color_manual(values = pal) +
  scale_size_area(limits = c(0, 1), breaks = c(0.01, 0.1, 0.5, 1), max_size=6)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = .8, vjust=.7))

ggsave("Figure_5.pdf", width=19, height=19, units= "cm")


#Remove all objects in the Environment tab
####11. Supplementary information, Figure S2 - Rarefaction and accumulation curves####
#Upload the data and data arrangements as previously done
metadata <- read.csv('metadata.csv', sep = ";", stringsAsFactors = T)
metadata$assemblage <- factor(metadata$assemblage, levels = c("Teulera_no_restaurat", "Teulera_restaurat", "Fornells_Cala_Roja", "Fornells_Miami"))
abundance <- read.csv('MOTUs_data.csv', sep = ';', stringsAsFactors = T)
metadata <- metadata[metadata$codi %in% colnames(abundance),]
abund <- abundance[,colnames(abundance) %in% metadata$codi]
row.names(abund) <- abundance$id
taxo <- abundance[,!colnames(abundance) %in% metadata$codi]
abund_t <- as.data.frame(t(abund)) 
pal <- c("#DD542C", "#C08E70", "#98C190", "#777539")

#Rarefaction curves of the number of MOTUs at increasing numbers of reads
#Find the minimum total number of reads across samples
raremin <- min(rowSums(abund_t)) 
raremin #sample with lower number of reads

#Create a vector with assemblage
assemblage <- metadata$assemblage

#Generate rarefraction curves + plot
#Supplementary Figure 2A
out1 <- rarecurve(abund_t, step = 20, sample = raremin, label = FALSE, xlim = c(0, 40000), groups = assemblage, col = pal, lwd = 2)
legend("topright", legend = levels(assemblage), col = pal, lty = 1)

#MOTU accumulation curves at increasing number of samples
abund_t <- as.data.frame(t(abund))
rownames(metadata) <- metadata$codi

#Calculate accumulation curves grouped by "assemblage"
accum <- accumcomp(abund_t, y=metadata, factor='assemblage', 
                     method='exact', conditioned=FALSE, plotit=FALSE)

#Convert accumulation curve output to long format for plotting
accum.long <- accumcomp.long(accum, ci=NA, label.freq=5)
head(accum.long)

#Reorder levels
accum.long$Grouping <- factor(accum.long$Grouping, levels = c("Teulera_no_restaurat", "Teulera_restaurat", "Fornells_Cala_Roja", "Fornells_Miami"))

#Do the plot
ggplot(data=accum.long, aes(x = Sites, y = Richness)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_line(aes(colour=Grouping), size=1) +
  labs(x = "Samples", y = "number of MOTUs") +  
  scale_colour_manual(values = pal) +
  theme(axis.text.x = element_text(angle = 360, vjust = 0.5, hjust=0.5))

ggsave("Supplementary_Figure_2B.pdf", width=22, height=12, units= "cm")

