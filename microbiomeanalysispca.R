#Adapted script from https://link.springer.com/chapter/10.1007/978-981-13-1534-3_10

library(zCompositions)
library(ALDEx2)
library(dplyr)
library(phyloseq)
setwd("~/Downloads/")

#load data
review<- read.csv("phyloseq_merged.csv", row.names=1, check.names = FALSE)
#transpose
review_t<-t(review)
#replace 0s with count multiplicative method w/ zCompositions package
review_r<-t(cmultRepl((review), method="CZM", output= "p-counts"))
#convert to proportions
review_p<-apply(review_r, 2, function(x){x/sum(x)})
#clr transform
names_add<- rownames(review_p)[ + order(apply(review_p, 1, sum), decreasing=T)]
review_pr<-review_p[names_add,]
review_clr<-t(apply(review_pr, 2, function(x){log(x)-mean(log(x))}))
#Singular Value Decomposition
review_PCX<- prcomp(review_clr)
sum(review_PCX$sdev[1:2]^2)/mvar(review_clr)
#create distance matrix and cluster data
dist<-dist(review_clr, method="euclidian")
hc<-hclust(dist, method="ward.D2")
re_order<-review_pr[,hc$order]
abund_PCX<- prcomp(abund_clr)

#permanova to determine if animal vs. plant is significant 
#p=0.001

review_dist<-dist(review_clr, method='euclidean')
permanovaoverall<-pairwise.adonis(review_dist, factors=meta_table$Category, perm = 999, p.adjust.m = 'bonferroni')

#permanova for meat v plant v dairy
pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
#1  Dairy vs NonAnimal  1  104833.6 159.8254 0.2267588   0.001      0.003   *
#  2     Dairy vs Animal  1  295079.1 478.7637 0.3806468   0.001      0.003   *
#  3 NonAnimal vs Animal  1  131229.3 151.0293 0.2005623   0.001      0.003   *
permanovaoverall<-pairwise.adonis(review_dist, factors=meta_table$dmoney, perm = 999, p.adjust.m = 'bonferroni')

#permanova in dairy 

#create otu table
otu<-read.csv("otutable.csv", row.names=1, check.names = FALSE)
tax<-read.csv("taxtable.csv", row.names=1, check.names = FALSE)
otu<-as.matrix(otu)
tax<-as.matrix(tax)
library("phyloseq")
OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax)
physeq = phyloseq(OTU, TAX)
sampledata = sample_data(data.frame(meta_table))
physeq1 = merge_phyloseq(physeq, sampledata)
#subset for only dairy information
dairy<-subset_samples(physeq1, Type=="Dairy")
dm<-as.data.frame(otu_table(dairy))
#clr 
names_add<- rownames(dm)[ + order(apply(dm, 1, sum), decreasing=T)]
review_pr<-dm[names_add,]review_clr<-t(apply(review_pr, 2, function(x){log(x)-mean(log(x))}))
samd<-as.data.frame(sample_data(dairy))
#permanova
pairwise.adonis(review_dist, factors=samd$Study, perm = 999, p.adjust.m = 'bonferroni')

#permanova in meat 
#subset animal info 
animal<-subset_samples(physeq1, Category=="Animal")
am<-as.data.frame(otu_table(animal))
#clr
names_add<- rownames(am)[ + order(apply(am, 1, sum), decreasing=T)]
review_pr<-am[names_add,]
review_clr<-t(apply(review_pr, 2, function(x){log(x)-mean(log(x))}))
samd<-as.data.frame(sample_data(animal))
#permanova
review_dist<-dist(review_clr, method='euclidean')
pairwise.adonis(review_dist, factors=samd$Study, perm = 999, p.adjust.m = 'bonferroni')

#pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
#1  Ben vs Gio  1  13253.01 15.28043 0.1903382   0.001      0.003   *
#  2 Ben vs Kang  1  35927.03 57.42052 0.1263599   0.001      0.003   *
#  3 Gio vs Kang  1  25546.93 43.31859 0.1058310   0.001      0.003   *


#permanova in plant 
#subset plant infoe
plant<-subset_samples(physeq1, Category=="NonAnimal")
pm<-as.data.frame(otu_table(plant))
names_add<- rownames(pm)[ + order(apply(pm, 1, sum), decreasing=T)]
#clr
review_pr<-pm[names_add,]
review_clr<-t(apply(review_pr, 2, function(x){log(x)-mean(log(x))}))
samd<-as.data.frame(sample_data(plant))

#permanova
review_dist<-dist(review_clr, method='euclidean')
pairwise.adonis(review_dist, factors=samd$Study, perm = 999, p.adjust.m = 'bonferroni')

#pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
#1  Einson vs Gio  1 15181.652 18.591297 0.19247383   0.001      0.006   *
#  2 Einson vs Wang  1 26559.311 28.872123 0.28622499   0.001      0.006   *
#  3  Einson vs Tan  1 60426.970 87.208027 0.38896033   0.001      0.006   *
#  4    Gio vs Wang  1 12383.405 21.426594 0.32749059   0.001      0.006   *
#  5     Gio vs Tan  1 20007.066 40.291129 0.26988294   0.001      0.006   *
#  6    Wang vs Tan  1  3849.186  6.996922 0.06361016   0.001      0.006   *


#aldex analysis
#upload metadata
meta_table<-read.csv("phyloseq_merged4.csv", row.names=1, check.names = FALSE)
# in a category if it says nonm animal it is plant based if not it is animal based
groups <- with(meta_table,ifelse(as.factor(Category)%in% c("NonAnimal"), c("PlantB"), c("AnimalB")))
#clr on transposed data
vdr <- aldex.clr(review_t, groups, mc.samples=128, verbose=FALSE)
#welchs t and wilcox rank sums test
vdr_t<-aldex.ttest(vdr, groups, paired.test = FALSE)
#calculate effect size
vdr_effect <- aldex.effect(vdr, groups, verbose=FALSE)
#merge the data into one data frame
vdr_all <- data.frame(vdr_t, vdr_effect)

# Identify significant taxa 
sig_by_both <- which(vdr_all$we.ep < 0.05 & vdr_all$wi.ep < 0.05)

#reach significance when the p-values are adjusted for multiple testing corrections using the Benjamini-Hochberg’s method.
sig_by_both_fdr <- which(vdr_all$we.eBH < 0.05 & vdr_all$wi.eBH < 0.05)
#export data
table <-xtable(vdr_all[sig_by_both_fdr, c(1:11)], caption="Table of significant taxa", digits=3,label="sig.table", align=c("l",rep("r",11) ))
print.xtable(table, type="html", file="Vdr_Table.html")

#make pca 

library(compositions)

row_review<-rownames(review_clr) #Make vector with sample names
pc_review<-as.data.frame(review_PCX$x[,1:2]) #Get PC1 and PC2 
metreview<-as.data.frame(bind_cols(pc_review,groups,meta_table$Study)) #Add metadata information
row.names(metreview)<-row_review

mvar_clr<-mvar(review_clr)

PCA <- ggplot(pc_review, aes(x=PC1,y=PC2, color=meta_table$Study, shape=groups))+
  geom_point(size=3)+ 
  geom_text(aes(label=row.names(metreview)),hjust=0, vjust=0) +
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
  theme(legend.position = 'bottom')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=13,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(review_PCX$sdev[1]^2/mvar_clr*100, digits=2), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(review_PCX$sdev[2]^2/mvar_clr*100, digits=2), "%", sep="")) +
  ggtitle("Food Post Harvest Bacteria")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')



#ANCOM taken/modified from https://github.com/FrederickHuangLin/ANCOM

library(readr)
library(tidyverse)

#reviewp<- read.csv("review_p.csv", row.names=1)
#otu_id = reviewp$`feature-id`
otu_data = data.frame(review_t, check.names = FALSE)
#rownames(reviewp) = otu_id

#preprocessing 
feature_table = otu_data; sample_var = "Sample.ID"; group_var = NULL
out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = FALSE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)

feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

# Step 2: ANCOM

main_var = "Category"; p_adj_method = "BH"; alpha = 0.05
adj_formula = NULL; rand_formula = NULL
t_start = Sys.time()
res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start # around 30s

ress=res$out
write.csv(ress, "/Users/susansinclair/Downloads/rres.csv")

# Step 3: Volcano Plot

# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
dat_ann = data.frame(x = min(res$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")

fig = res$fig +  
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)
fig  

#create phyloseq object from http://joey711.github.io/phyloseq/import-data.html#manual
otu<-read.csv("otutable.csv")
View(otu)
otu<-read.csv("otutable.csv", row.names=1, check.names = FALSE)
tax<-read.csv("taxtable.csv", row.names=1, check.names = FALSE)
otu<-as.matrix(otu)
tax<-as.matrix(tax)
library("phyloseq")
OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax)
physeq = phyloseq(OTU, TAX)
sampledata = sample_data(data.frame(meta_table))
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
plot(random_tree)
physeq1 = merge_phyloseq(physeq, sampledata, random_tree)
physeq1

#core microbiom based off of Leo Lahti, Sudarshan Shetty et al. https://microbiome.github.io/tutorials/Core.html

#install.packages("eulerr") # If not installed
library(eulerr)
library(microbiome)
#devtools::install_github('microsud/microbiomeutilities')
library(microbiomeutilities)

pseq<-physeq1
#convert to relatove abundance
pseq.rel <- microbiome::transform(pseq, "compositional")

category <- unique(as.character(meta(pseq.rel)$Category))


list_core <- c() # an empty object to store information

for (n in category){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq.rel, Category == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.01, # 0.001 in atleast 90% samples 
                         prevalence = 0.05)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

mycols <- c(nonCRC="#d6e2e9", CRC="#cbf3f0", H="#fcf5c7") 
plot(venn(list_core),
     fills = mycols)

print(list_core)

# format names
pseq.rel.f <- format_to_besthit(pseq.rel)
# check names
taxa_names(pseq.rel.f)[1:5]

list_core <- c() # an empty object to store information

for (n in category){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq.rel.f, Category == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection =0, # 0.001 in atleast 90% samples 
                         prevalence = 0.001)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

print(list_core)

#RA bubble plot for core microbiota adapted from https://github.com/joey711/phyloseq/issues/901
#transform into relative abundance and filter to a mean threshold

physeq3 = transform_sample_counts(physeq1, function(x) x / sum(x) )

#dataframe w/o condensed phyla
stom<-psmelt(physeq3)

#export dataframe
write.csv(stom, "stom.csv")

#average coremicrobiota families relative abundance in animal vs plant based groups in excel and import tab;e
data<-read.csv("stom.csv", check.names = FALSE)

#bubbleplot

bubbleplot_review<-ggplot(data, aes(x=Category, y=Family, size=Abundance, color=Category))+
  geom_point() +
  scale_size(range=c(.1,6), name = "Relative abundance (%)")+
  theme(axis.title = element_blank(), axis.text.x=element_text(color='black', size=8), axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=6)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text = element_text(size=10, face='bold'),
        panel.border = element_rect(color="black", fill=NA))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='viridis')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
bubbleplot_review









