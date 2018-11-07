---
title: "2016_cultivar_proj_up"
author: "Fang Liu"
date: "7/17/2018"
output: html_document
---


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Install required libraries

```{r echo=FALSE}
#install.packages("extrafont")
library(extrafont)
#install.packages('vegan')
library(vegan)
#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')
library(phyloseq)
#install.packages('ggplot2')
library(ggplot2)
#install.packages('reshape2')
library(reshape2)
#install.packages('dplyr')
library(dplyr)
#install.packages('plyr')
library(plyr)
#source('http://bioconductor.org/biocLite.R')
#biocLite('metagenomeSeq')
library('metagenomeSeq')
#install.packages("MInt")
library(MInt)
# R package with built in files stored in "/Library/Frameworks/R.framework/Versions/3.4/Resources/library/MInt/extdata/x.txt"

#install.packages("igraph")
library(igraph)
# call functions from a package explicitly using '::'
#detach(package:igraph)

# install.package("Hmisc")
library(Hmisc)
#install.packages("directlabels")
library(directlabels)
```

## Set working directory

```{r}
setwd('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up')
getwd()
```

## Upload mothur .shared and taxonomy file

```{r}
CV_phyloseq<-import_mothur(mothur_shared_file = "cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03_up.shared",mothur_constaxonomy_file = "cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.cons_up.taxonomy") 

CV_phyloseq # 175922 OTUs and 136 samples
sum(sample_sums(CV_phyloseq)) # 9945986
otu_table(CV_phyloseq)<-otu_table(CV_phyloseq)[,c(1,5:12,2:4,13,15:22,14,23:68,69,73:80,70:72,81,83:90,82,91:136)]

length(sample_names(CV_phyloseq)) # 136
Sample_ID<-sample_names(CV_phyloseq)

CV_meta<-read.csv("cultivar_meta.csv",head=TRUE,row.names = 1)
#CV_meta<-cbind(CV_meta,Treat=paste(CV_meta$Soil_type,CV_meta$Treatment,sep="_"))
CV_meta<-data.frame(CV_meta,Read_depth=sample_sums(CV_phyloseq))
identical(rownames(CV_meta),sample_names(CV_phyloseq))
dim(CV_meta)
head(CV_meta)
# change the default levels in a defined order
CV_meta$Treat<-factor(CV_meta$Treat,levels=c("Ag_Fresh", "Ag_Bulk", "Ag_WIL", "Ag_NNW", "Ag_W82", "Ag_DRT", "Ag_CNR", "Ag_SOJ","For_Fresh", "For_Bulk", "For_WIL", "For_NNW", "For_W82", "For_DRT", "For_CNR","For_SOJ"))
#Here the factor function just used to change the order of the levels of treatment.When display in the ggplot, the catogory will be in the order listed in the levels here instead of characterwise order.
CV_meta$Treatment<-factor(CV_meta$Treatment,levels = c('Fresh','Bulk','WIL','NNW','W82','DRT','CNR','SOJ'))
CV_meta$Class<-factor(CV_meta$Class,levels=c("Ag_Soil","Ag_WIL","Ag_NNW","Ag_W82","Ag_DRT","Ag_CNR","Ag_SOJ","For_Soil","For_WIL","For_NNW","For_W82","For_DRT","For_CNR","For_SOJ"))
CV_meta_phyloseq<-sample_data(CV_meta)

# cultivar_phyloseq
cultivar_phyloseq<-merge_phyloseq(CV_phyloseq,CV_meta_phyloseq)
cultivar_phyloseq
rank_names(cultivar_phyloseq)
colnames(tax_table(cultivar_phyloseq))<-c('Kingdom','Phylum','Class','Order','Family','Genus')
rank_names(cultivar_phyloseq)
r_cultivar_phyloseq<-transform_sample_counts(cultivar_phyloseq,function(x) x/sum(x))
sample_sums(r_cultivar_phyloseq)[1:5]

```

# Bar plot of sequencing depth across samples

```{r}
Read_depth<-sample_sums(cultivar_phyloseq)
summary(Read_depth)
tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Sequencing depth barplot.tiff', units="in", width=16, height=8, res=300)
CV_meta$Sample_ID<-factor(CV_meta$Sample_ID,levels = CV_meta$Sample_ID)
ggplot(data=CV_meta, aes(x=Sample_ID,y=Read_depth))+geom_bar(stat = "identity") +xlab("Sample_ID")+ylab("Sequencing depth") +
  scale_y_continuous(labels = scales::comma)+theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(family="Arial"), plot.title = element_text(family='Arial', hjust = .5,vjust = 3,size = 24,face="bold"), axis.text.x = element_text(family="Arial",size=8,angle=90),axis.text.y = element_text(family="Arial",size=14),axis.title.x = element_text(family="Arial",size=20,face="bold"),axis.title.y = element_text(family="Arial",size=20,face="bold"))
dev.off()
```


## Generate phyloseq objective for subsequent analysis


```{r}

# rarefied_cultivar_phyloseq

rarefied_cultivar_phyloseq<-rarefy_even_depth(cultivar_phyloseq,sample.size = min(sample_sums(cultivar_phyloseq)),rngseed = 1013,replace = FALSE,trimOTUs = TRUE,verbose = TRUE)
# min(sample_sums(cultivar_phyloseq)) is 19023
rarefied_cultivar_phyloseq #76725 (76725/175957=43.7% left) OTUs across 136 samples, with 6 taxonomic ranks
r_rarefied_cultivar_phyloseq<-transform_sample_counts(rarefied_cultivar_phyloseq,function(x) x/sum(x))
r_rarefied_cultivar_phyloseq

# Rhi_rarefied_cultivar_phyloseq

Rhi_rarefied_cultivar_phyloseq<-subset_samples(rarefied_cultivar_phyloseq,Compartment=="Rhizosphere") # 76725 OTUs across 102 samples
Rhi_rarefied_cultivar_phyloseq<-prune_taxa(taxa_sums(Rhi_rarefied_cultivar_phyloseq)>0,Rhi_rarefied_cultivar_phyloseq)
Rhi_rarefied_cultivar_phyloseq #63923 taxa across 102 samples
r_Rhi_rarefied_cultivar_phyloseq<-transform_sample_counts(Rhi_rarefied_cultivar_phyloseq,function(x) x/sum(x))

# Bulk_rarefied_cultivar_phyloseq
Bulk_rarefied_cultivar_phyloseq<-subset_samples(rarefied_cultivar_phyloseq,Compartment=="Bulk"|Compartment=="Fresh") 
Bulk_rarefied_cultivar_phyloseq
Bulk_rarefied_cultivar_phyloseq<-prune_taxa(taxa_sums(Bulk_rarefied_cultivar_phyloseq)>0,Bulk_rarefied_cultivar_phyloseq)
Bulk_rarefied_cultivar_phyloseq # 30059 taxa across 34 samples
r_Bulk_rarefied_cultivar_phyloseq<-transform_sample_counts(Bulk_rarefied_cultivar_phyloseq,function(x) x/sum(x))

# Ag_rarefied_cultivar_phyloseq

Ag_rarefied_cultivar_phyloseq<-subset_samples(rarefied_cultivar_phyloseq,Soil_type=="Agriculture")
Ag_rarefied_cultivar_phyloseq<-prune_taxa(taxa_sums(Ag_rarefied_cultivar_phyloseq)>0,Ag_rarefied_cultivar_phyloseq)
Ag_rarefied_cultivar_phyloseq # 43956 OTUs across 68 samples
r_Ag_rarefied_cultivar_phyloseq<-transform_sample_counts(Ag_rarefied_cultivar_phyloseq,function(x) x/sum(x))
r_Ag_rarefied_cultivar_phyloseq

# For_rarefied_cultivar_phyloseq

For_rarefied_cultivar_phyloseq<-subset_samples(rarefied_cultivar_phyloseq,Soil_type=="Forest")
For_rarefied_cultivar_phyloseq<-prune_taxa(taxa_sums(For_rarefied_cultivar_phyloseq)>0,For_rarefied_cultivar_phyloseq)
For_rarefied_cultivar_phyloseq # 36045 OTUs across 68 samples
r_For_rarefied_cultivar_phyloseq<-transform_sample_counts(For_rarefied_cultivar_phyloseq,function(x) x/sum(x))
r_For_rarefied_cultivar_phyloseq

# Ag_Rhi_rarefied_cultivar_phyloseq

Ag_Rhi_rarefied_cultivar_phyloseq<-subset_samples(rarefied_cultivar_phyloseq,Compartment=="Rhizosphere"&Soil_type=="Agriculture")
Ag_Rhi_rarefied_cultivar_phyloseq # 76864 OTUs across 51 samples
Ag_Rhi_rarefied_cultivar_phyloseq<-prune_taxa(taxa_sums(Ag_Rhi_rarefied_cultivar_phyloseq)>0,Ag_Rhi_rarefied_cultivar_phyloseq)
Ag_Rhi_rarefied_cultivar_phyloseq # 37014 OTUs across 51 samples
r_Ag_Rhi_rarefied_cultivar_phyloseq<-transform_sample_counts(Ag_Rhi_rarefied_cultivar_phyloseq,function(x) x/sum(x))

# For_Rhi_rarefied_cultivar_phyloseq

For_Rhi_rarefied_cultivar_phyloseq<-subset_samples(rarefied_cultivar_phyloseq,Compartment=="Rhizosphere"&Soil_type=="Forest")
For_Rhi_rarefied_cultivar_phyloseq 
For_Rhi_rarefied_cultivar_phyloseq<-prune_taxa(taxa_sums(For_Rhi_rarefied_cultivar_phyloseq)>0,For_Rhi_rarefied_cultivar_phyloseq)
For_Rhi_rarefied_cultivar_phyloseq # 29542 OTUs across 51 samples
r_For_Rhi_rarefied_cultivar_phyloseq<-transform_sample_counts(For_Rhi_rarefied_cultivar_phyloseq,function(x) x/sum(x))

```


## Generate phyloseq objectives on different taxonomy


```{r}
# genera_cultivar_phyloseq
genera_rarefied_cultivar_phyloseq<-tax_glom(rarefied_cultivar_phyloseq,taxrank = rank_names(rarefied_cultivar_phyloseq)[6],NArm=TRUE) 
genera_rarefied_cultivar_phyloseq  # 642 genus in 136 samples

genera_rarefied_cultivar_phyloseq<-filter_taxa(genera_rarefied_cultivar_phyloseq,function(x) sum(x)>50,TRUE) 
genera_rarefied_cultivar_phyloseq  # 348 genus

# r_genera_cultivar_phyloseq
r_genera_rarefied_cultivar_phyloseq<-transform_sample_counts(genera_rarefied_cultivar_phyloseq,function(x)x/sum(x))
r_genera_rarefied_cultivar_phyloseq


# family_cultivar_phyloseq
family_rarefied_cultivar_phyloseq<-tax_glom(rarefied_cultivar_phyloseq,taxrank = rank_names(rarefied_cultivar_phyloseq)[5],NArm=TRUE) 
family_rarefied_cultivar_phyloseq  # 250 family in 136 samples
family_rarefied_cultivar_phyloseq<-filter_taxa(family_rarefied_cultivar_phyloseq,function(x) sum(x)>50,TRUE) 
family_rarefied_cultivar_phyloseq  # 174  family
# r_family_cultivar_phylose
r_family_rarefied_cultivar_phyloseq<-transform_sample_counts(family_rarefied_cultivar_phyloseq,function(x)x/sum(x))
r_family_rarefied_cultivar_phyloseq

# order_cultivar_phyloseq
order_rarefied_cultivar_phyloseq<-tax_glom(rarefied_cultivar_phyloseq,taxrank = rank_names(rarefied_cultivar_phyloseq)[4],NArm=TRUE) 
order_rarefied_cultivar_phyloseq  # 123 order in 136 samples
order_rarefied_cultivar_phyloseq<-filter_taxa(order_rarefied_cultivar_phyloseq,function(x) sum(x)>50,TRUE) 
order_rarefied_cultivar_phyloseq  # 94 order left
# r_order_cultivar_phylose
r_order_rarefied_cultivar_phyloseq<-transform_sample_counts(order_rarefied_cultivar_phyloseq,function(x)x/sum(x))
r_order_rarefied_cultivar_phyloseq

# class_cultivar_phyloseq
class_rarefied_cultivar_phyloseq<-tax_glom(rarefied_cultivar_phyloseq,taxrank = rank_names(rarefied_cultivar_phyloseq)[3],NArm=TRUE) 
class_rarefied_cultivar_phyloseq  #  136 class in 136 samples
class_rarefied_cultivar_phyloseq<-filter_taxa(class_rarefied_cultivar_phyloseq,function(x) sum(x)>50,TRUE) 
class_rarefied_cultivar_phyloseq  # 61 left 
# r_class_cultivar_phylose
r_class_rarefied_cultivar_phyloseq<-transform_sample_counts(class_rarefied_cultivar_phyloseq,function(x)x/sum(x))
r_class_rarefied_cultivar_phyloseq


# phylum_cultivar_phyloseq
phylum_rarefied_cultivar_phyloseq<-tax_glom(rarefied_cultivar_phyloseq,taxrank = rank_names(rarefied_cultivar_phyloseq)[2],NArm=TRUE) 
phylum_rarefied_cultivar_phyloseq  #   phylum in 136 samples
phylum_rarefied_cultivar_phyloseq<-filter_taxa(phylum_rarefied_cultivar_phyloseq,function(x) sum(x)>50,TRUE) 
phylum_rarefied_cultivar_phyloseq  #  phylum
# r_phylum_cultivar_phylose
r_phylum_rarefied_cultivar_phyloseq<-transform_sample_counts(phylum_rarefied_cultivar_phyloseq,function(x)x/sum(x))
r_phylum_rarefied_cultivar_phyloseq

```

## Subset samples from WIL, NNW and W82 samples that collected at the same date

```{r}
# Ag_WIL_NNW_W82_rarefied_cultivar_phyloseq
Ag_WIL_NNW_W82_rarefied_cultivar_phyloseq<-subset_samples(Ag_Rhi_rarefied_cultivar_phyloseq,Treatment=="WIL"|Treatment=="NNW"|Treatment=="W82")
sample_names(Ag_WIL_NNW_W82_rarefied_cultivar_phyloseq) 
Ag_WIL_NNW_W82_rarefied_cultivar_phyloseq# 37014 OTUs across 27 samples
Ag_WIL_NNW_W82_rarefied_cultivar_phyloseq<-prune_taxa(taxa_sums(Ag_WIL_NNW_W82_rarefied_cultivar_phyloseq)>0,Ag_WIL_NNW_W82_rarefied_cultivar_phyloseq)
Ag_WIL_NNW_W82_rarefied_cultivar_phyloseq # 24934 OTUs across 27 samples
```

## Try different rarefaction at 50000 depth

```{r}
## Rarefy cultivar phyloseq using 50000 depth

# rarefied_50000_cultivar_phyloseq

dim(otu_table(cultivar_phyloseq))
less50000_samples<-as.matrix(otu_table(cultivar_phyloseq))[,sample_sums(cultivar_phyloseq)<50000]
more50000_phyloseq<-rarefy_even_depth(cultivar_phyloseq,sample.size = 50000,rngseed = 1013,replace = FALSE,trimOTUs = FALSE,verbose = TRUE)
identical(sample_names(more50000_phyloseq),names(which(sample_sums(cultivar_phyloseq)>=50000)))
rarefied_50000_cultivar_otu<-cbind(less50000_samples,as.matrix(otu_table(more50000_phyloseq)))
rarefied_50000_cultivar_otu_phyloseq<-otu_table(rarefied_50000_cultivar_otu,taxa_are_rows = TRUE)
rarefied_50000_cultivar_phyloseq<-merge_phyloseq(rarefied_50000_cultivar_otu_phyloseq,CV_meta_phyloseq,tax_table(cultivar_phyloseq))
rarefied_50000_cultivar_phyloseq # 175922 taxa across 136 samples
rarefied_50000_cultivar_phyloseq<-filter_taxa(rarefied_50000_cultivar_phyloseq,function(x) sum(x)>0,TRUE)
rarefied_50000_cultivar_phyloseq # 133045 taxa and 136 samples

sample_sums(rarefied_50000_cultivar_phyloseq)
```


## Remove all singletons before rarefaction

```{r}
# cultivar2_phyloseq (remove all singleton)

cultivar2_phyloseq<-filter_taxa(cultivar_phyloseq,function(x) sum(x)>1,prune = TRUE)
cultivar2_phyloseq # 63239 OTUs 
min(sample_sums(cultivar2_phyloseq)) # 18620

#rarefied_cultivar2_phyloseq

rarefied_cultivar2_phyloseq<-rarefy_even_depth(cultivar2_phyloseq,sample.size = min(sample_sums(cultivar2_phyloseq)),rngseed = 1013,replace = FALSE,trimOTUs = TRUE,verbose = TRUE)
rarefied_cultivar2_phyloseq # 45490 taxa across 136 samples
r_rarefied_cultivar2_phyloseq<-transform_sample_counts(rarefied_cultivar2_phyloseq,function(x) x/sum(x))

#Rhi_rarefied_cultivar2_phyloseq

Rhi_rarefied_cultivar2_phyloseq<-subset_samples(rarefied_cultivar2_phyloseq,Compartment=="Rhizosphere") 
Rhi_rarefied_cultivar2_phyloseq<-prune_taxa(taxa_sums(Rhi_rarefied_cultivar2_phyloseq)>0,Rhi_rarefied_cultivar2_phyloseq)
Rhi_rarefied_cultivar2_phyloseq #39678 taxa across 102 samples
r_Rhi_rarefied_cultivar2_phyloseq<-transform_sample_counts(Rhi_rarefied_cultivar2_phyloseq,function(x) x/sum(x))

# Ag_Rhi_rarefied_cultivar2_phyloseq
Ag_Rhi_rarefied_cultivar2_phyloseq<-subset_samples(rarefied_cultivar2_phyloseq,Compartment=="Rhizosphere"&Soil_type=="Agriculture")
Ag_Rhi_rarefied_cultivar2_phyloseq 
Ag_Rhi_rarefied_cultivar2_phyloseq<-prune_taxa(taxa_sums(Ag_Rhi_rarefied_cultivar2_phyloseq)>0,Ag_Rhi_rarefied_cultivar2_phyloseq)
Ag_Rhi_rarefied_cultivar2_phyloseq # 23509 OTUs across 51 samples
# r_Ag_Rhi_rarefied_cultivar2_phyloseq
r_Ag_Rhi_rarefied_cultivar2_phyloseq<-transform_sample_counts(Ag_Rhi_rarefied_cultivar2_phyloseq,function(x) x/sum(x))

# Ag_WIL_NNW_W82_rarefied_cultivar2_phyloseq
Ag_WIL_NNW_W82_rarefied_cultivar2_phyloseq<-subset_samples(Ag_Rhi_rarefied_cultivar2_phyloseq,Treatment=="WIL"|Treatment=="NNW"|Treatment=="W82")
sample_names(Ag_WIL_NNW_W82_rarefied_cultivar2_phyloseq) 
Ag_WIL_NNW_W82_rarefied_cultivar2_phyloseq
Ag_WIL_NNW_W82_rarefied_cultivar2_phyloseq<-prune_taxa(taxa_sums(Ag_WIL_NNW_W82_rarefied_cultivar2_phyloseq)>0,Ag_WIL_NNW_W82_rarefied_cultivar2_phyloseq)
Ag_WIL_NNW_W82_rarefied_cultivar2_phyloseq # 17732 OTUs across 27 samples

# For_Rhi_rarefied_cultivar2_phyloseq

For_Rhi_rarefied_cultivar2_phyloseq<-subset_samples(rarefied_cultivar2_phyloseq,Compartment=="Rhizosphere"&Soil_type=="Forest")
For_Rhi_rarefied_cultivar2_phyloseq 
For_Rhi_rarefied_cultivar2_phyloseq<-prune_taxa(taxa_sums(For_Rhi_rarefied_cultivar2_phyloseq)>0,For_Rhi_rarefied_cultivar2_phyloseq)
For_Rhi_rarefied_cultivar2_phyloseq # 18830 OTUs across 51 samples
# r_For_Rhi_rarefied_cultivar2_phyloseq
r_For_Rhi_rarefied_cultivar2_phyloseq<-transform_sample_counts(For_Rhi_rarefied_cultivar2_phyloseq,function(x) x/sum(x))

```


## Diversity index


```{r}
Diversity<-data.frame(estimate_richness(rarefied_cultivar_phyloseq,split=TRUE,measures = NULL),sample_data(rarefied_cultivar_phyloseq))
Diversity
identical(rownames(Diversity),sample_names(rarefied_cultivar_phyloseq))

# have a quick look at the color I used
treats_colors<-c('#a1b7b2','#b1b4cc','#ffee32','#07aeba','#f936f6','#b6cc0e','#ed5567','#3a44ff','#662f3b','#66434a','#316022','#fc0505',"#bc883a",'#05fff2','#b105fc','#210000')
#treats_colors<-c('#fc6000','#78b517','#ce0404','#ce5103','#cec002','#61ce01','#00c7ce','#000ace','#620581','#108ca5','#c300ce','#565648','#056835','#f71d78','#1dd6f7','#f25252')
treats=c(levels(CV_meta$Treat))

# pie plot of treat color
tiff("/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Treat_corlor_pie.tiff",units = 'in',width=10,height = 8,res=300,compression = 'lzw')
pie(rep (1,16),col=treats_colors,labels=paste(treats,treats_colors,sep=":"),radius = 1.0)
dev.off()

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Shannon diversity with fresh and bulk combined to soil.tiff', units="in", width=7, height=5, res=300)

ggplot(Diversity,aes(x=Treat,y=Shannon,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors,name="Treatment")+labs(title="Shannon Diversity",x="Treatment",y="Shannon Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Simpson diversity with fresh and bulk combined to soil.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Treat,y=Simpson,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values =treats_colors,name="Treatment")+labs(title="Simpson Diversity",x="Treatment",y="Simpson Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()


#tiff('Inverse Simpson diversity .tiff', units="in", width=7, height=5, res=300)
#ggplot(Diversity,aes(x=Treat,y=InvSimpson,fill=Treat))+geom_boxplot(size=0.5)+
#geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = c('#fc6000','#78b517','#ce0404','#ce5103','#cec002','#61ce01','#00c7ce','#000ace','#620581','#108ca5','#c300ce','#565648','#056835','#f71d78','#1dd6f7','#f25252'))+labs(title="Inverse Simpson Diversity",x="Treatment",y="InvSimpson Index ")+theme(plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=10,angle=60,hjust = 1),
#    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
#dev.off()

# Soil_type and compartment effects
mod_full<-lm(Shannon~Soil_type*Compartment,data=Diversity)
anova(mod_full)
# Both factors are significant, but no interaction effects on Shannon diversity

TukeyHSD(aov(Shannon~Soil_type*Compartment,data=Diversity),"Compartment",ordered=FALSE)

# All Agriculture samples  - comparison between all treatments in agriculture soil
TukeyHSD(aov(Shannon~Treat,data = Diversity[1:68,]),"Treat",ordered = FALSE)

# All Forest samples - comparison between all treatments in forest soil
TukeyHSD(aov(Shannon~Treat,data = Diversity[69:136,]),"Treat",ordered = FALSE)

# Soil_type vs Cultivar Rhi - compare rhizosphere diversity between soil types and cultivars
mod_Rhi<-lm(Shannon~Soil_type*Treatment,data=Diversity[Diversity$Compartment=="Rhizosphere",])
# Here Treatment within rhizospehre samples indicate the cultivar factors.
anova(mod_Rhi )
# There is significant interaction effects between soil type and cultivars on Shannon diversity.

# Cultivars effects on AgRhi
Ag_mod2<-lm(Shannon~Treatment,data=Diversity[Diversity$Compartment=="Rhizosphere"&Diversity$Soil_type=="Agriculture",])
anova(Ag_mod2)

Ag_aov<-aov(Shannon~Treatment,data=Diversity[Diversity$Compartment=="Rhizosphere"&Diversity$Soil_type=="Agriculture",])
summary(Ag_aov)
TukeyHSD(Ag_aov)

# Cultivars effects on ForRhi
Forest_mod2<-lm(Shannon~Treatment,data=Diversity[Diversity$Compartment=="Rhizosphere"&Diversity$Soil_type=="Forest",])
anova(Forest_mod2)
For_aov<-aov(Shannon~Treatment,data=Diversity[Diversity$Compartment=="Rhizosphere"&Diversity$Soil_type=="Forest",])
summary(For_aov)
TukeyHSD(For_aov) # drought tolerant cultivar diversity is significantly smaller than other cultivars.

## Correlation between Shannon diversity and network average density 

#---------Top 50 network--------

df<-data.frame(Shannon=aggregate(Shannon~Class,data=Diversity,mean),net6_edge_density=c(0.7649,0.4947,0.068,0.0842,0.1188,0.1361,0.1554,0.6302,0.2490,0.1961,0.0667,0.0774,0.1069,0.0976))
df$soil_type<-sub("_.*","",df$Shannon.Class)
rownames(df)<-df$Shannon.Class
lm<-lm(net6_edge_density~df$Shannon.Shannon,data=df)
summary(lm)
#Coefficients:
#                   Estimate Std. Error t value Pr(>|t|)  
#(Intercept)        -0.90948    0.60889  -1.494   0.1611  
#df$Shannon.Shannon  0.18128    0.09631   1.882   0.0843 .
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.2083 on 12 degrees of freedom
#Multiple R-squared:  0.228,	Adjusted R-squared:  0.1636 
#F-statistic: 3.543 on 1 and 12 DF,  p-value: 0.08426

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Correlation between network edge6 density and microbial community diversity_Arial.tiff', units="in", width=8, height=6, res=300)
ggplot(df,aes(Shannon.Shannon,net6_edge_density,label=rownames(df)))+geom_point(size=1.5)+ylim(0,0.8)+geom_text(family="Arial",size=3)+annotate(geom="text",x=5.7,y=0.6,label="Y = -0.91 + 0.18X \n Adj R^2=0.16, p=0.08",family="Arial",size=4 )+geom_smooth(method = "lm",se=FALSE,color="black",size=0.5) +ggtitle("\n Correlation between microbiome diversity and \n and network edge density - Top50")+xlab("Shannon diversity")+ylab("Network edge density")+ theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(family="Arial"), plot.title = element_text(family='Arial', hjust = .5,vjust = 3,size = 24,face="bold"), axis.text.x = element_text(family="Arial",size=12),axis.text.y = element_text(family="Arial",size=16),axis.title.x = element_text(family="Arial",size=18,face="bold"),axis.title.y = element_text(family="Arial",size=18,face="bold"))
dev.off()


#--------Global network based correlation ----------------

df2<-data.frame(Shannon=aggregate(Shannon~Class,data=Diversity,mean),net4_edge_density=c(0.0235,0.0158,0.0064,0.0093,0.0076,0.0081,0.0073,0.0149,0.0088,0.0077,0.0074,0.0073,0.0068,0.0075))

rownames(df2)<-df$Shannon.Class
lm2<-lm(net4_edge_density~Shannon.Shannon,data=df2)
summary(lm2)

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Correlation between network edge4 density and microbial community diversity.tiff', units="in", width=8, height=6, res=300)
ggplot(df2,aes(Shannon.Shannon,net4_edge_density,label=rownames(df2)))+geom_point()+ylim(0,0.03)+geom_text(family="Arial",size=3)+annotate(geom="text",x=5.7,y=0.0175,label="Y = -0.018 + 0.0045X \n AdjR2=0.26, p=0.036",family="Arial",size=4 )+geom_smooth(method = "lm",se=FALSE,color="black",size=0.5) +ggtitle("\n Correlation between microbiome diversity and \n and network edge density - Global network")+xlab("Shannon diversity")+ylab("Network edge density")+ theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(family="Arial"), plot.title = element_text(family='Arial', hjust = .5,vjust = 3,size = 24,face="bold"), axis.text.x = element_text(family="Arial",size=12),axis.text.y = element_text(family="Arial",size=16),axis.title.x = element_text(family="Arial",size=18,face="bold"),axis.title.y = element_text(family="Arial",size=18,face="bold"))
dev.off()
```




## Permanova test (permutational analysis of variance based on disdance matrix)


```{r}

## Here I want to test one thing. How Rhizosphere effect varied between soybean growing in forest compared with agriculture soil across soybean cultivars
r_Ag_SOJ_vs_soil_phyloseq<-subset_samples(r_rarefied_cultivar_phyloseq,Class=="Ag_SOJ"|Class=="Ag_Soil")
r_Ag_SOJ_vs_soil_phyloseq<-prune_taxa(taxa_sums(r_Ag_SOJ_vs_soil_phyloseq)>0,r_Ag_SOJ_vs_soil_phyloseq)
adonis2(t(otu_table(r_Ag_SOJ_vs_soil_phyloseq))~Compartment,data=data.frame(sample_data(r_Ag_SOJ_vs_soil_phyloseq)),permutations=999,by="margin") # 44.94% variance
r_Ag_WIL_vs_soil_phyloseq<-subset_samples(r_rarefied_cultivar_phyloseq,Class=="Ag_WIL"|Class=="Ag_Soil")
r_Ag_WIL_vs_soil_phyloseq<-prune_taxa(taxa_sums(r_Ag_WIL_vs_soil_phyloseq)>0,r_Ag_WIL_vs_soil_phyloseq)
adonis2(t(otu_table(r_Ag_WIL_vs_soil_phyloseq))~Compartment,data=data.frame(sample_data(r_Ag_WIL_vs_soil_phyloseq)),permutations=999,by="margin") #42.98%
r_Ag_NNW_vs_soil_phyloseq<-subset_samples(r_rarefied_cultivar_phyloseq,Class=="Ag_NNW"|Class=="Ag_Soil")
r_Ag_NNW_vs_soil_phyloseq<-prune_taxa(taxa_sums(r_Ag_NNW_vs_soil_phyloseq)>0,r_Ag_NNW_vs_soil_phyloseq)
adonis2(t(otu_table(r_Ag_NNW_vs_soil_phyloseq))~Compartment,data=data.frame(sample_data(r_Ag_NNW_vs_soil_phyloseq)),permutations=999,by="margin") #48.83%
r_Ag_DRT_vs_soil_phyloseq<-subset_samples(r_rarefied_cultivar_phyloseq,Class=="Ag_DRT"|Class=="Ag_Soil")
r_Ag_DRT_vs_soil_phyloseq<-prune_taxa(taxa_sums(r_Ag_DRT_vs_soil_phyloseq)>0,r_Ag_DRT_vs_soil_phyloseq)
adonis2(t(otu_table(r_Ag_DRT_vs_soil_phyloseq))~Compartment,data=data.frame(sample_data(r_Ag_DRT_vs_soil_phyloseq)),permutations=999,by="margin") #53.87%
r_Ag_CNR_vs_soil_phyloseq<-subset_samples(r_rarefied_cultivar_phyloseq,Class=="Ag_CNR"|Class=="Ag_Soil")
r_Ag_CNR_vs_soil_phyloseq<-prune_taxa(taxa_sums(r_Ag_CNR_vs_soil_phyloseq)>0,r_Ag_CNR_vs_soil_phyloseq)
adonis2(t(otu_table(r_Ag_CNR_vs_soil_phyloseq))~Compartment,data=data.frame(sample_data(r_Ag_CNR_vs_soil_phyloseq)),permutations=999,by="margin") #48.90%
r_Ag_W82_vs_soil_phyloseq<-subset_samples(r_rarefied_cultivar_phyloseq,Class=="Ag_W82"|Class=="Ag_Soil")
r_Ag_W82_vs_soil_phyloseq<-prune_taxa(taxa_sums(r_Ag_W82_vs_soil_phyloseq)>0,r_Ag_W82_vs_soil_phyloseq)
adonis2(t(otu_table(r_Ag_W82_vs_soil_phyloseq))~Compartment,data=data.frame(sample_data(r_Ag_W82_vs_soil_phyloseq)),permutations=999,by="margin") #46.70%


r_For_SOJ_vs_soil_phyloseq<-subset_samples(r_rarefied_cultivar_phyloseq,Class=="For_SOJ"|Class=="For_Soil")
r_For_SOJ_vs_soil_phyloseq<-prune_taxa(taxa_sums(r_For_SOJ_vs_soil_phyloseq)>0,r_For_SOJ_vs_soil_phyloseq)
adonis2(t(otu_table(r_For_SOJ_vs_soil_phyloseq))~Compartment,data=data.frame(sample_data(r_For_SOJ_vs_soil_phyloseq)),permutations=999,by="margin") # 57.11%

r_For_WIL_vs_soil_phyloseq<-subset_samples(r_rarefied_cultivar_phyloseq,Class=="For_WIL"|Class=="For_Soil")
r_For_WIL_vs_soil_phyloseq<-prune_taxa(taxa_sums(r_For_WIL_vs_soil_phyloseq)>0,r_For_WIL_vs_soil_phyloseq)
adonis2(t(otu_table(r_For_WIL_vs_soil_phyloseq))~Compartment,data=data.frame(sample_data(r_For_WIL_vs_soil_phyloseq)),permutations=999,by="margin") #52.11%

r_For_NNW_vs_soil_phyloseq<-subset_samples(r_rarefied_cultivar_phyloseq,Class=="For_NNW"|Class=="For_Soil")
r_For_NNW_vs_soil_phyloseq<-prune_taxa(taxa_sums(r_For_NNW_vs_soil_phyloseq)>0,r_For_NNW_vs_soil_phyloseq)
adonis2(t(otu_table(r_For_NNW_vs_soil_phyloseq))~Compartment,data=data.frame(sample_data(r_For_NNW_vs_soil_phyloseq)),permutations=999,by="margin") #59.29%

r_For_DRT_vs_soil_phyloseq<-subset_samples(r_rarefied_cultivar_phyloseq,Class=="For_DRT"|Class=="For_Soil")
r_For_DRT_vs_soil_phyloseq<-prune_taxa(taxa_sums(r_For_DRT_vs_soil_phyloseq)>0,r_For_DRT_vs_soil_phyloseq)
adonis2(t(otu_table(r_For_DRT_vs_soil_phyloseq))~Compartment,data=data.frame(sample_data(r_For_DRT_vs_soil_phyloseq)),permutations=999,by="margin") #55.66%
r_For_CNR_vs_soil_phyloseq<-subset_samples(r_rarefied_cultivar_phyloseq,Class=="For_CNR"|Class=="For_Soil")
r_For_CNR_vs_soil_phyloseq<-prune_taxa(taxa_sums(r_For_CNR_vs_soil_phyloseq)>0,r_For_CNR_vs_soil_phyloseq)
adonis2(t(otu_table(r_For_CNR_vs_soil_phyloseq))~Compartment,data=data.frame(sample_data(r_For_CNR_vs_soil_phyloseq)),permutations=999,by="margin") #57.09.90%
r_For_W82_vs_soil_phyloseq<-subset_samples(r_rarefied_cultivar_phyloseq,Class=="For_W82"|Class=="For_Soil")
r_For_W82_vs_soil_phyloseq<-prune_taxa(taxa_sums(r_For_W82_vs_soil_phyloseq)>0,r_For_W82_vs_soil_phyloseq)
adonis2(t(otu_table(r_For_W82_vs_soil_phyloseq))~Compartment,data=data.frame(sample_data(r_For_W82_vs_soil_phyloseq)),permutations=999,by="margin") #55.74%


## r_rarefied_cultivar_phyloseq (in fact the results are exactly the same when using rarefied_cultivar_phyloseq)

#1) #>>> 1) Overall Soil_type, compartment, and read depth

set.seed(1013)
adonis2(t(otu_table(r_rarefied_cultivar_phyloseq))~Soil_type+Compartment+Read_depth,data=CV_meta,permutations=999,by="margin")
#             Df SumOfSqs        F Pr(>F)    
#Soil_type     1  21.0950 259.7105  0.001 ***
#Compartment   2   2.0970  12.9086  0.001 ***
#Read_depth    1   0.1463   1.8015  0.136    
#Residual    131  10.6405  

adonis2(t(otu_table(r_rarefied_cultivar_phyloseq))~Soil_type*Compartment+Read_depth,data=CV_meta,permutations=999,by='margin')
#                       Df SumOfSqs       F Pr(>F)    
#Read_depth              1   0.0987  1.4327  0.215    
#Soil_type:Compartment   2   1.7530 12.7223  0.001 ***
#Residual              129   8.8875


adonis2(t(otu_table(r_rarefied_cultivar_phyloseq))~Soil_type:Compartment+Read_depth,data=CV_meta,permutations=999,by='margin') # Soil_type, Compartment and interaction together could explain 25.7201/(25.7201+0.0987+8.8875)=74.11% variance.

# In this way, Soil_type:Compartment is not interaction effects, instead it is Soil_type, Compartment as well as their interactions.

#                       Df SumOfSqs       F Pr(>F)    
#Read_depth              1   0.0987  1.4327  0.201    
#Soil_type:Compartment   5  25.7201 74.6647  0.001 ***
#Residual              129   8.8875

#>>> 2) Soil type impacts on Bulk and Fresh samples

adonis2(t(otu_table(r_Bulk_rarefied_cultivar_phyloseq))~Soil_type,data=CV_meta[CV_meta$Compartment=="Bulk"|CV_meta$Compartment=="Fresh",],permutations=999,by="margin") #81.46% variance
#          Df SumOfSqs      F Pr(>F)    
#Soil_type  1   6.9740 140.59  0.001 ***
#Residual  32   1.5874

#>>> 3) Soil type impacts on Rhizosphere samples

adonis2(t(otu_table(r_Rhi_rarefied_cultivar_phyloseq))~Soil_type,data=CV_meta[CV_meta$Compartment=="Rhizosphere",],permutations=999,by="margin") # 70.62% variance

#           Df SumOfSqs      F Pr(>F)    
#Soil_type   1  18.7208 240.31  0.001 ***
#Residual  100   7.7901 


#>>>>---Rhizosphere effects for agriculture and forest samples respectively

# r_Ag_rarefied_cultivar_phyloseq
identical(rownames(t(otu_table(r_Ag_rarefied_cultivar_phyloseq))),rownames(CV_meta[CV_meta$Soil_type=="Agriculture",]))
adonis2(t(otu_table(r_Ag_rarefied_cultivar_phyloseq))~Compartment,data=CV_meta[CV_meta$Soil_type=="Agriculture",],permutations=999,by="margin") #28.21%

#            Df SumOfSqs      F Pr(>F)    
#Compartment  2   2.2074 12.774  0.001 ***
#Residual    65   5.6164 

# r_For_rarefied_cultivar_phyloseq
adonis2(t(otu_table(r_For_rarefied_cultivar_phyloseq))~Compartment,data=data.frame(sample_data(r_For_rarefied_cultivar_phyloseq)),permutations=999,by="margin") #38.56% variance

#            Df SumOfSqs      F Pr(>F)    
#Compartment  2   2.1157 20.405  0.001 ***
#Residual    65   3.3698
```

## r_Rhi_rarefied_cultivar_phyloseq

```{r}

#1) Overall soil_type vs cultivars impacts on rhizosphere microbiome

set.seed(1013)
adonis2(t(otu_table(r_Rhi_rarefied_cultivar_phyloseq))~Soil_type+Treatment,data=data.frame(sample_data(r_Rhi_rarefied_cultivar_phyloseq)),permutations=999,by='margin') 

#          Df SumOfSqs        F Pr(>F)    
#Soil_type  1  18.7208 264.1470  0.001 ***
#Treatment  5   1.0572   2.9835  0.004 ** 
#Residual  95   6.7329


adonis2(t(otu_table(r_Rhi_rarefied_cultivar_phyloseq))~Soil_type:Treatment,data=data.frame(sample_data(r_Rhi_rarefied_cultivar_phyloseq)),permutations=999,by='margin')
# Here is the summary of the effects of both Soil_type, Treatment and their interactions
#                    Df SumOfSqs      F Pr(>F)    
#Soil_type:Treatment 11   20.470 27.724  0.001 ***
#Residual            90    6.041  

adonis2(t(otu_table(r_Rhi_rarefied_cultivar_phyloseq))~Soil_type*Treatment,data=data.frame(sample_data(r_Rhi_rarefied_cultivar_phyloseq)),permutations=999,by='margin')
# Here is the pure interaction effects
#                    Df SumOfSqs      F Pr(>F)  
#Soil_type:Treatment  5   0.6919 2.0617  0.035 *
#Residual            90   6.0410  

# 2) >>>>>  Host genotypes impacts on Ag_Rhi VS Forest_Rhi

## r_Ag_Rhi_rarefied_cultivar_phyloseq
set.seed(1013)
adonis2(t(otu_table(r_Ag_Rhi_rarefied_cultivar_phyloseq))~Treatment,data=data.frame(sample_data(r_Ag_Rhi_rarefied_cultivar_phyloseq)),permutations=999,by="margin")
#          Df SumOfSqs     F Pr(>F)    
#Treatment  5   1.1296 2.711  0.001 ***
#Residual  45   3.7501

Variantion<-adonis2(t(otu_table(r_Ag_Rhi_rarefied_cultivar_phyloseq))~Treatment,data=data.frame(sample_data(r_Ag_Rhi_rarefied_cultivar_phyloseq)),permutations=999,by=NULL)$SumOfSqs[1]/sum(adonis2(t(otu_table(r_Ag_Rhi_rarefied_cultivar_phyloseq))~Treatment,data=data.frame(sample_data(r_Ag_Rhi_rarefied_cultivar_phyloseq)),permutations=999,by=NULL)$SumOfSqs)
Variantion # 23.15% variation could be explained by the cultivars

#CAP analysis

cap_Ag_Rhi<-capscale(t(otu_table(r_Ag_Rhi_rarefied_cultivar_phyloseq))~Treatment,data=data.frame(sample_data(r_Ag_Rhi_rarefied_cultivar_phyloseq)),distance = 'bray',dfun = vegdist)
anova(cap_Ag_Rhi) #23.15% 

#         Df SumOfSqs     F Pr(>F)    
#Model     5   1.1296 2.711  0.001 ***
#Residual 45   3.7501

## r_For_Rhi_rarefied_cultivar_phyloseq
set.seed(1013)
adonis2(t(otu_table(r_For_Rhi_rarefied_cultivar_phyloseq))~Treatment,data=data.frame(sample_data(r_For_Rhi_rarefied_cultivar_phyloseq)),permutations=999,by=NULL)

#          Df SumOfSqs     F Pr(>F)    
#Treatment  5   1.1296 2.711  0.001 ***
#Residual  45   3.7501

Variantion<-adonis2(t(otu_table(r_For_Rhi_rarefied_cultivar_phyloseq))~Treatment,data=data.frame(sample_data(r_For_Rhi_rarefied_cultivar_phyloseq)),permutations=999,by=NULL)$SumOfSqs[1]/sum(adonis2(t(otu_table(r_For_Rhi_rarefied_cultivar_phyloseq))~Treatment,data=data.frame(sample_data(r_For_Rhi_rarefied_cultivar_phyloseq)),permutations=999,by=NULL)$SumOfSqs)
Variantion #21.29%

# CAP analysis
cap_For_Rhi<-capscale(t(otu_table(r_For_Rhi_rarefied_cultivar_phyloseq))~Treatment,data=data.frame(sample_data(r_For_Rhi_rarefied_cultivar_phyloseq)),distance = 'bray',dfun = vegdist)
anova(cap_For_Rhi) # 21.29 variance 

#         Df SumOfSqs      F Pr(>F)    
#Model     5  0.61955 2.4339  0.001 ***
#Residual 45  2.29092

```

## Collection time VS cultivars impacts on soybean rhizosphere microbiome assembly 

##  1) Across all rhizosphere samples
```{r}

adonis2(t(otu_table(r_Rhi_rarefied_cultivar_phyloseq))~Treatment+Time,data=data.frame(sample_data(r_Rhi_rarefied_cultivar_phyloseq)),permutations=999,by='margin') # Cultivar explain 1.3% variance and nonsignificant. Time could explain 62.23% variance.

#          Df SumOfSqs       F Pr(>F)    
#Treatment  2   0.3457  1.5265  0.162    
#Time      13  16.0551 10.9066  0.001 ***
#Residual  83   9.3985

adonis2(t(otu_table(r_Rhi_rarefied_cultivar_phyloseq))~Treatment*Time,data=data.frame(sample_data(r_Rhi_rarefied_cultivar_phyloseq)),permutations=999,by='margin')

#               Df SumOfSqs F Pr(>F)
#Treatment:Time  0   0.0000         
#Residual       83   9.3985
```


##   Cultivar vs Time impacts  

#2) r_Ag_Rhi_rarefied_cultivar_phyloseq
```{r}
# i) >>>>>>>  Collection Time effects

adonis2(t(otu_table(r_Ag_Rhi_rarefied_cultivar_phyloseq))~Time,data=data.frame(sample_data(r_Ag_Rhi_rarefied_cultivar_phyloseq)),permutations=999,by=NULL) #35.38% variance

#          Df SumOfSqs      F Pr(>F)    
#Model    10   1.7265 2.1902  0.001 ***
#Residual 40   3.1532

# ii) >>>>>>>  Cultivar impacts

adonis2(t(otu_table(r_Ag_Rhi_rarefied_cultivar_phyloseq))~Treatment,data=data.frame(sample_data(r_Ag_Rhi_rarefied_cultivar_phyloseq)),permutations=999,by=NULL) #35.38% variance

#          Df SumOfSqs     F Pr(>F)    
#Model     5   1.1296 2.711  0.001 ***
#Residual 45   3.7501 

## iii) >>>> Cultivar vs Collection time effects

identical(rownames(t(otu_table(r_Ag_Rhi_rarefied_cultivar_phyloseq))),rownames(data.frame(sample_data(r_Ag_Rhi_rarefied_cultivar_phyloseq))))
adonis2(t(otu_table(r_Ag_Rhi_rarefied_cultivar_phyloseq))~Treatment+Time,data=data.frame(sample_data(r_Ag_Rhi_rarefied_cultivar_phyloseq)),permutations=999,by='margin') # Both of cultivar and Time effects are significant, but time explained more variance.

#           Df SumOfSqs      F Pr(>F)    
#Treatment  1  0.17758 2.3275  0.003 ** 
#Time       6  0.77446 1.6918  0.001 ***
#Residual  39  2.97559

## iv) >>>> Conditional impacts by cultivars and Time by partial out the impacts of each other

# >>>> i) Cultivar impacts after partial out time

Ag_cultivar_vs_time_df<-data.frame(Cultivar=data.frame(sample_data(r_Ag_Rhi_rarefied_cultivar_phyloseq))$Treatment,Time=data.frame(sample_data(r_Ag_Rhi_rarefied_cultivar_phyloseq))$Time)
Ag_cultivar_vs_time_df

cap_treatment<-capscale(t(otu_table(r_Ag_Rhi_rarefied_cultivar_phyloseq))~Treatment+Condition(Time),data=data.frame(sample_data(r_Ag_Rhi_rarefied_cultivar_phyloseq)),distance = 'bray',dfun = vegdist)
anova(cap_treatment) # When partial out the impacts from collection time, the impacts from cultivar is still significant. But can only explain 5.6% of variance. 

#          Df SumOfSqs      F Pr(>F)    
#Model     1  0.17758 2.3275  0.001 ***
#Residual 39  2.97559 

# >>>> ii) Time impacts after partial out cultivars

cap_time<-capscale(t(otu_table(r_Ag_Rhi_rarefied_cultivar_phyloseq))~Time+Condition(Treatment),data=data.frame(sample_data(r_Ag_Rhi_rarefied_cultivar_phyloseq)),distance = 'bray',dfun = vegdist)
anova(cap_time) # Significant and 20.65% variance

#Df SumOfSqs      F Pr(>F)    
#Model     6  0.77446 1.6918  0.001 ***
#Residual 39  2.97559

##########################
#  Based on above analysis, our conclusion about collection time impacts are that they are significant and larger than cultivar because it totally overlapped with cultivar impacts. To evaluate the imapcts of time, I divided samples based on cultivars and analized the impacs of collecting time on each cultivar. For some of them, collecting time has significant impacts, but some not.
#########################

## To illustrate the pure cultivar impact by using samples collected at the same day between cultivars.

# 08-08 samples

r_rarefied_cultivar_08_08_phyloseq<-subset_samples(r_rarefied_cultivar_phyloseq,Time=="08_08"&Compartment=="Rhizosphere")
r_rarefied_cultivar_08_08_phyloseq<-filter_taxa(r_rarefied_cultivar_08_08_phyloseq,function(x) sum(x)>0,TRUE)
sample_names(r_rarefied_cultivar_08_08_phyloseq)

adonis2(t(otu_table(r_rarefied_cultivar_08_08_phyloseq))~Treatment,data=data.frame(sample_data(r_rarefied_cultivar_08_08_phyloseq)), permutations = 999,by="margin") # 23.4% variance could be explained by cultivars

#           Df SumOfSqs      F Pr(>F)   
#Treatment  1  0.17758 2.1388  0.007 **
#Residual   7  0.58122

#PCoA plot using samples from the same collecting time

r_rarefied_cultivar_08_08_phyloseq_bray<-vegdist(t(otu_table(r_rarefied_cultivar_08_08_phyloseq)),method="bray",binary = FALSE)
r_rarefied_cultivar_08_08_phyloseq_PCoA<-ordinate(r_rarefied_cultivar_08_08_phyloseq,method = "PCoA", r_rarefied_cultivar_08_08_phyloseq_bray)


tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/PCoA rarefied_cultivar_08_08.tiff', units="in", width=7, height=5, res=300)
plot_ordination(r_rarefied_cultivar_08_08_phyloseq,r_rarefied_cultivar_08_08_phyloseq_PCoA,type="samples",color="Treat")+geom_point(size=4)+scale_color_manual(values = c("#ffee32","#07aeba"))+ ggtitle("\n Soybean rhizosphere microbiome \n collected on 08/08/2016")+labs(x="PCoA1 [25.6%]",y="PCoA2 [19.3%]") +theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="bold",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title=element_text(size = 15,face="bold",family = "Arial"),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))
dev.off()
```

## Time effects for each for the cultivars
# Based on adonis and PCoA 

```{r}
## 1) >>>> Ag_WIL_rarefied_cultivar_phyloseq

Ag_WIL_rarefied_cultivar_phyloseq<-subset_samples(rarefied_cultivar_phyloseq,Treat=="Ag_WIL")
r_Ag_WIL_rarefied_cultivar_phyloseq<-transform_sample_counts(Ag_WIL_rarefied_cultivar_phyloseq,function(x) x/sum(x))
sample_data(r_Ag_WIL_rarefied_cultivar_phyloseq)
r_Ag_WIL_rarefied_cultivar_phyloseq<-prune_taxa(taxa_sums(r_Ag_WIL_rarefied_cultivar_phyloseq)>0,r_Ag_WIL_rarefied_cultivar_phyloseq)
r_Ag_WIL_rarefied_cultivar_phyloseq

adonis2(t(otu_table(r_Ag_WIL_rarefied_cultivar_phyloseq))~Time,data=data.frame(sample_data(r_Ag_WIL_rarefied_cultivar_phyloseq)),permutations = 999) # 14.9% variance, here shows that time effecs is not significant.

#         Df SumOfSqs      F Pr(>F)  
#Time      1  0.10719 1.4026  0.058 .
#Residual  8  0.61136  


## >>> PCoA 
r_Ag_WIL_rarefied_cultivar_phyloseq_bray<-vegdist(t(otu_table(r_Ag_WIL_rarefied_cultivar_phyloseq)),method="bray",binary = FALSE)
r_Ag_WIL_rarefied_cultivar_phyloseq_PCoA<-ordinate(r_Ag_WIL_rarefied_cultivar_phyloseq,method = "PCoA", r_Ag_WIL_rarefied_cultivar_phyloseq_bray)

plot_ordination(r_Ag_WIL_rarefied_cultivar_phyloseq,r_Ag_WIL_rarefied_cultivar_phyloseq_PCoA,type="samples",color="Time")+
  geom_point(size=4)+scale_shape_manual(values = c(17,16))+scale_color_manual(values = c("#838B8B","#e83535"))+ ggtitle("Soybean rhizosphere microbiome of r_Ag_WIL_rarefied_cultivar phyloseq") +theme(plot.title=(element_text(size=15,family = "Arial",face="bold",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title=element_text(size = 15,face="bold",family = "Arial"),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))

## 2) >>>> Ag_DRT_rarefied_cultivar_phyloseq

Ag_DRT_rarefied_cultivar_phyloseq<-subset_samples(rarefied_cultivar_phyloseq,Treat=="Ag_DRT")
r_Ag_DRT_rarefied_cultivar_phyloseq<-transform_sample_counts(Ag_DRT_rarefied_cultivar_phyloseq,function(x) x/sum(x))
sample_data(r_Ag_DRT_rarefied_cultivar_phyloseq)
r_Ag_DRT_rarefied_cultivar_phyloseq<-prune_taxa(taxa_sums(r_Ag_DRT_rarefied_cultivar_phyloseq)>0,r_Ag_DRT_rarefied_cultivar_phyloseq)
r_Ag_DRT_rarefied_cultivar_phyloseq

adonis2(t(otu_table(r_Ag_DRT_rarefied_cultivar_phyloseq))~Time,data=data.frame(sample_data(r_Ag_DRT_rarefied_cultivar_phyloseq)),permutations = 999) # 20.57% variance, here shows that time effecs is not significant.

#         Df SumOfSqs     F Pr(>F)  
#Time      1  0.12376 1.554  0.076 .
#Residual  6  0.47786 


## >>> PCoA 
r_Ag_DRT_rarefied_cultivar_phyloseq_bray<-vegdist(t(otu_table(r_Ag_DRT_rarefied_cultivar_phyloseq)),method="bray",binary = FALSE)
r_Ag_DRT_rarefied_cultivar_phyloseq_PCoA<-ordinate(r_Ag_DRT_rarefied_cultivar_phyloseq,method = "PCoA", r_Ag_DRT_rarefied_cultivar_phyloseq_bray)

plot_ordination(r_Ag_DRT_rarefied_cultivar_phyloseq,r_Ag_DRT_rarefied_cultivar_phyloseq_PCoA,type="samples",color="Time")+
  geom_point(size=4)+scale_shape_manual(values = c(17,16))+scale_color_manual(values = c("#838B8B","#e83535"))+ ggtitle("Soybean rhizosphere microbiome of r_Ag_DRT_rarefied_cultivar phyloseq") +theme(plot.title=(element_text(size=15,family = "Arial",face="bold",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title=element_text(size = 15,face="bold",family = "Arial"),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))

## 3) >>>> Ag_CNR_rarefied_cultivar_phyloseq

Ag_CNR_rarefied_cultivar_phyloseq<-subset_samples(rarefied_cultivar_phyloseq,Treat=="Ag_CNR")
r_Ag_CNR_rarefied_cultivar_phyloseq<-transform_sample_counts(Ag_CNR_rarefied_cultivar_phyloseq,function(x) x/sum(x))
sample_data(r_Ag_CNR_rarefied_cultivar_phyloseq)
r_Ag_CNR_rarefied_cultivar_phyloseq<-prune_taxa(taxa_sums(r_Ag_CNR_rarefied_cultivar_phyloseq)>0,r_Ag_CNR_rarefied_cultivar_phyloseq)
r_Ag_CNR_rarefied_cultivar_phyloseq

adonis2(t(otu_table(r_Ag_CNR_rarefied_cultivar_phyloseq))~Time,data=data.frame(sample_data(r_Ag_CNR_rarefied_cultivar_phyloseq)),permutations = 999) # 14.81 variance, here shows that time effecs is not significant.

#          Df SumOfSqs      F Pr(>F)
#Time      1  0.08012 1.2169  0.109
#Residual  7  0.46085  


## >>> PCoA 
r_Ag_CNR_rarefied_cultivar_phyloseq_bray<-vegdist(t(otu_table(r_Ag_CNR_rarefied_cultivar_phyloseq)),method="bray",binary = FALSE)
r_Ag_CNR_rarefied_cultivar_phyloseq_PCoA<-ordinate(r_Ag_CNR_rarefied_cultivar_phyloseq,method = "PCoA", r_Ag_CNR_rarefied_cultivar_phyloseq_bray)

plot_ordination(r_Ag_CNR_rarefied_cultivar_phyloseq,r_Ag_CNR_rarefied_cultivar_phyloseq_PCoA,type="samples",color="Time")+
  geom_point(size=4)+scale_shape_manual(values = c(17,16))+scale_color_manual(values = c("#838B8B","#e83535"))+ ggtitle("Soybean rhizosphere microbiome of r_Ag_CNR_rarefied_cultivar phyloseq") +theme(plot.title=(element_text(size=15,family = "Arial",face="bold",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title=element_text(size = 15,face="bold",family = "Arial"),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))

## 4) >>>> Ag_NNW_rarefied_cultivar_phyloseq

Ag_NNW_rarefied_cultivar_phyloseq<-subset_samples(rarefied_cultivar_phyloseq,Treat=="Ag_NNW")
r_Ag_NNW_rarefied_cultivar_phyloseq<-transform_sample_counts(Ag_NNW_rarefied_cultivar_phyloseq,function(x) x/sum(x))
sample_data(r_Ag_NNW_rarefied_cultivar_phyloseq)
r_Ag_NNW_rarefied_cultivar_phyloseq<-prune_taxa(taxa_sums(r_Ag_NNW_rarefied_cultivar_phyloseq)>0,r_Ag_NNW_rarefied_cultivar_phyloseq)
r_Ag_NNW_rarefied_cultivar_phyloseq

adonis2(t(otu_table(r_Ag_NNW_rarefied_cultivar_phyloseq))~Time,data=data.frame(sample_data(r_Ag_NNW_rarefied_cultivar_phyloseq)),permutations = 999) # 21.02 variance, here shows that time effecs is significant.

#          Df SumOfSqs      F Pr(>F)
#Time      1  0.13500 1.8631  0.009 **
#Residual  7  0.50723 


## >>> PCoA 
r_Ag_NNW_rarefied_cultivar_phyloseq_bray<-vegdist(t(otu_table(r_Ag_NNW_rarefied_cultivar_phyloseq)),method="bray",binary = FALSE)
r_Ag_NNW_rarefied_cultivar_phyloseq_PCoA<-ordinate(r_Ag_NNW_rarefied_cultivar_phyloseq,method = "PCoA", r_Ag_NNW_rarefied_cultivar_phyloseq_bray)

plot_ordination(r_Ag_NNW_rarefied_cultivar_phyloseq,r_Ag_NNW_rarefied_cultivar_phyloseq_PCoA,type="samples",color="Time")+
  geom_point(size=4)+scale_shape_manual(values = c(17,16))+scale_color_manual(values = c("#838B8B","#e83535"))+ ggtitle("Soybean rhizosphere microbiome of r_Ag_NNW_rarefied_cultivar phyloseq") +theme(plot.title=(element_text(size=15,family = "Arial",face="bold",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title=element_text(size = 15,face="bold",family = "Arial"),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))

## 5) >>>> Ag_SOJ_rarefied_cultivar_phyloseq

Ag_SOJ_rarefied_cultivar_phyloseq<-subset_samples(rarefied_cultivar_phyloseq,Treat=="Ag_SOJ")
r_Ag_SOJ_rarefied_cultivar_phyloseq<-transform_sample_counts(Ag_SOJ_rarefied_cultivar_phyloseq,function(x) x/sum(x))
sample_data(r_Ag_SOJ_rarefied_cultivar_phyloseq)
r_Ag_SOJ_rarefied_cultivar_phyloseq<-prune_taxa(taxa_sums(r_Ag_SOJ_rarefied_cultivar_phyloseq)>0,r_Ag_SOJ_rarefied_cultivar_phyloseq)
r_Ag_SOJ_rarefied_cultivar_phyloseq

adonis2(t(otu_table(r_Ag_SOJ_rarefied_cultivar_phyloseq))~Time,data=data.frame(sample_data(r_Ag_SOJ_rarefied_cultivar_phyloseq)),permutations = 999) # 26.26 variance, here shows that time effecs is significant.

#          Df SumOfSqs      F Pr(>F)
#Time      1  0.13500 1.8631  0.009 **
#Residual  7  0.50723 


## >>> PCoA 
r_Ag_SOJ_rarefied_cultivar_phyloseq_bray<-vegdist(t(otu_table(r_Ag_SOJ_rarefied_cultivar_phyloseq)),method="bray",binary = FALSE)
r_Ag_SOJ_rarefied_cultivar_phyloseq_PCoA<-ordinate(r_Ag_SOJ_rarefied_cultivar_phyloseq,method = "PCoA", r_Ag_SOJ_rarefied_cultivar_phyloseq_bray)

plot_ordination(r_Ag_SOJ_rarefied_cultivar_phyloseq,r_Ag_SOJ_rarefied_cultivar_phyloseq_PCoA,type="samples",color="Time")+
  geom_point(size=4)+scale_shape_manual(values = c(17,16))+scale_color_manual(values = c("#838B8B","#e83535"))+ ggtitle("Soybean rhizosphere microbiome of r_Ag_SOJ_rarefied_cultivar phyloseq") +theme(plot.title=(element_text(size=15,family = "Arial",face="bold",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title=element_text(size = 15,face="bold",family = "Arial"),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))


## 6) >>>> Ag_W82_rarefied_cultivar_phyloseq

Ag_W82_rarefied_cultivar_phyloseq<-subset_samples(rarefied_cultivar_phyloseq,Treat=="Ag_W82")
r_Ag_W82_rarefied_cultivar_phyloseq<-transform_sample_counts(Ag_W82_rarefied_cultivar_phyloseq,function(x) x/sum(x))
sample_data(r_Ag_W82_rarefied_cultivar_phyloseq)
r_Ag_W82_rarefied_cultivar_phyloseq<-prune_taxa(taxa_sums(r_Ag_W82_rarefied_cultivar_phyloseq)>0,r_Ag_W82_rarefied_cultivar_phyloseq)
r_Ag_W82_rarefied_cultivar_phyloseq

adonis2(t(otu_table(r_Ag_W82_rarefied_cultivar_phyloseq))~Time,data=data.frame(sample_data(r_Ag_W82_rarefied_cultivar_phyloseq)),permutations = 999) # 26.38 variance, here shows that time effecs is significant.

#          Df SumOfSqs      F Pr(>F)
#Time      1  0.22662 2.1497  0.029 *
#Residual  6  0.63251


## >>> PCoA 
r_Ag_W82_rarefied_cultivar_phyloseq_bray<-vegdist(t(otu_table(r_Ag_W82_rarefied_cultivar_phyloseq)),method="bray",binary = FALSE)
r_Ag_W82_rarefied_cultivar_phyloseq_PCoA<-ordinate(r_Ag_W82_rarefied_cultivar_phyloseq,method = "PCoA", r_Ag_W82_rarefied_cultivar_phyloseq_bray)

plot_ordination(r_Ag_W82_rarefied_cultivar_phyloseq,r_Ag_W82_rarefied_cultivar_phyloseq_PCoA,type="samples",color="Time")+
  geom_point(size=4)+scale_shape_manual(values = c(17,16))+scale_color_manual(values = c("#838B8B","#e83535"))+ ggtitle("Soybean rhizosphere microbiome of r_Ag_W82_rarefied_cultivar phyloseq") +theme(plot.title=(element_text(size=15,family = "Arial",face="bold",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title=element_text(size = 15,face="bold",family = "Arial"),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))

```


##   Cultivar vs Time impacts  

#3) r_For_Rhi_rarefied_cultivar_phyloseq

```{r}
# i) >>>>>>>  Collection Time effects
For_cultivar_vs_time_df<-data.frame(Cultivar=data.frame(sample_data(r_For_Rhi_rarefied_cultivar_phyloseq))$Treatment,Time=data.frame(sample_data(r_For_Rhi_rarefied_cultivar_phyloseq))$Time)

adonis2(t(otu_table(r_For_Rhi_rarefied_cultivar_phyloseq))~Time,data=data.frame(sample_data(r_For_Rhi_rarefied_cultivar_phyloseq)),permutations=999,by=NULL) #30.99% variance

#          Df SumOfSqs      F Pr(>F)    
#Model    10   1.7265 2.1902  0.001 ***
#Residual 40   3.1532

# ii) >>>>>>>  Cultivar impacts

adonis2(t(otu_table(r_For_Rhi_rarefied_cultivar_phyloseq))~Treatment,data=data.frame(sample_data(r_For_Rhi_rarefied_cultivar_phyloseq)),permutations=999,by=NULL) #35.38% variance

#          Df SumOfSqs     F Pr(>F)    
#Model     5   1.1296 2.711  0.001 ***
#Residual 45   3.7501 

## iii) >>>> Cultivar vs Collection time effects

identical(rownames(t(otu_table(r_For_Rhi_rarefied_cultivar_phyloseq))),rownames(data.frame(sample_data(r_For_Rhi_rarefied_cultivar_phyloseq))))
adonis2(t(otu_table(r_For_Rhi_rarefied_cultivar_phyloseq))~Treatment+Time,data=data.frame(sample_data(r_For_Rhi_rarefied_cultivar_phyloseq)),permutations=999,by='margin') # Both of cultivar and Time effects are significant, but time explained more variance.

#           Df SumOfSqs      F Pr(>F)    
#Treatment  1  0.17758 2.3275  0.003 ** 
#Time       6  0.77446 1.6918  0.001 ***
#Residual  39  2.97559

## iv) >>>> Conditional impacts by cultivars and Time by partial out the impacts of each other

# >>>> i) Cultivar impacts after partial out time

For_cultivar_vs_time_df<-data.frame(Cultivar=data.frame(sample_data(r_For_Rhi_rarefied_cultivar_phyloseq))$Treatment,Time=data.frame(sample_data(r_For_Rhi_rarefied_cultivar_phyloseq))$Time)
For_cultivar_vs_time_df

cap_treatment<-capscale(t(otu_table(r_For_Rhi_rarefied_cultivar_phyloseq))~Treatment+Condition(Time),data=data.frame(sample_data(r_For_Rhi_rarefied_cultivar_phyloseq)),distance = 'bray',dfun = vegdist)
anova(cap_treatment) # When partial out the impacts from collection time, the impacts from cultivar is still significant. But can only explain 5.6% of variance. 

#          Df SumOfSqs      F Pr(>F)    
#Model     1  0.17758 2.3275  0.001 ***
#Residual 39  2.97559 

# >>>> ii) Time impacts after partial out cultivars

cap_time<-capscale(t(otu_table(r_For_Rhi_rarefied_cultivar_phyloseq))~Time+Condition(Treatment),data=data.frame(sample_data(r_For_Rhi_rarefied_cultivar_phyloseq)),distance = 'bray',dfun = vegdist)
anova(cap_time) # Significant and 20.65% variance

#Df SumOfSqs      F Pr(>F)    
#Model     6  0.77446 1.6918  0.001 ***
#Residual 39  2.97559

##########################
#  Based on above analysis, our conclusion about collection time impacts are that they are significant and larger than cultivar because it totally overlapped with cultivar impacts. To evaluate the imapcts of time, I divided samples based on cultivars and analized the impacs of collecting time on each cultivar. For some of them, collecting time has significant impacts, but some not.
#########################

For_08_26_rarefied_cultivar_phyloseq<-subset_samples(For_Rhi_rarefied_cultivar_phyloseq,Time=="08_26")
For_08_26_rarefied_cultivar_phyloseq<-prune_taxa(taxa_sums(For_08_26_rarefied_cultivar_phyloseq)>0, For_08_26_rarefied_cultivar_phyloseq)
For_08_26_rarefied_cultivar_phyloseq
r_For_08_26_rarefied_cultivar_phyloseq<-transform_sample_counts(For_08_26_rarefied_cultivar_phyloseq,function(x) x/sum(x))


adonis2(t(otu_table(For_08_26_rarefied_cultivar_phyloseq))~Treatment,data=data.frame(sample_data(For_08_26_rarefied_cultivar_phyloseq)), permutations = 999,by="margin") # 31.30% variance could be explained by cultivars

#          Df SumOfSqs     F Pr(>F)  
#Treatment  1  0.16813 3.189  0.018 *
#Residual   7  0.36905

#PCoA plot using samples from the same collecting time

r_For_08_26_rarefied_cultivar_phyloseq_bray<-vegdist(t(otu_table(r_For_08_26_rarefied_cultivar_phyloseq)),method="bray",binary = FALSE)
r_For_08_26_rarefied_cultivar_phyloseq_PCoA<-ordinate(r_For_08_26_rarefied_cultivar_phyloseq,method = "PCoA", r_For_08_26_rarefied_cultivar_phyloseq_bray)


tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/PCoA rarefied_cultivar_For_08_26.tiff', units="in", width=7, height=5, res=300)
plot_ordination(r_For_08_26_rarefied_cultivar_phyloseq,r_For_08_26_rarefied_cultivar_phyloseq_PCoA,type="samples",color="Treat")+geom_point(size=4)+scale_color_manual(values = c("#05fff2","#b105fc"))+ ggtitle("\n Soybean rhizosphere microbiome \n collected on 08/26/2016")+labs(x="PCoA1 [39.0%]",y="PCoA2 [21.0%]") +theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="bold",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title=element_text(size = 15,face="bold",family = "Arial"),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))
dev.off()
```


## heatmap using only LDA2_OTU_ID


```{r}
LDA2_OTU_ID<-read.csv(file="uniq_LDA2_OTUID_Ag_and_For.csv",header = TRUE) # 1766 OTUs
LDA_rarefied_cultivar_phyloseq<-prune_taxa(as.character(LDA2_OTU_ID$OUT_ID),rarefied_cultivar_phyloseq)
plot_heatmap(LDA_rarefied_cultivar_phyloseq, sample.label="Sample_ID",taxa.label = "Genus")
```


## Constrained analysis of principle coordinate  Partial distance_based redundance analysis
## ------------------CAP analysis----------------------

```{r}
# r_rarefied_cultivar_phyloseq

rarefied_cultivar_cap_compartment<-capscale(t(otu_table(r_rarefied_cultivar_phyloseq))~Compartment+Condition(Soil_type+Read_depth+Time),data=CV_meta,distance = 'bray',dfun = vegdist)
rarefied_cultivar_cap_compartment
anova(rarefied_cultivar_cap_compartment)

rarefied_cultivar_cap_soil<-capscale(t(otu_table(r_rarefied_cultivar_phyloseq))~Soil_type+Condition(Treatment+Read_depth+Time),data=CV_meta,distance = 'bray',dfun = vegdist)
rarefied_cultivar_cap_soil
anova(rarefied_cultivar_cap_soil)

rarefied_cultivar_cap_treatment<-capscale(t(otu_table(r_rarefied_cultivar_phyloseq))~Treatment+Condition(Soil_type+Read_depth+Time),data=CV_meta,distance = 'bray',dfun = vegdist)
rarefied_cultivar_cap_treatment
anova(rarefied_cultivar_cap_treatment)


# r_Rhi_rarefied_cultivar_phyloseq

Rhi_rarefied_cultivar_cap_soil<-capscale(t(otu_table(r_Rhi_rarefied_cultivar_phyloseq))~Soil_type+Condition(Treatment+Read_depth+Time),data=data.frame(sample_data(r_Rhi_rarefied_cultivar_phyloseq)) ,distance = 'bray',dfun = vegdist)
Rhi_rarefied_cultivar_cap_soil
anova(Rhi_rarefied_cultivar_cap_soil)


Rhi_rarefied_cultivar_cap_treatment<-capscale(t(otu_table(r_Rhi_rarefied_cultivar_phyloseq))~Treatment+Condition(Soil_type+Read_depth+Time),data=data.frame(sample_data(r_Rhi_rarefied_cultivar_phyloseq)) ,distance = 'bray',dfun = vegdist)
Rhi_rarefied_cultivar_cap_treatment
anova(Rhi_rarefied_cultivar_cap_treatment)

Rhi_rarefied_cultivar_cap_time<-capscale(t(otu_table(r_Rhi_rarefied_cultivar_phyloseq))~Time+Condition(Soil_type+Read_depth+Treatment),data=data.frame(sample_data(r_Rhi_rarefied_cultivar_phyloseq)) ,distance = 'bray',dfun = vegdist)
Rhi_rarefied_cultivar_cap_time
anova(Rhi_rarefied_cultivar_cap_time)

# r_Ag_Rhi_rarefied_cultivar_phyloseq

Ag_Rhi_rarefied_cultivar_cap_treatment_p<-capscale(t(otu_table(r_Ag_Rhi_rarefied_cultivar_phyloseq))~Treatment+Condition(Read_depth+Time),data=data.frame(sample_data(r_Ag_Rhi_rarefied_cultivar_phyloseq)),distance = 'bray',dfun = vegdist)
Ag_Rhi_rarefied_cultivar_cap_treatment_p
anova(Ag_Rhi_rarefied_cultivar_cap_treatment_p)

Ag_Rhi_rarefied_cultivar_cap_treatment<-capscale(t(otu_table(r_Ag_Rhi_rarefied_cultivar_phyloseq))~Treatment,data=data.frame(sample_data(r_Ag_Rhi_rarefied_cultivar_phyloseq)),distance = 'bray',dfun = vegdist)
Ag_Rhi_rarefied_cultivar_cap_treatment
anova(Ag_Rhi_rarefied_cultivar_cap_treatment)

Ag_Rhi_rarefied_cultivar_cap_time_p<-capscale(t(otu_table(r_Ag_Rhi_rarefied_cultivar_phyloseq))~Time+Condition(Treatment+Read_depth),data=data.frame(sample_data(r_Ag_Rhi_rarefied_cultivar_phyloseq)),distance = 'bray',dfun = vegdist)
Ag_Rhi_rarefied_cultivar_cap_time_p
anova(Ag_Rhi_rarefied_cultivar_cap_time_p)

Ag_Rhi_rarefied_cultivar_cap_time<-capscale(t(otu_table(r_Ag_Rhi_rarefied_cultivar_phyloseq))~Time,data=data.frame(sample_data(r_Ag_Rhi_rarefied_cultivar_phyloseq)),distance = 'bray',dfun = vegdist)
Ag_Rhi_rarefied_cultivar_cap_time
anova(Ag_Rhi_rarefied_cultivar_cap_time)



# r_For_Rhi_rarefied_cultivar_phyloseq

For_Rhi_rarefied_cultivar_cap_treatment_p<-capscale(t(otu_table(r_For_Rhi_rarefied_cultivar_phyloseq))~Treatment+Condition(Read_depth+Time),data=data.frame(sample_data(r_For_Rhi_rarefied_cultivar_phyloseq)),distance = 'bray',dfun = vegdist)
For_Rhi_rarefied_cultivar_cap_treatment_p
anova(For_Rhi_rarefied_cultivar_cap_treatment_p)

For_Rhi_rarefied_cultivar_cap_treatment<-capscale(t(otu_table(r_For_Rhi_rarefied_cultivar_phyloseq))~Treatment,data=data.frame(sample_data(r_For_Rhi_rarefied_cultivar_phyloseq)),distance = 'bray',dfun = vegdist)
For_Rhi_rarefied_cultivar_cap_treatment
anova(For_Rhi_rarefied_cultivar_cap_treatment)

For_Rhi_rarefied_cultivar_cap_time_p<-capscale(t(otu_table(r_For_Rhi_rarefied_cultivar_phyloseq))~Time+Condition(Treatment+Read_depth),data=data.frame(sample_data(r_For_Rhi_rarefied_cultivar_phyloseq)),distance = 'bray',dfun = vegdist)
For_Rhi_rarefied_cultivar_cap_time_p
anova(For_Rhi_rarefied_cultivar_cap_time_p)

For_Rhi_rarefied_cultivar_cap_time<-capscale(t(otu_table(r_For_Rhi_rarefied_cultivar_phyloseq))~Time,data=data.frame(sample_data(r_For_Rhi_rarefied_cultivar_phyloseq)),distance = 'bray',dfun = vegdist)
For_Rhi_rarefied_cultivar_cap_time
anova(For_Rhi_rarefied_cultivar_cap_time)

```



## PCoA plot


```{r}
#cultivar_phyloseq

cultivar_phyloseq_bray<-vegdist(t(otu_table(cultivar_phyloseq)),method="bray",binary = FALSE)
cultivar_phyloseq_PCoA<-ordinate(cultivar_phyloseq,method = "PCoA", cultivar_phyloseq_bray)

plot_ordination(cultivar_phyloseq,cultivar_phyloseq_PCoA,type="samples",color="Treatment",shape = "Soil_type")+
  geom_point(size=4)+scale_shape_manual(values = c(17,16))+scale_color_manual(values = c("#838B8B","#0A0A0A","#007700","#bc1616","#020de2","#ab0ef9","#0ec8d6","#f77300"))+ ggtitle("Soybean rhizosphere microbiome of cultivar phyloseq")+labs(x="PCoA1 [57.6%]",y="PCoA2 [8.8%]") +theme(plot.title=(element_text(size=15,family = "Arial",face="bold",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title=element_text(size = 15,face="bold",family = "Arial"),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))

# pie plot for cultivar color for PCoA plot
tiff("/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Pie_plot_for_cultivar_colors.tiff",units = 'in',width = 10,height = 8,res = 300)
cultivars<-c("WIL","NNW","W82","DRT","CNR","SOJ")
cultivars_colors<-c("#007700","#bc1616","#020de2","#ab0ef9","#0ec8d6","#f77300")
pie(rep (1,6),col=cultivars_colors,labels=paste(cultivars,cultivars_colors,sep=":"),radius = 1.0)
dev.off()

tiff("/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Pie_plot_for_treatment_colors.tiff",units = 'in',width = 10,height = 8,res = 300)
treatments<-c("Fresh","Bulk","WIL","NNW","W82","DRT","CNR","SOJ")
treatments_colors<-c("#838B8B","#0A0A0A","#007700","#bc1616","#020de2","#ab0ef9","#0ec8d6","#f77300")
pie(rep (1,8),col=treatments_colors,labels=paste(treatments,treatments_colors,sep=":"),radius = 1.0)
dev.off()

# rarefied_phyloseq

rarefied_cultivar_phyloseq_bray<-vegdist(t(otu_table(rarefied_cultivar_phyloseq)),method="bray",binary = FALSE)
rarefied_cultivar_phyloseq_PCoA<-ordinate(rarefied_cultivar_phyloseq,method = "PCoA", rarefied_cultivar_phyloseq_bray)
ls(rarefied_cultivar_phyloseq_PCoA)
rarefied_cultivar_phyloseq_PCoA$trace 
#The amount of variance explained is equal to the trace of the matrix (sum of the diagonals of the decomposed correlation matrix).- http://www2.sas.com/proceedings/sugi30/203-30.pdf

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/PCoA rarefied_cultivar_19023.tiff', units="in", width=7, height=5, res=300)
#ggtitle("PCoA of soybean rhizosphere and bulk microbiome - rarefied_19023")+
plot_ordination(rarefied_cultivar_phyloseq,rarefied_cultivar_phyloseq_PCoA,type="samples",color="Treatment",shape = "Soil_type")+
  geom_point(size=4)+scale_shape_manual(values = c(17,16))+scale_color_manual(values = c("#838B8B","#0A0A0A","#007700","#bc1616","#020de2","#ab0ef9","#0ec8d6","#f77300"))+labs(x="PCoA1 [64.6%]",y="PCoA2 [7.2%]") +theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="bold",vjust=3,hjust = .5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title=element_text(size = 15,face="bold",family = "Arial"),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin=unit(c(1,1,1,1),"cm"))
dev.off()

# rarefied_50000_cultivar_phyloseq

rarefied_50000_cultivar_phyloseq_bray<-vegdist(t(otu_table(rarefied_50000_cultivar_phyloseq)),method="bray",binary = FALSE)
rarefied_50000_cultivar_phyloseq_PCoA<-ordinate(rarefied_50000_cultivar_phyloseq,method = "PCoA", rarefied_50000_cultivar_phyloseq_bray)
ls(rarefied_50000_cultivar_phyloseq_PCoA)
rarefied_50000_cultivar_phyloseq_PCoA$trace 
#The amount of variance explained is equal to the trace of the matrix (sum of the diagonals of the decomposed correlation matrix).- http://www2.sas.com/proceedings/sugi30/203-30.pdf

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/PCoA rarefied_cultivar_50000.tiff', units="in", width=7, height=5, res=300)
#+ ggtitle("PCoA of soybean rhizosphere and bulk microbiome_rarefied_50000")
plot_ordination(rarefied_50000_cultivar_phyloseq,rarefied_50000_cultivar_phyloseq_PCoA,type="samples",color="Treatment",shape = "Soil_type")+
  geom_point(size=4)+scale_shape_manual(values = c(17,16))+scale_color_manual(values = c("#838B8B","#0A0A0A","#007700","#bc1616","#020de2","#ab0ef9","#0ec8d6","#f77300"))+labs(x="PCoA1 [65.5%]",y="PCoA2 [7.2%]") +theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="bold",vjust=3,hjust = .5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title=element_text(size = 15,face="bold",family = "Arial"),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin=unit(c(1,1,1,1),"cm"))
dev.off()
#P$layers<-P$layers[-1] (when using hollow shapes to avoid double layers)

# rarefied_cultivar2_phyloseq

rarefied_cultivar2_phyloseq_bray<-vegdist(t(otu_table(rarefied_cultivar2_phyloseq)),method="bray",binary = FALSE)
rarefied_cultivar2_phyloseq_PCoA<-ordinate(rarefied_cultivar2_phyloseq,method = "PCoA", rarefied_cultivar2_phyloseq_bray)
ls(rarefied_cultivar2_phyloseq_PCoA)
rarefied_cultivar2_phyloseq_PCoA$trace 
#The amount of variance explained is equal to the trace of the matrix (sum of the diagonals of the decomposed correlation matrix).- http://www2.sas.com/proceedings/sugi30/203-30.pdf


tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/PCoA rarefied_cultivar_19023_singleton_remove.tiff', units="in", width=7, height=5, res=300)
#ggtitle("PCoA of soybean rhizosphere and bulk microbiome - rarefied_18138") +
plot_ordination(rarefied_cultivar2_phyloseq,rarefied_cultivar2_phyloseq_PCoA,type="samples",color="Treatment",shape = "Soil_type")+
  geom_point(size=4)+scale_shape_manual(values = c(17,16))+scale_color_manual(values = c("#838B8B","#0A0A0A","#007700","#bc1616","#020de2","#ab0ef9","#0ec8d6","#f77300"))+ labs(x="PCoA1 [65.4%]",y="PCoA2 [7.2%]") +theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="bold",vjust=3,hjust = .5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title=element_text(size = 15,face="bold",family = "Arial"),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin=unit(c(1,1,1,1),"cm"))
dev.off()

# Ag_WIL_NNW_W82_rarefied_cultivar_phyloseq

Ag_WIL_NNW_W82_rarefied_cultivar_phyloseq_bray<-vegdist(t(otu_table(Ag_WIL_NNW_W82_rarefied_cultivar_phyloseq)),method="bray",binary = FALSE)
Ag_WIL_NNW_W82_rarefied_cultivar_phyloseq_PCoA<-ordinate(Ag_WIL_NNW_W82_rarefied_cultivar_phyloseq,method = "PCoA", Ag_WIL_NNW_W82_rarefied_cultivar_phyloseq_bray)
ls(Ag_WIL_NNW_W82_rarefied_cultivar_phyloseq_PCoA)
Ag_WIL_NNW_W82_rarefied_cultivar_phyloseq_PCoA$trace 
#The amount of variance explained is equal to the trace of the matrix (sum of the diagonals of the decomposed correlation matrix).- http://www2.sas.com/proceedings/sugi30/203-30.pdf


tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/PCoA Ag_WIL_NNW_W82_rarefied_cultivar_axes2_and_3.tiff', units="in", width=7, height=5, res=300)
#+ ggtitle("PCoA of soybean rhizosphere and bulk microbiome - Ag_WIL_4_6_rarefied_19023") 
plot_ordination(Ag_WIL_NNW_W82_rarefied_cultivar_phyloseq,axes=c(2,3),Ag_WIL_NNW_W82_rarefied_cultivar_phyloseq_PCoA,type="samples",color="Time",shape = "Treatment")+
  geom_point(size=4)+scale_shape_manual(values = c(17,16,4))+labs(x="PCoA1 [10.5%]",y="PCoA2 [8.9]")+ theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="bold",vjust=3,hjust = .5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title.x=element_text(size = 15,face="bold",family = "Arial",margin = margin(t=10,b=0,l=0,r=0)),axis.title.y=element_text(size = 15,face="bold",family = "Arial",margin = margin(t=0,b=0,l=0,r=10)),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin=unit(c(1,1,1,1),"cm"))
dev.off()


tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/PCoA Ag_WIL_NNW_W82_rarefied_cultivar_axes1_and_2.tiff', units="in", width=7, height=5, res=300)
#+ ggtitle("PCoA of soybean rhizosphere and bulk microbiome - Ag_WIL_4_6_rarefied_19023") 
plot_ordination(Ag_WIL_NNW_W82_rarefied_cultivar_phyloseq,axes=c(1,2),Ag_WIL_NNW_W82_rarefied_cultivar_phyloseq_PCoA,type="samples",color="Time",shape = "Treatment")+
  geom_point(size=4)+scale_shape_manual(values = c(17,16,4))+labs(x="PCoA1 [19.4%]",y="PCoA2 [10.5%]")+theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="bold",vjust=3,hjust = .5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title.x=element_text(size = 15,face="bold",family = "Arial",margin = margin(t=10,b=0,l=0,r=0)),axis.title.y=element_text(size = 15,face="bold",family = "Arial",margin = margin(t=0,b=0,l=0,r=10)),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin=unit(c(1,1,1,1),"cm"))
dev.off()

# Ag_WIL_NNW_W82_rarefied_cultivar2_phyloseq

Ag_WIL_NNW_W82_rarefied_cultivar2_phyloseq_bray<-vegdist(t(otu_table(Ag_WIL_NNW_W82_rarefied_cultivar2_phyloseq)),method="bray",binary = FALSE)
Ag_WIL_NNW_W82_rarefied_cultivar2_phyloseq_PCoA<-ordinate(Ag_WIL_NNW_W82_rarefied_cultivar2_phyloseq,method = "PCoA", Ag_WIL_NNW_W82_rarefied_cultivar2_phyloseq_bray)
ls(Ag_WIL_NNW_W82_rarefied_cultivar2_phyloseq_PCoA)
Ag_WIL_NNW_W82_rarefied_cultivar2_phyloseq_PCoA$trace 
#The amount of variance explained is equal to the trace of the matrix (sum of the diagonals of the decomposed correlation matrix).- http://www2.sas.com/proceedings/sugi30/203-30.pdf


tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/PCoA Ag_WIL_NNW_W82_rarefied_cultivar2_axes2_and_3.tiff', units="in", width=7, height=5, res=300)
#+ ggtitle("PCoA of soybean rhizosphere and bulk microbiome - Ag_WIL_4_6_rarefied_19023") 
plot_ordination(Ag_WIL_NNW_W82_rarefied_cultivar2_phyloseq,axes=c(2,3),Ag_WIL_NNW_W82_rarefied_cultivar2_phyloseq_PCoA,type="samples",color="Time",shape = "Treatment")+
  geom_point(size=4)+scale_shape_manual(values = c(17,16,4))+labs(x="PCoA1 [10.5%]",y="PCoA2 [8.9]")+ theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="bold",vjust=3,hjust = .5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title.x=element_text(size = 15,face="bold",family = "Arial",margin = margin(t=10,b=0,l=0,r=0)),axis.title.y=element_text(size = 15,face="bold",family = "Arial",margin = margin(t=0,b=0,l=0,r=10)),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin=unit(c(1,1,1,1),"cm"))
dev.off()


tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/PCoA Ag_WIL_NNW_W82_rarefied_cultivar2_axes1_and_2.tiff', units="in", width=7, height=5, res=300)
#+ ggtitle("PCoA of soybean rhizosphere and bulk microbiome - Ag_WIL_4_6_rarefied_19023") 
plot_ordination(Ag_WIL_NNW_W82_rarefied_cultivar2_phyloseq,axes=c(1,2),Ag_WIL_NNW_W82_rarefied_cultivar2_phyloseq_PCoA,type="samples",color="Time",shape = "Treatment")+
  geom_point(size=4)+scale_shape_manual(values = c(17,16,4))+labs(x="PCoA1 [19.4%]",y="PCoA2 [10.5%]")+theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="bold",vjust=3,hjust = .5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title.x=element_text(size = 15,face="bold",family = "Arial",margin = margin(t=10,b=0,l=0,r=0)),axis.title.y=element_text(size = 15,face="bold",family = "Arial",margin = margin(t=0,b=0,l=0,r=10)),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin=unit(c(1,1,1,1),"cm"))
dev.off()

# Rhi_rarefied_cultivar_phyloseq

Rhi_rarefied_cultivar_phyloseq_bray<-vegdist(t(otu_table(Rhi_rarefied_cultivar_phyloseq)),method="bray",binary = FALSE)
Rhi_rarefied_cultivar_phyloseq_PCoA<-ordinate(Rhi_rarefied_cultivar_phyloseq,method = "PCoA", Rhi_rarefied_cultivar_phyloseq_bray)

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/PCoA_rarefied_Rhi_cultivar.tiff', units="in", width=7, height=5, res=300)
plot_ordination(Rhi_rarefied_cultivar_phyloseq,Rhi_rarefied_cultivar_phyloseq_PCoA,type="samples",color="Treatment",shape = "Soil_type")+
  geom_point(size=4)+scale_shape_manual(values = c(17,16))+scale_color_manual(values = c("#007700","#bc1616","#020de2","#ab0ef9","#0ec8d6","#f77300"))+ ggtitle("Rhi_rarefied_cultivar_phyloseq_weighted Bray-Curtis")+labs(x='PCoA1 [70.9%]',y='PCoA2 [4.2%]') +theme(plot.title=(element_text(size=15,family = "Arial",face="bold",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title=element_text(size = 15,face="bold",family = "Arial"),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin=unit(c(1,1,1,1),"cm"))
dev.off()

# Rhi_rarefied_cultivar2_phyloseq

Rhi_rarefied_cultivar2_phyloseq_bray<-vegdist(t(otu_table(Rhi_rarefied_cultivar2_phyloseq)),method="bray",binary = FALSE)
Rhi_rarefied_cultivar2_phyloseq_PCoA<-ordinate(Rhi_rarefied_cultivar2_phyloseq,method = "PCoA", Rhi_rarefied_cultivar2_phyloseq_bray)

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/PCoA_rarefied_Rhi_cultivar2.tiff', units="in", width=7, height=5, res=300)
plot_ordination(Rhi_rarefied_cultivar2_phyloseq,Rhi_rarefied_cultivar2_phyloseq_PCoA,type="samples",color="Treatment",shape = "Soil_type")+
  geom_point(size=4)+scale_shape_manual(values = c(17,16))+scale_color_manual(values = c("#007700","#bc1616","#020de2","#ab0ef9","#0ec8d6","#f77300"))+ ggtitle("Rhi_rarefied_cultivar2_phyloseq_weighted Bray-Curtis")+labs(x="PCoA1 [71.7%]",y="PCoA2 [4.3%]") +theme(plot.title=(element_text(size=15,family = "Arial",face="bold",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title=element_text(size = 15,face="bold",family = "Arial"),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin=unit(c(1,1,1,1),"cm"))
dev.off()

# For_Rhi_rarefied_cultivar_phyloseq

For_Rhi_rarefied_cultivar_phyloseq_bray<-vegdist(t(otu_table(For_Rhi_rarefied_cultivar_phyloseq)),method="bray",binary = FALSE)
For_Rhi_rarefied_cultivar_phyloseq_PCoA<-ordinate(For_Rhi_rarefied_cultivar_phyloseq,method = "PCoA", For_Rhi_rarefied_cultivar_phyloseq_bray)

ls(For_Rhi_rarefied_cultivar_phyloseq_PCoA)
For_Rhi_rarefied_cultivar_phyloseq_PCoA$trace #2.91

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/PCoA_rarefied_For_Rhi_cultivar.tiff', units="in", width=7, height=5, res=300)
#+ggtitle(" \n For_Rhi_rarefied_cultivar_phyloseq_weighted Bray-Curtis \n")
plot_ordination(For_Rhi_rarefied_cultivar_phyloseq,For_Rhi_rarefied_cultivar_phyloseq_PCoA,type="samples",color="Treatment") +ggtitle("For Rhizosphere PCoA")+
  geom_point(size=4)+ scale_color_manual(values = c("#007700","#bc1616","#020de2","#ab0ef9","#0ec8d6","#f77300"))+labs(x="PCoA1 [17.6%]",y="PCoA2 [9.5%]") +theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="bold",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title.x=element_text(size = 15,face="bold",family = "Arial",margin = margin(t=10,b=0,l=0,r=0)),axis.title.y=element_text(size = 15,face="bold",family = "Arial",margin = margin(t=0,b=0,l=0,r=10)),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1),"cm"))
dev.off()


# For_Rhi_rarefied_cultivar2_phyloseq

For_Rhi_rarefied_cultivar2_phyloseq_bray<-vegdist(t(otu_table(For_Rhi_rarefied_cultivar2_phyloseq)),method="bray",binary = FALSE)
For_Rhi_rarefied_cultivar2_phyloseq_PCoA<-ordinate(For_Rhi_rarefied_cultivar2_phyloseq,method = "PCoA", For_Rhi_rarefied_cultivar2_phyloseq_bray)

ls(For_Rhi_rarefied_cultivar2_phyloseq_PCoA)
For_Rhi_rarefied_cultivar2_phyloseq_PCoA$trace #2.67

plot_ordination(For_Rhi_rarefied_cultivar2_phyloseq,For_Rhi_rarefied_cultivar2_phyloseq_PCoA,type="samples",color="Treatment")+
  geom_point(size=4)+ scale_color_manual(values = c("#007700","#bc1616","#020de2","#ab0ef9","#0ec8d6","#f77300"))+ggtitle("For_Rhi_rarefied_cultivar2_phyloseq_weighted Bray-Curtis")+labs(x="PCoA1 [17.6%]",y="PCoA2 [9.4%]") +theme(plot.title=(element_text(size=15,family = "Arial",face="bold",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title.x=element_text(size = 15,face="bold",family = "Arial",margin = margin(t=10,b=0,l=0,r=0)),axis.title.x=element_text(size = 15,face="bold",family = "Arial",margin = margin(t=0,b=0,l=0,r=10)),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))



# Ag_Rhi_rarefied_cultivar_phyloseq

Ag_Rhi_rarefied_cultivar_phyloseq_bray<-vegdist(t(otu_table(Ag_Rhi_rarefied_cultivar_phyloseq)),method="bray",binary = FALSE)
Ag_Rhi_rarefied_cultivar_phyloseq_PCoA<-ordinate(Ag_Rhi_rarefied_cultivar_phyloseq,method = "PCoA", Ag_Rhi_rarefied_cultivar_phyloseq_bray)

ls(Ag_Rhi_rarefied_cultivar_phyloseq_PCoA)
Ag_Rhi_rarefied_cultivar_phyloseq_PCoA$trace #4.875

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/PCoA_rarefied_Ag_Rhi_cultivar.tiff', units="in", width=7, height=5, res=300)
#+ggtitle( "\n Ag_Rhi_rarefied_cultivar_phyloseq_weighted Bray-Curtis\n ")
plot_ordination(Ag_Rhi_rarefied_cultivar_phyloseq,Ag_Rhi_rarefied_cultivar_phyloseq_PCoA,type="samples",color="Treatment")+ggtitle("Ag Rhizosphere PCoA")+
  geom_point(size=4)+ scale_color_manual(values = c("#007700","#bc1616","#020de2","#ab0ef9","#0ec8d6","#f77300"))+labs(x="PCoA1 [14.6%]",y="PCoA2 [8.7%]")+theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="bold",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title.x=element_text(size = 15,family = "Arial",margin = margin(t=10,b=0,l=0,r=0)),axis.title.y=element_text(size = 15,family = "Arial",margin = margin(t=0,b=0,l=0,r=10)),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1),"cm"))
dev.off()
```


## CAP plot


```{r}

#rarefied_cultivar_phyloseq

rarefied_cultivar_phyloseq_cap<-ordinate(rarefied_cultivar_phyloseq,'CAP','bray',~Soil_type+Treatment)

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/CAP_Rhi_cultivar.tiff', units="in", width=7, height=5, res=300)
plot_ordination(rarefied_cultivar_phyloseq,rarefied_cultivar_phyloseq_cap,type="samples",color="Treatment",shape = "Soil_type")+
  geom_point(size=4)+scale_shape_manual(values = c(17,16))+scale_color_manual(values = c("#838B8B","#0A0A0A","#007700","#bc1616","#020de2","#ab0ef9","#0ec8d6","#f77300"))+ ggtitle("CAP plot of soybean rhizosphere and bulk microbiome \n ~Soil_type+Treatment") +theme(plot.title=(element_text(size=15,family = "Arial",face="bold",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title=element_text(size = 15,face="bold",family = "Arial"),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin=unit(c(1,1,1,1),"cm"))
dev.off()

# Forest rhizosphere CAP

#1) Treatment
For_Rhi_rarefied_cultivar_phyloseq_cap<-ordinate(For_Rhi_rarefied_cultivar_phyloseq,'CAP','bray',~Treatment)
ls(For_Rhi_rarefied_cultivar_phyloseq_cap)
For_Rhi_rarefied_cultivar_phyloseq_cap$adjust

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/CAP_For_Rhi_cultivar.tiff', units="in", width=7, height=5, res=300)
#+ggtitle("Forest soil ~ Treatment \n CAP plot of soybean rhizosphere microbiome")
plot_ordination(For_Rhi_rarefied_cultivar_phyloseq,For_Rhi_rarefied_cultivar_phyloseq_cap,type="samples",color="Treatment")+ggtitle("For Rhizosphere CAP")+geom_point(size=4)+scale_color_manual(values = c("#007700","#bc1616","#020de2","#ab0ef9","#0ec8d6","#f77300")) +theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="bold",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title.x=element_text(size = 15,face="bold",family = "Arial",margin=margin(t=10,b=0,l=0,r=0)),axis.title.y=element_text(size = 15,face="bold",family = "Arial",margin=margin(t=0,b=0,l=0,r=10)),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1),"cm"))
dev.off()

anova(For_Rhi_rarefied_cultivar_phyloseq_cap)

#2) Time
For_Rhi_rarefied_cultivar_phyloseq_cap<-ordinate(For_Rhi_rarefied_cultivar_phyloseq,'CAP','bray',~Time)
ls(For_Rhi_rarefied_cultivar_phyloseq_cap)
For_Rhi_rarefied_cultivar_phyloseq_cap$adjust

plot_ordination(For_Rhi_rarefied_cultivar_phyloseq,For_Rhi_rarefied_cultivar_phyloseq_cap,type="samples",color="Treatment",label = "Time")+
  geom_point(size=4)+scale_color_manual(values = c("#838B8B","#e83535","#10cc8e","#0A0A0A","#020de2","#ab0ef9","#75891b","#f9830c"))+ ggtitle("For_Rhi_rarefied_cultivar_phyloseq_weighted Bray-Curtis CAP plot Axis1 and Axis2") +theme(title=(element_text(size=15,family = "Arial",face="bold")),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title=element_text(size = 15,face="bold",family = "Arial"),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))

#3) Treatment +Time
#For_Rhi_rarefied_cultivar_phyloseq_cap<-ordinate(For_Rhi_rarefied_cultivar_phyloseq,'CAP','bray',~Treatment+Time)
#ls(For_Rhi_rarefied_cultivar_phyloseq_cap)
#For_Rhi_rarefied_cultivar_phyloseq_cap$adjust

#plot_ordination(For_Rhi_rarefied_cultivar_phyloseq,For_Rhi_rarefied_cultivar_phyloseq_cap,type="samples",color="Treatment",label = "Time")+
#  geom_point(size=4)+scale_color_manual(values = c("#838B8B","#e83535","#10cc8e","#0A0A0A","#020de2","#ab0ef9","#75891b","#f9830c"))+ ggtitle("For_Rhi_rarefied_cultivar_phyloseq_weighted Bray-Curtis CAP plot Axis1 and Axis2") +theme(title=(element_text(size=15,family = "Arial",face="bold")),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
#    axis.title=element_text(size = 15,face="bold",family = "Arial"),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))

# Agriculture rhizosphere CAP

#1) Treatment
Ag_Rhi_rarefied_cultivar_phyloseq_cap<-ordinate(Ag_Rhi_rarefied_cultivar_phyloseq,'CAP','bray',~Treatment)
ls(Ag_Rhi_rarefied_cultivar_phyloseq_cap)
Ag_Rhi_rarefied_cultivar_phyloseq_cap$adjust

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/CAP_Ag_Rhi_cultivar.tiff', units="in", width=7, height=5, res=300)
#+ ggtitle("Agriculture soil ~ Treatment \n CAP plot of soybean rhizosphere microbiome ")
plot_ordination(Ag_Rhi_rarefied_cultivar_phyloseq,Ag_Rhi_rarefied_cultivar_phyloseq_cap,type="samples",color="Treatment")+ggtitle("Ag Rhizosphere CAP")+
  geom_point(size=4)+scale_color_manual(values = c("#007700","#bc1616","#020de2","#ab0ef9","#0ec8d6","#f77300")) +theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="bold",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title.x=element_text(size = 15,face="bold",family = "Arial",margin=margin(t=10,b=0,l=0,r=0)),axis.title.y=element_text(size = 15,face="bold",family = "Arial",margin=margin(t=0,b=0,l=0,r=10)),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin=unit(c(1,1,1,1),"cm"))
dev.off()

#2) Time
Ag_Rhi_rarefied_cultivar_phyloseq_cap<-ordinate(Ag_Rhi_rarefied_cultivar_phyloseq,'CAP','bray',~Time)
ls(Ag_Rhi_rarefied_cultivar_phyloseq_cap)
Ag_Rhi_rarefied_cultivar_phyloseq_cap$adjust

plot_ordination(Ag_Rhi_rarefied_cultivar_phyloseq,Ag_Rhi_rarefied_cultivar_phyloseq_cap,type="samples",color="Time",label = "Treatment")+
  geom_point(size=4)+scale_color_manual(values = c("#8569D5", "#00CDCD", "orange","#DA5724", "#508578", "#CD9BCD",
   "#AD6F3B", "#673770","#D14285", "#652926", "#C84248"))+ ggtitle("Ag_Rhi_rarefied_cultivar_phyloseq_weighted Bray-Curtis CAP plot Axis1 and Axis2") +theme(title=(element_text(size=15,family = "Arial",face="bold")),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title=element_text(size = 15,face="bold",family = "Arial"),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))

#3) pure Treatment 
Ag_Rhi_rarefied_cultivar_phyloseq_cap<-ordinate(Ag_Rhi_rarefied_cultivar_phyloseq,'CAP','bray',~Treatment+Condition(Time))
#ls(Ag_Rhi_rarefied_cultivar_phyloseq_cap)
#Ag_Rhi_rarefied_cultivar_phyloseq_cap$adjust

plot_ordination(Ag_Rhi_rarefied_cultivar_phyloseq,Ag_Rhi_rarefied_cultivar_phyloseq_cap,type="samples",color="Treatment",label = "Time")+
geom_point(size=4)+scale_color_manual(values = c("#838B8B","#e83535","#10cc8e","#0A0A0A","#020de2","#ab0ef9","#75891b","#f9830c"))+ ggtitle("Ag_Rhi_rarefied_cultivar_phyloseq_weighted Bray-Curtis CAP plot Axis1 and Axis2") +theme(title=(element_text(size=15,family = "Arial",face="bold")),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
axis.title=element_text(size = 15,face="bold",family = "Arial"),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))

#4) pure Time
Ag_Rhi_rarefied_cultivar_phyloseq_cap<-ordinate(Ag_Rhi_rarefied_cultivar_phyloseq,'CAP','bray',~Time+Condition(Treatment))
#ls(Ag_Rhi_rarefied_cultivar_phyloseq_cap)
#Ag_Rhi_rarefied_cultivar_phyloseq_cap$adjust

plot_ordination(Ag_Rhi_rarefied_cultivar_phyloseq,Ag_Rhi_rarefied_cultivar_phyloseq_cap,type="samples",color="Treatment",label = "Time")+
geom_point(size=4)+scale_color_manual(values = c("#838B8B","#e83535","#10cc8e","#0A0A0A","#020de2","#ab0ef9","#75891b","#f9830c"))+ ggtitle("Ag_Rhi_rarefied_cultivar_phyloseq_weighted Bray-Curtis CAP plot Axis1 and Axis2") +theme(title=(element_text(size=15,family = "Arial",face="bold")),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
axis.title=element_text(size = 15,face="bold",family = "Arial"),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))


```


## rarefaction plot


```{r}
long_cultivar_otu<-as.data.frame(otu_table(cultivar_phyloseq))
dim(long_cultivar_otu)
short_cultivar_otu<-t(long_cultivar_otu)
dim(short_cultivar_otu)

# rarefaction with minimum read depth

cultivar_otu_rarefy<-rarefy(short_cultivar_otu,sample=seq(0,min(sample_sums(cultivar_phyloseq)),by=100))
dim(cultivar_otu_rarefy)
rownames(cultivar_otu_rarefy)<-rownames(short_cultivar_otu)
cultivar_otu_rarefy[1:10,1:10]
cultivar_otu_rarefy<-data.frame(subsample=seq(0,min(sample_sums(cultivar_phyloseq)),by=100),t(cultivar_otu_rarefy))
cultivar_otu_rarefy[1:10,1:10]
melt_cultivar_otu_rarefy<-melt(cultivar_otu_rarefy,id.vars = 'subsample')
head(melt_cultivar_otu_rarefy)
write.csv(melt_cultivar_otu_rarefy,file="melt_cultivar_otu_rarefy.csv")
melt_cultivar_otu_rarefy<-read.csv(file="melt_cultivar_otu_rarefy.csv")

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Rarefaction plot with minimum depth.tiff', units="in", width=16, height=8, res=600)
ggplot(melt_cultivar_otu_rarefy,aes(x=subsample,y=value,colour=variable))+geom_line(aes(x=subsample,y=value,group=variable),size=0.4)+
geom_dl(aes(label = variable), method = list(dl.combine( "last.points"), cex = 0.2))+ylab("Bactria OTU numbers")+xlab("Sequencing depth")+theme(panel.background = element_rect(fill="white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),title=(element_text(size=15,family = "Arial",face="bold")),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title=element_text(size = 15,face="bold",family = "Arial"),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=5,family = "Arial")),legend.position="none")+ggtitle("Rarefaction curve with minimum read depth")
dev.off()

# rarefaction with maximum read depth
# This rarefaction table is generated from mothur using raw OTU table and rarefied with interval of 100
rarefaction_df<-read.csv('cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.0.03.pick.groups.rarefaction.csv',header = TRUE)
head(rarefaction_df)
rarefaction_df<-rarefaction_df[c(1,which(rarefaction_df$numsampled%%100==0)),]
rarefaction_df_up<-rarefaction_df[,-c(grep("lci|hci",colnames(rarefaction_df),perl=TRUE))]
head(rarefaction_df_up)

melt_cultivar_otu_rarefy_max<-melt(rarefaction_df_up,id.vars = 'numsampled')
head(melt_cultivar_otu_rarefy_max)

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Rarefaction plot with maximum depth.tiff', units="in", width=16, height=8, res=600)
ggplot(melt_cultivar_otu_rarefy_max,aes(x=numsampled,y=value,colour=variable))+geom_line(aes(x=numsampled,y=value,group=variable),size=0.4)+
geom_dl(aes(label = variable), method = list(dl.combine("last.points"), cex = 0.2))+ylab("Bactria OTU numbers")+xlab("Sequencing depth")+theme(panel.background = element_rect(fill="white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), title=(element_text(size=15,family = "Arial",face="bold")),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),
    axis.title=element_text(size = 15,face="bold",family = "Arial"),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=5,family = "Arial")),legend.position="none")+ggtitle("Rarefaction curve with minimum read depth")
dev.off()
```


## Stacked barplot at Phylum level


```{r}
## prepare dataset

phylum_cultivar_phyloseq<-tax_glom(rarefied_cultivar_phyloseq,taxrank = rank_names(rarefied_cultivar_phyloseq)[2],NArm=TRUE) 
r_phylum_cultivar_phyloseq<-transform_sample_counts(phylum_cultivar_phyloseq,function(x) x/sum(x))
filter_r_phylum_cultivar_phyloseq<- filter_taxa(r_phylum_cultivar_phyloseq, function(x) sum(x>0.01)>0.2*length(x),prune = TRUE)
filter_r_phylum_cultivar_phyloseq
sample_sums(filter_r_phylum_cultivar_phyloseq)
tax_table(filter_r_phylum_cultivar_phyloseq)[,2]
Others<-1-sample_sums(filter_r_phylum_cultivar_phyloseq)
phylum_cultivar<-rbind(otu_table(filter_r_phylum_cultivar_phyloseq),Others)
phylum_cultivar[,1:10]
rownames(phylum_cultivar)<-c(tax_table(filter_r_phylum_cultivar_phyloseq)[,2],"Others")
phylum_cultivar
phylum_cultivar<-data.frame(Treat=as.factor(CV_meta$Treat),t(phylum_cultivar))
phylum_cultivar<-phylum_cultivar%>% group_by(Treat) %>%
  summarise_all(funs(mean))
phylum_cultivar
rowSums(phylum_cultivar[,2:13])
phylum_cultivar_melt<-melt(phylum_cultivar)
colnames(phylum_cultivar_melt)<-c('Treat','Phylum','Relative_abundance')
phylum_cultivar_melt

# quick look at the phylum colors
#phylum_color<-c("#8569D5", "#00CDCD", "orange","#DA5724", "#508578", "#CD9BCD", "#AD6F3B", "#673770","#D14285", "#652926", "#C84248","#CBD588")
# match the phylum color definition with that of phylogenetic tree and network.

phylum_color<-c("#f2701a", "#b6cc0e", "#ed5567","#07aeba", "#3a44ff", "#316022", "#f936f6", "#9b9696","#85f785", "#ffee32", "#723434","#990101")

phylum_level<-c("Verrucomicrobia","Proteobacteria","Acidobacteria","Bacteroidetes","Actinobacteria","Firmicutes","TM7", "Bacteria_unclassified", "Chloroflexi","Planctomycetes", "Gemmatimonadetes","Others")

tiff("/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/phylum_color_pie_plot.tiff",units = 'in',width = 15,height = 10,res = 300,compression = 'lzw')
pie(rep (1,12),col=phylum_color,labels=paste(phylum_level,phylum_color,sep=":"),radius = 1.0)
dev.off()

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Microbial community composition at Phylum level.tiff', units="in", width=7, height=5, res=300)

## plot barplot

ggplot(phylum_cultivar_melt, aes(x =Treat , y = Relative_abundance,fill=Phylum))+geom_bar(stat='identity',colour="grey",size=0.1)+labs(x="Treatment",y="Relative Abundance")+scale_fill_manual(values=phylum_color)+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), title=(element_text(size=15,family = "Arial",face="bold")),axis.text.x = element_text(angle=60,hjust = 1),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=11,family="Arial"),
    axis.title.x=element_text(size = 15,face="bold",family = "Arial",margin=margin(t=10,b=0,l=0,r=0)),axis.title.y=element_text(size = 15,face="bold",family = "Arial",margin=margin(t=0,b=0,l=0,r=10)),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))

dev.off()
```


## Heatmap

```{r}
# rarefied_cultivar_phyloseq

prune_rarefied_cultivar_phyloseq<- prune_taxa(names(sort(taxa_sums(rarefied_cultivar_phyloseq),TRUE)[1:200]), rarefied_cultivar_phyloseq)
prune_rarefied_cultivar_phyloseq
plot_heatmap(prune_rarefied_cultivar_phyloseq, sample.label="Sample_ID",sample.order = sample_names(rarefied_cultivar_phyloseq)[c(1:12,64:68,13:63,69:80,132:136,81:131)],taxa.label = "Phylum")

#rarefied_cultivar_phyloseq_Protoebacteria
rarefied_cultivar_phyloseq_Protoebacteria<-subset_taxa(rarefied_cultivar_phyloseq,Phylum=="Proteobacteria")

prune_rarefied_cultivar_phyloseq_Proteobacteria<- prune_taxa(names(sort(taxa_sums(rarefied_cultivar_phyloseq_Protoebacteria),TRUE)[1:100]), rarefied_cultivar_phyloseq_Protoebacteria)
prune_rarefied_cultivar_phyloseq_Proteobacteria
plot_heatmap(prune_rarefied_cultivar_phyloseq_Proteobacteria, sample.label="Sample_ID",sample.order = sample_names(rarefied_cultivar_phyloseq_Protoebacteria)[c(1:12,64:68,13:63,69:80,132:136,81:131)],taxa.label = "Genus")

```


## Boxplot for differential abundant genus, family, order,class and phylum


```{r}
# Pseudoxanthomonas genus
Pseudoxanthomonas_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Pseudoxanthomonas",]
Pseudoxanthomonas_genus_ab<-data.frame(Pseudoxanthomonas=t(Pseudoxanthomonas_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Pseudoxanthomonas_genus_ab)<-c('Pseudoxanthomonas','Treat')
Pseudoxanthomonas_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Pseudoxanthomonas genus.tiff', units="in", width=6, height=5, res=300)

ggplot(Pseudoxanthomonas_genus_ab,aes(x=Treat,y=Pseudoxanthomonas,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Pseudoxanthomonas \n",x="Treatment",y=" Relative abundance ")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Opitutus genus
Opitutus_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Opitutus",]
Opitutus_genus_ab<-data.frame(Opitutus=t(Opitutus_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Opitutus_genus_ab)<-c('Opitutus','Treat')
Opitutus_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Opitutus genus.tiff', units="in", width=6, height=5, res=300)

ggplot(Opitutus_genus_ab,aes(x=Treat,y=Opitutus,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values =treats_colors)+labs(title="Opitutus \n",x="Treatment",y=" Relative abundance ")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Staphylococcus genus
# Singulisphaera genus
Singulisphaera_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Singulisphaera",]
Singulisphaera_genus_ab<-data.frame(Singulisphaera=t(Singulisphaera_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Singulisphaera_genus_ab)<-c('Singulisphaera','Treat')
Singulisphaera_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Singulisphaera genus.tiff', units="in", width=6, height=5, res=300)

ggplot(Singulisphaera_genus_ab,aes(x=Treat,y=Singulisphaera,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Singulisphaera \n",x="Treatment",y=" Relative abundance ")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Lysobacter genus
Lysobacter_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Lysobacter",]
Lysobacter_genus_ab<-data.frame(Lysobacter=t(Lysobacter_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Lysobacter_genus_ab)<-c('Lysobacter','Treat')
Lysobacter_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Lysobacter genus.tiff', units="in", width=6, height=5, res=300)

ggplot(Lysobacter_genus_ab,aes(x=Treat,y=Lysobacter,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Lysobacter \n",x="Treatment",y=" Relative abundance ")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()


# Gemmatimonas genus
Gemmatimonas_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Gemmatimonas",]
Gemmatimonas_genus_ab<-data.frame(Gemmatimonas=t(Gemmatimonas_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Gemmatimonas_genus_ab)<-c('Gemmatimonas','Treat')
Gemmatimonas_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Gemmatimonas genus.tiff', units="in", width=6, height=5, res=300)

ggplot(Gemmatimonas_genus_ab,aes(x=Treat,y=Gemmatimonas,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values =treats_colors)+labs(title="Gemmatimonas \n",x="Treatment",y=" Relative abundance ")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()


# Gemmata genus
Gemmata_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Gemmata",]
Gemmata_genus_ab<-data.frame(Gemmata=t(Gemmata_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Gemmata_genus_ab)<-c('Gemmata','Treat')
Gemmata_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Gemmata genus.tiff', units="in", width=6, height=5, res=300)

ggplot(Gemmata_genus_ab,aes(x=Treat,y=Gemmata,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values =treats_colors)+labs(title="Gemmata \n",x="Treatment",y=" Relative abundance ")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()


# Agromyces genus
Agromyces_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Agromyces",]
Agromyces_genus_ab<-data.frame(Agromyces=t(Agromyces_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Agromyces_genus_ab)<-c('Agromyces','Treat')
Agromyces_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Agromyces genus.tiff', units="in", width=6, height=5, res=300)

ggplot(Agromyces_genus_ab,aes(x=Treat,y=Agromyces,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Agromyces \n",x="Treatment",y=" Relative abundance ")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()


# Nocardioides genus
Nocardioides_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Nocardioides",]
Nocardioides_genus_ab<-data.frame(Nocardioides=t(Nocardioides_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Nocardioides_genus_ab)<-c('Nocardioides','Treat')
Nocardioides_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Nocardioides genus.tiff', units="in", width=6, height=5, res=300)

ggplot(Nocardioides_genus_ab,aes(x=Treat,y=Nocardioides,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Nocardioides \n",x="Treatment",y=" Relative abundance ")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()


# Sphingomonas genus
Sphingomonas_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Sphingomonas",]
Sphingomonas_genus_ab<-data.frame(Sphingomonas=t(Sphingomonas_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Sphingomonas_genus_ab)<-c('Sphingomonas','Treat')
Sphingomonas_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Sphingomonas genus.tiff', units="in", width=6, height=5, res=300)

ggplot(Sphingomonas_genus_ab,aes(x=Treat,y=Sphingomonas,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Sphingomonas \n",x="Treatment",y=" Relative abundance ")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()


# Dyella genus
Dyella_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Dyella",]
Dyella_genus_ab<-data.frame(Dyella=t(Dyella_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Dyella_genus_ab)<-c('Dyella','Treat')
Dyella_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Dyella genus.tiff', units="in", width=6, height=5, res=300)

ggplot(Dyella_genus_ab,aes(x=Treat,y=Dyella,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Dyella \n",x="Treatment",y=" Relative abundance ")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Phenylobacterium genus
Phenylobacterium_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Phenylobacterium",]
Phenylobacterium_genus_ab<-data.frame(Phenylobacterium=t(Phenylobacterium_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Phenylobacterium_genus_ab)<-c('Phenylobacterium','Treat')
Phenylobacterium_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Phenylobacterium genus.tiff', units="in", width=6, height=5, res=300)


ggplot(Phenylobacterium_genus_ab,aes(x=Treat,y=Phenylobacterium,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Phenylobacterium \n",x="Treatment",y=" Relative abundance ")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Flavisolibacter genus
Flavisolibacter_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Flavisolibacter",]
Flavisolibacter_genus_ab<-data.frame(Flavisolibacter=t(Flavisolibacter_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Flavisolibacter_genus_ab)<-c('Flavisolibacter','Treat')
Flavisolibacter_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Flavisolibacter genus.tiff', units="in", width=6, height=5, res=300)


ggplot(Flavisolibacter_genus_ab,aes(x=Treat,y=Flavisolibacter,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Flavisolibacter \n",x="Treatment",y=" Relative abundance ")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()


# Streptomyces genus
Streptomyces_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Streptomyces",]
Streptomyces_genus_ab<-data.frame(Streptomyces=t(Streptomyces_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Streptomyces_genus_ab)<-c('Streptomyces','Treat')
Streptomyces_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Streptomyces genus.tiff', units="in", width=6, height=5, res=300)

ggplot(Streptomyces_genus_ab,aes(x=Treat,y=Streptomyces,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Streptomyces \n",x="Treatment",y=" Relative abundance ")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()


# Bradyrhizobium genus
Bradyrhizobium_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Bradyrhizobium",]
Bradyrhizobium_genus_ab<-data.frame(Bradyrhizobium=t(Bradyrhizobium_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Bradyrhizobium_genus_ab)<-c('Bradyrhizobium','Treat')
Bradyrhizobium_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Bradyrhizobium genus.tiff', units="in", width=6, height=5, res=300)


ggplot(Bradyrhizobium_genus_ab,aes(x=Treat,y=Bradyrhizobium,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Bradyrhizobium \n",x="Treatment",y=" Relative abundance ")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Massilia genus
# Pasteuria genus
Pasteuria_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Pasteuria",]
Pasteuria_genus_ab<-data.frame(Pasteuria=t(Pasteuria_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Pasteuria_genus_ab)<-c('Pasteuria','Treat')
Pasteuria_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Pasteuria genus.tiff', units="in", width=6, height=5, res=300)


ggplot(Pasteuria_genus_ab,aes(x=Treat,y=Pasteuria,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Pasteuria \n",x="Treatment",y=" Relative abundance ")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Ohtaekwangia genus
Ohtaekwangia_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Ohtaekwangia",]
Ohtaekwangia_genus_ab<-data.frame(Ohtaekwangia=t(Ohtaekwangia_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Ohtaekwangia_genus_ab)<-c('Ohtaekwangia','Treat')
Ohtaekwangia_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Ohtaekwangia genus.tiff', units="in", width=6, height=5, res=300)

ggplot(Ohtaekwangia_genus_ab,aes(x=Treat,y=Ohtaekwangia,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Ohtaekwangia \n",x="Treatment",y=" Relative abundance ")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()



# Mucilaginibacter genus
Mucilaginibacter_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Mucilaginibacter",]
Mucilaginibacter_genus_ab<-data.frame(Mucilaginibacter=t(Mucilaginibacter_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Mucilaginibacter_genus_ab)<-c('Mucilaginibacter','Treat')
Mucilaginibacter_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Mucilaginibacter genus.tiff', units="in", width=6, height=5, res=300)


ggplot(Mucilaginibacter_genus_ab,aes(x=Treat,y=Mucilaginibacter,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Mucilaginibacter \n",x="Treatment",y=" Relative abundance ")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Mycobacterium genus
Mycobacterium_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Mycobacterium",]
Mycobacterium_genus_ab<-data.frame(Mycobacterium=t(Mycobacterium_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Mycobacterium_genus_ab)<-c('Mycobacterium','Treat')
Mycobacterium_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Mycobacterium genus.tiff', units="in", width=6, height=5, res=300)

ggplot(Mycobacterium_genus_ab,aes(x=Treat,y=Mycobacterium,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Mycobacterium \n",x="Treatment",y=" Relative abundance ")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Pseudomonas genus

Pseudomonas_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Pseudomonas",]
Pseudomonas_genus_ab<-data.frame(Pseudomonas=t(Pseudomonas_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Pseudomonas_genus_ab)<-c('Pseudomonas','Treat')
Pseudomonas_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Pseudomonas genus.tiff', units="in", width=6, height=5, res=300)
ggplot(Pseudomonas_genus_ab,aes(x=Treat,y=Pseudomonas,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Pseudomonas \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Amaricoccus genus

Amaricoccus_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Amaricoccus",]
Amaricoccus_genus_ab<-data.frame(Amaricoccus=t(Amaricoccus_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Amaricoccus_genus_ab)<-c('Amaricoccus','Treat')
Amaricoccus_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Amaricoccus genus.tiff', units="in", width=6, height=5, res=300)
ggplot(Amaricoccus_genus_ab,aes(x=Treat,y=Amaricoccus,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Amaricoccus \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()


# Pantoea genus
Pantoea_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Pantoea",]
Pantoea_genus_ab<-data.frame(Pantoea=t(Pantoea_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Pantoea_genus_ab)<-c('Pantoea','Treat')
Pantoea_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Pantoea genus.tiff', units="in", width=6, height=5, res=300)
ggplot(Pantoea_genus_ab,aes(x=Treat,y=Pantoea,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Pantoea \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Glycomyces genus
Glycomyces_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Glycomyces",]
Glycomyces_genus_ab<-data.frame(Glycomyces=t(Glycomyces_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Glycomyces_genus_ab)<-c('Glycomyces','Treat')
Glycomyces_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Glycomyces genus.tiff', units="in", width=6, height=5, res=300)
ggplot(Glycomyces_genus_ab,aes(x=Treat,y=Glycomyces,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Glycomyces \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()


# Burkholderia genus
Burkholderia_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Burkholderia",]
Burkholderia_genus_ab<-data.frame(Burkholderia=t(Burkholderia_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Burkholderia_genus_ab)<-c('Burkholderia','Treat')
Burkholderia_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Burkholderia genus.tiff', units="in", width=6, height=5, res=300)
ggplot(Burkholderia_genus_ab,aes(x=Treat,y=Burkholderia,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Burkholderia \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Gp3 genus
Gp3_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Gp3",]
Gp3_genus_ab<-data.frame(Gp3=t(Gp3_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Gp3_genus_ab)<-c('Gp3','Treat')
Gp3_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Gp3 genus.tiff', units="in", width=6, height=5, res=300)
ggplot(Gp3_genus_ab,aes(x=Treat,y=Gp3,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Gp3 \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Gp6 genus
Gp6_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Gp6",]
Gp6_genus_ab<-data.frame(Gp6=t(Gp6_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Gp6_genus_ab)<-c('Gp6','Treat')
Gp6_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Gp6 genus.tiff', units="in", width=6, height=5, res=300)
ggplot(Gp6_genus_ab,aes(x=Treat,y=Gp6,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Gp6 \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Gp16 genus
Gp16_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Gp16",]
Gp16_genus_ab<-data.frame(Gp16=t(Gp16_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Gp16_genus_ab)<-c('Gp16','Treat')
Gp16_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Gp16 genus.tiff', units="in", width=6, height=5, res=300)
ggplot(Gp16_genus_ab,aes(x=Treat,y=Gp16,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Gp16 \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()


# Gp1 genus
Gp1_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Gp1",]
Gp1_genus_ab<-data.frame(Gp1=t(Gp1_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Gp1_genus_ab)<-c('Gp1','Treat')
Gp1_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Gp1 genus.tiff', units="in", width=6, height=5, res=300)
ggplot(Gp1_genus_ab,aes(x=Treat,y=Gp1,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Gp1 \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()


# Gp2 genus
Gp2_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Gp2",]
Gp2_genus_ab<-data.frame(Gp2=t(Gp2_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Gp2_genus_ab)<-c('Gp2','Treat')
Gp2_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Gp2 genus.tiff', units="in", width=6, height=5, res=300)
ggplot(Gp2_genus_ab,aes(x=Treat,y=Gp2,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Gp2 \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Gp4 genus
Gp4_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Gp4",]
Gp4_genus_ab<-data.frame(Gp4=t(Gp4_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Gp4_genus_ab)<-c('Gp4','Treat')
Gp4_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Gp4 genus.tiff', units="in", width=6, height=5, res=300)
ggplot(Gp4_genus_ab,aes(x=Treat,y=Gp4,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Gp4 \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()


# TM7_genus_incertae_sedis genus
TM7_genus_incertae_sedis_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="TM7_genus_incertae_sedis",]
TM7_genus_incertae_sedis_genus_ab<-data.frame(TM7_genus_incertae_sedis=t(TM7_genus_incertae_sedis_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(TM7_genus_incertae_sedis_genus_ab)<-c('TM7_genus_incertae_sedis','Treat')
TM7_genus_incertae_sedis_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of TM7_genus_incertae_sedis genus.tiff', units="in", width=6, height=5, res=300)
ggplot(TM7_genus_incertae_sedis_genus_ab,aes(x=Treat,y=TM7_genus_incertae_sedis,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="TM7_genus_incertae_sedis \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Kribbella genus
Kribbella_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Kribbella",]
Kribbella_genus_ab<-data.frame(Kribbella=t(Kribbella_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Kribbella_genus_ab)<-c('Kribbella','Treat')
Kribbella_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Kribbella genus.tiff', units="in", width=6, height=5, res=300)
ggplot(Kribbella_genus_ab,aes(x=Treat,y=Kribbella,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Kribbella \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Planctomyces genus
Planctomyces_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Planctomyces",]
Planctomyces_genus_ab<-data.frame(Planctomyces=t(Planctomyces_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Planctomyces_genus_ab)<-c('Planctomyces','Treat')
Planctomyces_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Planctomyces genus.tiff', units="in", width=6, height=5, res=300)
ggplot(Planctomyces_genus_ab,aes(x=Treat,y=Planctomyces,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Planctomyces \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Ktedonobacter genus
Ktedonobacter_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Ktedonobacter",]
Ktedonobacter_genus_ab<-data.frame(Ktedonobacter=t(Ktedonobacter_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Ktedonobacter_genus_ab)<-c('Ktedonobacter','Treat')
Ktedonobacter_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Ktedonobacter genus.tiff', units="in", width=6, height=5, res=300)
ggplot(Ktedonobacter_genus_ab,aes(x=Treat,y=Ktedonobacter,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Ktedonobacter \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()


# Spartobacteria_genera_incertae_sedis genus
Spartobacteria_genera_incertae_sedis_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Spartobacteria_genera_incertae_sedis",]
Spartobacteria_genera_incertae_sedis_genus_ab<-data.frame(Spartobacteria_genera_incertae_sedis=t(Spartobacteria_genera_incertae_sedis_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Spartobacteria_genera_incertae_sedis_genus_ab)<-c('Spartobacteria_genera_incertae_sedis','Treat')
Spartobacteria_genera_incertae_sedis_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Spartobacteria_genera_incertae_sedis genus.tiff', units="in", width=6, height=5, res=300)
ggplot(Spartobacteria_genera_incertae_sedis_genus_ab,aes(x=Treat,y=Spartobacteria_genera_incertae_sedis,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Spartobacteria_genera_incertae_sedis \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Enterobacter genus
Enterobacter_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Enterobacter",]
Enterobacter_genus_ab<-data.frame(Enterobacter=t(Enterobacter_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Enterobacter_genus_ab)<-c('Enterobacter','Treat')
Enterobacter_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Enterobacter genus.tiff', units="in", width=6, height=5, res=300)
ggplot(Enterobacter_genus_ab,aes(x=Treat,y=Enterobacter,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Enterobacter \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()


# Nocardioides genus
Nocardioides_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Nocardioides",]
Nocardioides_genus_ab<-data.frame(Nocardioides=t(Nocardioides_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Nocardioides_genus_ab)<-c('Nocardioides','Treat')
Nocardioides_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Nocardioides genus.tiff', units="in", width=6, height=5, res=300)
ggplot(Nocardioides_genus_ab,aes(x=Treat,y=Nocardioides,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Nocardioides \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Novosphingobium genus
Novosphingobium_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Novosphingobium",]
Novosphingobium_genus_ab<-data.frame(Novosphingobium=t(Novosphingobium_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Novosphingobium_genus_ab)<-c('Novosphingobium','Treat')
Novosphingobium_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Novosphingobium genus.tiff', units="in", width=6, height=5, res=300)
ggplot(Novosphingobium_genus_ab,aes(x=Treat,y=Novosphingobium,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Novosphingobium \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

#Bacillus genus
Bacillus_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Bacillus",]
Bacillus_genus_ab<-data.frame(Bacillus=t(Bacillus_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Bacillus_genus_ab)<-c('Bacillus','Treat')
Bacillus_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Bacillus genus.tiff', units="in", width=6, height=5, res=300)
ggplot(Bacillus_genus_ab,aes(x=Treat,y=Bacillus,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Bacillus \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

#Rhizobium genus
Rhizobium_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Rhizobium",]
Rhizobium_genus_ab<-data.frame(Rhizobium=t(Rhizobium_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Rhizobium_genus_ab)<-c('Rhizobium','Treat')
Rhizobium_genus_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Rhizobium genus.tiff', units="in", width=6, height=5, res=300)
ggplot(Rhizobium_genus_ab,aes(x=Treat,y=Rhizobium,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Rhizobium \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()


# Stenotrophomonas genus

Stenotrophomonas_genus_ab<-otu_table(r_genera_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_genera_rarefied_cultivar_phyloseq))$Genus=="Stenotrophomonas",]
Stenotrophomonas_genus_ab<-data.frame(Stenotrophomonas=t(Stenotrophomonas_genus_ab),Treat=sample_data(r_genera_rarefied_cultivar_phyloseq)$Treat)
colnames(Stenotrophomonas_genus_ab)<-c('Stenotrophomonas','Treat')
Stenotrophomonas_genus_ab


tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Stenotrophomonas genus.tiff', units="in", width=6, height=5, res=300)
ggplot(Stenotrophomonas_genus_ab,aes(x=Treat,y=Stenotrophomonas,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Stenotrophomona \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()


#Bacteroidetes phylum
Bacteroidetes_phylum_ab<-otu_table(r_phylum_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_phylum_rarefied_cultivar_phyloseq))$Phylum=="Bacteroidetes",]
Bacteroidetes_phylum_ab<-data.frame(Bacteroidetes=t(Bacteroidetes_phylum_ab),Treat=sample_data(r_phylum_rarefied_cultivar_phyloseq)$Treat)
colnames(Bacteroidetes_phylum_ab)<-c('Bacteroidetes','Treat')
Bacteroidetes_phylum_ab


tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Bacteroidetes phylum.tiff', units="in", width=6, height=5, res=300)
ggplot(Bacteroidetes_phylum_ab,aes(x=Treat,y=Bacteroidetes,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Bacteroidetes \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()


# Actinobacteria phylum
Actinobacteria_phylum_ab<-otu_table(r_phylum_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_phylum_rarefied_cultivar_phyloseq))$Phylum=="Actinobacteria",]
Actinobacteria_phylum_ab<-data.frame(Actinobacteria=t(Actinobacteria_phylum_ab),Treat=sample_data(r_phylum_rarefied_cultivar_phyloseq)$Treat)
colnames(Actinobacteria_phylum_ab)<-c('Actinobacteria','Treat')
Actinobacteria_phylum_ab

tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Actinobacteria phylum.tiff', units="in", width=6, height=5, res=300)
ggplot(Actinobacteria_phylum_ab,aes(x=Treat,y=Actinobacteria,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Actinobacteria \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Acidobacteria phylum
Acidobacteria_phylum_ab<-otu_table(r_phylum_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_phylum_rarefied_cultivar_phyloseq))$Phylum=="Acidobacteria",]
Acidobacteria_phylum_ab<-data.frame(Acidobacteria=t(Acidobacteria_phylum_ab),Treat=sample_data(r_phylum_rarefied_cultivar_phyloseq)$Treat)
colnames(Acidobacteria_phylum_ab)<-c('Acidobacteria','Treat')
Acidobacteria_phylum_ab


tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Acidobacteria phylum.tiff', units="in", width=6, height=5, res=300)
ggplot(Acidobacteria_phylum_ab,aes(x=Treat,y=Acidobacteria,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Acidobacteria \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()


# Proteobacteria phylum
Proteobacteria_phylum_ab<-otu_table(r_phylum_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_phylum_rarefied_cultivar_phyloseq))$Phylum=="Proteobacteria",]
Proteobacteria_phylum_ab<-data.frame(Acidobacteria=t(Proteobacteria_phylum_ab),Treat=sample_data(r_phylum_rarefied_cultivar_phyloseq)$Treat)
colnames(Proteobacteria_phylum_ab)<-c('Proteobacteria','Treat')
Proteobacteria_phylum_ab


tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Proteobacteria phylum.tiff', units="in", width=6, height=5, res=300)
ggplot(Proteobacteria_phylum_ab,aes(x=Treat,y=Proteobacteria,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Proteobacteria \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Verrucomicrobia phylum
Verrucomicrobia_phylum_ab<-otu_table(r_phylum_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_phylum_rarefied_cultivar_phyloseq))$Phylum=="Verrucomicrobia",]
Verrucomicrobia_phylum_ab<-data.frame(Acidobacteria=t(Verrucomicrobia_phylum_ab),Treat=sample_data(r_phylum_rarefied_cultivar_phyloseq)$Treat)
colnames(Verrucomicrobia_phylum_ab)<-c('Verrucomicrobia','Treat')
Verrucomicrobia_phylum_ab


tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Verrucomicrobia phylum.tiff', units="in", width=6, height=5, res=300)
ggplot(Verrucomicrobia_phylum_ab,aes(x=Treat,y=Verrucomicrobia,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Verrucomicrobia \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Firmicutes phylum
Firmicutes_phylum_ab<-otu_table(r_phylum_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_phylum_rarefied_cultivar_phyloseq))$Phylum=="Firmicutes",]
Firmicutes_phylum_ab<-data.frame(Firmicutes=t(Firmicutes_phylum_ab),Treat=sample_data(r_phylum_rarefied_cultivar_phyloseq)$Treat)
colnames(Firmicutes_phylum_ab)<-c('Firmicutes','Treat')
Firmicutes_phylum_ab


tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Firmicutes phylum.tiff', units="in", width=6, height=5, res=300)
ggplot(Firmicutes_phylum_ab,aes(x=Treat,y=Firmicutes,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Firmicutes \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Alphaproteobacteria class

Alphaproteobacteria_class_ab<-otu_table(r_class_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_class_rarefied_cultivar_phyloseq))$Class=="Alphaproteobacteria",]
Alphaproteobacteria_class_ab<-data.frame(Alphaproteobacteria=t(Alphaproteobacteria_class_ab),Treat=sample_data(r_class_rarefied_cultivar_phyloseq)$Treat)
colnames(Alphaproteobacteria_class_ab)<-c('Alphaproteobacteria','Treat')
Alphaproteobacteria_class_ab


tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Alphaproteobacteria class.tiff', units="in", width=6, height=5, res=300)
ggplot(Alphaproteobacteria_class_ab,aes(x=Treat,y=Alphaproteobacteria,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Alphaproteobacteria \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Gammaproteobacteria class

Gammaproteobacteria_class_ab<-otu_table(r_class_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_class_rarefied_cultivar_phyloseq))$Class=="Gammaproteobacteria",]
Gammaproteobacteria_class_ab<-data.frame(Gammaproteobacteria=t(Gammaproteobacteria_class_ab),Treat=sample_data(r_class_rarefied_cultivar_phyloseq)$Treat)
colnames(Gammaproteobacteria_class_ab)<-c('Gammaproteobacteria','Treat')
Gammaproteobacteria_class_ab


tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Gammaproteobacteria class.tiff', units="in", width=6, height=5, res=300)
ggplot(Gammaproteobacteria_class_ab,aes(x=Treat,y=Gammaproteobacteria,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Gammaproteobacteria \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Betaproteobacteria class

Betaproteobacteria_class_ab<-otu_table(r_class_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_class_rarefied_cultivar_phyloseq))$Class=="Betaproteobacteria",]
Betaproteobacteria_class_ab<-data.frame(Betaproteobacteria=t(Betaproteobacteria_class_ab),Treat=sample_data(r_class_rarefied_cultivar_phyloseq)$Treat)
colnames(Betaproteobacteria_class_ab)<-c('Betaproteobacteria','Treat')
Betaproteobacteria_class_ab


tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Betaproteobacteria class.tiff', units="in", width=6, height=5, res=300)
ggplot(Betaproteobacteria_class_ab,aes(x=Treat,y=Betaproteobacteria,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Betaproteobacteria \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Enterobacteriaceae family

Enterobacteriaceae_family_ab<-otu_table(r_family_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_family_rarefied_cultivar_phyloseq))$Family=="Enterobacteriaceae",]

Enterobacteriaceae_family_ab<-data.frame(Enterobacteriaceae=t(Enterobacteriaceae_family_ab),Treat=sample_data(r_family_rarefied_cultivar_phyloseq)$Treat)
colnames(Enterobacteriaceae_family_ab)<-c('Enterobacteriaceae','Treat')
Enterobacteriaceae_family_ab


tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Enterobacteriaceae family.tiff', units="in", width=6, height=5, res=300)
ggplot(Enterobacteriaceae_family_ab,aes(x=Treat,y=Enterobacteriaceae,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Enterobacteriaceae \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))

dev.off()

# Chitinophagaceae family

Chitinophagaceae_family_ab<-otu_table(r_family_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_family_rarefied_cultivar_phyloseq))$Family=="Chitinophagaceae",]

Chitinophagaceae_family_ab<-data.frame(Chitinophagaceae=t(Chitinophagaceae_family_ab),Treat=sample_data(r_family_rarefied_cultivar_phyloseq)$Treat)
colnames(Chitinophagaceae_family_ab)<-c('Chitinophagaceae','Treat')
Chitinophagaceae_family_ab


tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Chitinophagaceae family.tiff', units="in", width=6, height=5, res=300)
ggplot(Chitinophagaceae_family_ab,aes(x=Treat,y=Chitinophagaceae,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Chitinophagaceae \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))

dev.off()

# Enterobacteriaceae family

Enterobacteriaceae_family_ab<-otu_table(r_family_rarefied_cultivar_phyloseq)[data.frame(tax_table(r_family_rarefied_cultivar_phyloseq))$Family=="Enterobacteriaceae",]

Enterobacteriaceae_family_ab<-data.frame(Enterobacteriaceae=t(Enterobacteriaceae_family_ab),Treat=sample_data(r_family_rarefied_cultivar_phyloseq)$Treat)
colnames(Enterobacteriaceae_family_ab)<-c('Enterobacteriaceae','Treat')
Enterobacteriaceae_family_ab


tiff('/Users/fangliu/Documents/2016_cultivar_project/R_analysis_up/R_output_image/Relative abundance of Enterobacteriaceae family.tiff', units="in", width=6, height=5, res=300)
ggplot(Enterobacteriaceae_family_ab,aes(x=Treat,y=Enterobacteriaceae,fill=Treat))+geom_boxplot(size=0.5)+
geom_jitter(shape=16, position=position_jitter(0.2),aes(fill=Treat))+scale_fill_manual(values = treats_colors)+labs(title="Enterobacteriaceae \n",x="Treatment",y="Relative abundance")+theme(panel.background = element_rect(fill='white',color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.text.x = element_text(family="Arial",size=13,angle=60,hjust = 1),legend.position = "none",
    axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = -0.5),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))

dev.off()


```



