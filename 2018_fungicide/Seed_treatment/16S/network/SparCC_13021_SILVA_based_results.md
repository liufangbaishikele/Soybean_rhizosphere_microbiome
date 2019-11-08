---
title: "Seed_16S_SparCC"
author: "Fang Liu"
date: "8/16/2019"
output: html_document
editor_options: 
  chunk_output_type: console
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## set up working environment and load libraries

```{r}
setwd("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC")
library(phyloseq)
library(igraph)
library(plyr)
```

## Read into a whole .shared and cons.taxonomy file to define label factor 
## This step is used to generate labels for each genus, which could be used to merge with four groups of vertex dataframe based on Genus name.

```{r}
# Read in the original strigolactone output OTU table and taxonomy table

Seed_13021<-import_mothur(mothur_shared_file = "/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.subsample.0.03.pick.shared",mothur_constaxonomy_file = "/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.subsample.0.03.pick.0.03.cons.taxonomy")
Seed_13021 # 21355 taxa and 150 samples
colnames(tax_table(Seed_13021))<-c('Kingdom','Phylum','Class','Order','Family','Genus')

# Read in the meta data
Seed_meta<-read.csv('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/Seed_meta.csv',row.names = 1,header = TRUE)
head(Seed_meta)

Seed_meta_phyloseq<-sample_data(Seed_meta)

#>>>>> Seed_phyloseq_up >>>>> with sample data integrated

Seed_13021_up<-merge_phyloseq(Seed_13021,Seed_meta_phyloseq)
Seed_13021_up
sample_names(Seed_13021_up)[1:5]
Seed_13021_up
identical(sample_names(Seed_13021),sample_names(Seed_13021_up))

# creat label vector

meta<-data.frame(label=1:length(taxa_names(Seed_13021_up)),OtuID=taxa_names(Seed_13021_up),Phylum=tax_table(Seed_13021_up)[,2])
head(meta)
```

#  Subset OTUs to creat input for SparCC analysis
---------------------------------------------------

```{r}
#Subset the samples based on treatment

#---Soil---

Soil<-subset_samples(Seed_13021_up,Treat=="Soil")
Soil
#Extract top 1000 OTUs across 8 samples
Soil_up<-prune_taxa(labels(sort(taxa_sums(Soil),decreasing = TRUE)[1:500]),Soil)
Soil_up #1000 taxa and 15 samples
min(taxa_sums(Soil_up)) # 67
Soil_up_OTU_count<-data.frame(otu_table(Soil_up))
#write.csv(Soil_up_OTU_count,file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Soil_up_OTU_count.csv")
table(Soil_up_OTU_count==0)[2]/sum(table(Soil_up_OTU_count==0)) #0.05493333

#---Bulk_1wk---

Bulk_1wk<-subset_samples(Seed_13021_up,Treat=="Bulk_1wk")
Bulk_1wk
#Extract top 1000 OTUs across 8 samples
Bulk_1wk_up<-prune_taxa(labels(sort(taxa_sums(Bulk_1wk),decreasing = TRUE)[1:500]),Bulk_1wk)
Bulk_1wk_up #1000 taxa and 15 samples
min(taxa_sums(Bulk_1wk_up)) # 70
Bulk_1wk_up_OTU_count<-data.frame(otu_table(Bulk_1wk_up))
#write.csv(Bulk_1wk_up_OTU_count,file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Bulk_1wk_up_OTU_count.csv")
table(Bulk_1wk_up_OTU_count==0)[2]/sum(table(Bulk_1wk_up_OTU_count==0)) #0.05546667

#---Bulk_3wk---

Bulk_3wk<-subset_samples(Seed_13021_up,Treat=="Bulk_3wk")
Bulk_3wk
#Extract top 1000 OTUs across 8 samples
Bulk_3wk_up<-prune_taxa(labels(sort(taxa_sums(Bulk_3wk),decreasing = TRUE)[1:500]),Bulk_3wk)
Bulk_3wk_up #1000 taxa and 15 samples
min(taxa_sums(Bulk_3wk_up)) # 70
Bulk_3wk_up_OTU_count<-data.frame(otu_table(Bulk_3wk_up))
#write.csv(Bulk_3wk_up_OTU_count,file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Bulk_3wk_up_OTU_count.csv")
table(Bulk_3wk_up_OTU_count==0)[2]/sum(table(Bulk_3wk_up_OTU_count==0)) #0.06713333

#---Bulk_4wk---

Bulk_4wk<-subset_samples(Seed_13021_up,Treat=="Bulk_4wk")
Bulk_4wk
#Extract top 1000 OTUs across 8 samples
Bulk_4wk_up<-prune_taxa(labels(sort(taxa_sums(Bulk_4wk),decreasing = TRUE)[1:500]),Bulk_4wk)
Bulk_4wk_up #1000 taxa and 15 samples
min(taxa_sums(Bulk_4wk_up)) # 49
Bulk_4wk_up_OTU_count<-data.frame(otu_table(Bulk_4wk_up))
#write.csv(Bulk_4wk_up_OTU_count,file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Bulk_4wk_up_OTU_count.csv")
table(Bulk_4wk_up_OTU_count==0)[2]/sum(table(Bulk_4wk_up_OTU_count==0)) #0.06713333

#---Rhi_1wk---

Rhi_1wk<-subset_samples(Seed_13021_up,Treat=="Rhi_1wk")
Rhi_1wk
#Extract top 1000 OTUs across 8 samples
Rhi_1wk_up<-prune_taxa(labels(sort(taxa_sums(Rhi_1wk),decreasing = TRUE)[1:500]),Rhi_1wk)
Rhi_1wk_up #1000 taxa and 15 samples
min(taxa_sums(Rhi_1wk_up)) # 31
Rhi_1wk_up_OTU_count<-data.frame(otu_table(Rhi_1wk_up))
#write.csv(Rhi_1wk_up_OTU_count,file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Rhi_1wk_up_OTU_count.csv")
table(Rhi_1wk_up_OTU_count==0)[2]/sum(table(Rhi_1wk_up_OTU_count==0)) #0.1121

# --- Rhi_3wk ---

Rhi_3wk<-subset_samples(Seed_13021_up,Treat=="Rhi_3wk")
Rhi_3wk
#Extract top 1000 OTUs across 8 samples
Rhi_3wk_up<-prune_taxa(labels(sort(taxa_sums(Rhi_3wk),decreasing = TRUE)[1:500]),Rhi_3wk)
Rhi_3wk_up #1000 taxa and 15 samples
min(taxa_sums(Rhi_3wk_up)) # 71
Rhi_3wk_up_OTU_count<-data.frame(otu_table(Rhi_3wk_up))
#write.csv(Rhi_3wk_up_OTU_count,file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Rhi_3wk_up_OTU_count.csv")
table(Rhi_3wk_up_OTU_count==0)[2]/sum(table(Rhi_3wk_up_OTU_count==0)) #0.044

# --- Rhi_4wk ---

Rhi_4wk<-subset_samples(Seed_13021_up,Treat=="Rhi_4wk")
Rhi_4wk
#Extract top 1000 OTUs across 15 samples
Rhi_4wk_up<-prune_taxa(labels(sort(taxa_sums(Rhi_4wk),decreasing = TRUE)[1:500]),Rhi_4wk)
Rhi_4wk_up #1000 taxa and 15 samples
min(taxa_sums(Rhi_4wk_up)) # 43
Rhi_4wk_up_OTU_count<-data.frame(otu_table(Rhi_4wk_up))
#write.csv(Rhi_4wk_up_OTU_count,file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Rhi_4wk_up_OTU_count.csv")
table(Rhi_4wk_up_OTU_count==0)[2]/sum(table(Rhi_4wk_up_OTU_count==0)) #0.05013333

# --- Endo_1wk ---

Endo_1wk<-subset_samples(Seed_13021_up,Treat=="Endo_1wk")
Endo_1wk
#Extract top 1000 OTUs across 15 samples
Endo_1wk_up<-prune_taxa(labels(sort(taxa_sums(Endo_1wk),decreasing = TRUE)[1:500]),Endo_1wk)
Endo_1wk_up #1000 taxa and 15 samples
min(taxa_sums(Endo_1wk_up)) # 35
Endo_1wk_up_OTU_count<-data.frame(otu_table(Endo_1wk_up))
write.csv(Endo_1wk_up_OTU_count,file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Endo_1wk_up_OTU_count.csv")
table(Endo_1wk_up_OTU_count==0)[2]/sum(table(Endo_1wk_up_OTU_count==0)) #0.1417333

# --- Endo_3wk ---

Endo_3wk<-subset_samples(Seed_13021_up,Treat=="Endo_3wk")
Endo_3wk
#Extract top 1000 OTUs across 15 samples
Endo_3wk_up<-prune_taxa(labels(sort(taxa_sums(Endo_3wk),decreasing = TRUE)[1:500]),Endo_3wk)
Endo_3wk_up #1000 taxa and 15 samples
min(taxa_sums(Endo_3wk_up)) # 18
Endo_3wk_up_OTU_count<-data.frame(otu_table(Endo_3wk_up))
#write.csv(Endo_3wk_up_OTU_count,file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Endo_3wk_up_OTU_count.csv")
table(Endo_3wk_up_OTU_count==0)[2]/sum(table(Endo_3wk_up_OTU_count==0)) #0.3006

# --- Endo_4wk ---

Endo_4wk<-subset_samples(Seed_13021_up,Treat=="Endo_4wk")
Endo_4wk
#Extract top 1000 OTUs across 15 samples
Endo_4wk_up<-prune_taxa(labels(sort(taxa_sums(Endo_4wk),decreasing = TRUE)[1:500]),Endo_4wk)
Endo_4wk_up #1000 taxa and 15 samples
min(taxa_sums(Endo_4wk_up)) # 20
Endo_4wk_up_OTU_count<-data.frame(otu_table(Endo_4wk_up))
#write.csv(Endo_4wk_up_OTU_count,file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Endo_4wk_up_OTU_count.csv")
table(Endo_4wk_up_OTU_count==0)[2]/sum(table(Endo_4wk_up_OTU_count==0)) #0.1938667

#--- Seed ---

SEED_all<-import_mothur(mothur_shared_file = "/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared",mothur_constaxonomy_file = "/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.subsample.0.03.pick.0.03.cons.taxonomy")

SEED_all<-merge_phyloseq(SEED_all,Seed_meta_phyloseq)

SEED<-subset_samples(SEED_all,Compartment=="SEED")
SEED_up<-filter_taxa(SEED,function(x) sum(x)>1,prune = TRUE)
SEED_up # 841 taxa and 15 samples
SEED_up<-prune_samples(sample_names(SEED_up)[-4],SEED_up) # Here removed 

sample_sums(SEED_up)
#Extract top 1000 OTUs across 15 samples
SEED_top50<-prune_taxa(labels(sort(taxa_sums(SEED_up),decreasing = TRUE)[1:50]),SEED_up)
SEED_top50 #50 taxa and 15 samples
min(taxa_sums(SEED_top50)) # 38
SEED_top50_OTU_count<-data.frame(otu_table(SEED_top50))
#write.csv(SEED_top50_OTU_count,file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/SEED_top50_OTU_count.csv")
table(SEED_top50_OTU_count==0)[2]/sum(table(SEED_top50_OTU_count==0)) #0.2871429 
```

## Generate SparCC network with output corr and pvalue matrixed

# >>>>> Soil_up >>>>>>>

```{r}
## import SparCC and p_value matrix into R 
Soil_up_SparCC_cor<-read.csv("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Soil_up_OTU_cor_sparcc.csv",header = TRUE,row.names = 1)
Soil_up_SparCC_cor<-as.matrix(Soil_up_SparCC_cor)
dim(Soil_up_SparCC_cor)
Soil_up_SparCC_cor[1:5,1:5]

# Read in corresponding p_value dataframe
Soil_up_SparCC_pvalue<-read.csv("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Soil_up_OTU_two_sided_pvalue.csv",header = TRUE,row.names = 1)
Soil_up_SparCC_pvalue<-as.matrix(Soil_up_SparCC_pvalue)
dim(Soil_up_SparCC_pvalue)
Soil_up_SparCC_pvalue[1:5,1:5]
```


## Generate vertices information and add relative abundance information to vertices and color information


```{r}
r_Soil_up_phyloseq<-transform_sample_counts(Soil_up,function(x) x/sum(x))

Soil_up_vertex<-data.frame(OtuID=taxa_names(r_Soil_up_phyloseq),tax_table(r_Soil_up_phyloseq),size=taxa_sums(r_Soil_up_phyloseq)/length(sample_sums(r_Soil_up_phyloseq))) # Here size are the average relative abundance of each OTU across all samples within same treatment
dim(Soil_up_vertex)
Soil_up_vertex[1:5,]
Soil_up_vertex<-Soil_up_vertex[,c(1,3,7,8)]
head(Soil_up_vertex)

Bacteria_phylum_vs_color_checklist<-data.frame(phylum=c('Proteobacteria','Actinobacteria','Bacteria_unclassified','Acidobacteria','Nitrospirae','Chloroflexi','Candidatus_Saccharibacteria','Firmicutes','Bacteroidetes','Verrucomicrobia','candidate_division_WPS-1','Gemmatimonadetes','Planctomycetes','BRC1','Latescibacteria','Chlamydiae','Others','Patescibacteria','Rokubacteria','Cyanobacteria'),colors=c('#ff9082','#98d66d','#fcef5d','#c45cc0','#45b1f9','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#bc8c38','#bcba6d','#91749e','#b2798d','#343837','#b23301','#235931'))
phylum_colors<-as.character(Bacteria_phylum_vs_color_checklist$colors)
phylum_labels<-as.character(Bacteria_phylum_vs_color_checklist$phylum)

#Others = ("candidate_division_WPS-2","Deinococcus-Thermus","Elusimicrobia","Hydrogenedentes","Microgenomates","Parcubacteria","Spirochaetes","SR1","Synergistetes","Tenericutes","Armatimonadetes")

tiff("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Phylum_colour_pie.tiff",units = 'in',width=10,height = 8,res=300,compression = 'lzw')
pie(rep (1,20),col=phylum_colors,labels=paste(phylum_labels,phylum_colors,sep=":"),radius = 1.0)
dev.off()

Soil_up_vertex$color<-factor(Soil_up_vertex$Phylum,levels=phylum_labels,labels=phylum_colors)
dim(Soil_up_vertex)
Soil_up_vertex[1:10,]
# I want to add label variable to Soil_up_vertex
Soil_up_vertex<-join (Soil_up_vertex,meta,by="OtuID") 

# Originally, I used merge function, but it turned out to be problematic, because after merge it will The rows are by default lexicographically sorted on the common columns. To resolve this problem, I changed merge function to joint function.
dim(Soil_up_vertex)
Soil_up_vertex[1:5,]
```


## Creat links and have a look at edge adn vertex information

```{r}
# read SparCC correlation matrix and p_value matrix to igraph
Soil_up_net<-graph_from_adjacency_matrix(Soil_up_SparCC_cor,weighted = TRUE)
Soil_up_net  #IGRAPH 8cfb97f DNW- 500 250000 -- 
summary(V(Soil_up_net)$name==as.character(Soil_up_vertex$OtuID)) 
V(Soil_up_net)$name[1:5]
E(Soil_up_net)$weight[1:5] # here the weight are just the SparCC coefficience

# Add vertex information to the network
V(Soil_up_net)$vertex_color<-as.character(Soil_up_vertex$color)
V(Soil_up_net)$size<-Soil_up_vertex$size
V(Soil_up_net)$Phylum<-as.character(Soil_up_vertex$Phylum)
V(Soil_up_net)$label<-as.character(Soil_up_vertex$label)

# Add edge information (p_value) to the network
E(Soil_up_net)$p_value<-Soil_up_SparCC_pvalue
E(Soil_up_net)$p_value[1:10]

## Add edge width to the network
E(Soil_up_net)$width<-E(Soil_up_net)$weight*10
E(Soil_up_net)$width[1:20]
```


------------------------------------------------------------------------------------------------------------------------------------------
In order to simplify the network, including remove of loop and mutual edges. If using simplify function on Soil_up_net, this will loose p_value and edge color information (which indicate if the interaction are negarive or positive)

To solve this problem, I will write out the edge and vertices in dataframe format and read into these two dataframe and modify the edge information globally before use ``graph_from_data_frame`` to create network object.
-----------------------------------------------------------------------------------------------------------------------------------------


## write out edges and vertices in dataframe format, and edit on the network

```{r}

# -- generate edge dataframe from whole network---

Soil_up_edges<-igraph::as_data_frame(Soil_up_net,what = "edges")
#Soil_up_edges
dim(Soil_up_edges) # 250000 rows and 5 columns
head(Soil_up_edges)
# Here the weight indicate the strength of interaction (SparCC value), p_value is the significance of interactions nad width is formulated based on weight * 40.

# -- generate vertices dataframe from whole network---

Soil_up_vertices<-igraph::as_data_frame(Soil_up_net,what = "vertices")
dim(Soil_up_vertices)
head(Soil_up_vertices)

# remove edges that from and to the same nodes(selfloops)

Soil_up_edges<-Soil_up_edges[which(Soil_up_edges$from!=Soil_up_edges$to),]
dim(Soil_up_edges) # 249500 rows and 5 columns
head(Soil_up_edges)

# Creat a new network using the above vertices and edge data frame
Soil_up_net1 <- graph_from_data_frame(d=Soil_up_edges, vertices=Soil_up_vertices, directed=T)#here if I used directed=FALSE, I will lost my p_value information as well as weight information. 
Soil_up_net1 # IGRAPH 526a39b DNW- 800 639200 --
E(Soil_up_net1)$weight[1:10]

# convert directed network to undirected network

Soil_up_net2<-as.undirected(Soil_up_net1 ,mode = 'mutual',edge.attr.comb = "mean") # as a result, weight, p_value and width of edges are all averaged for mutual edges.
Soil_up_net2 # IGRAPH 1bce66c UNW- 500 124750 - -
E(Soil_up_net2)$weight[1:10]
E(Soil_up_net2)$p_value[1:10]
E(Soil_up_net2)$width[1:10]


## simplify the network from Soil_up_net2

#1)  remove non_significant edges with alpha=0.05

summary(E(Soil_up_net2)$p_value<0.05) #TRUE 74374 and FALSE 530076   

Soil_up_net3<-delete.edges(Soil_up_net2,which(E(Soil_up_net2)$p_value>=0.05))
Soil_up_net3 # IGRAPH f25c247 UNW- 800 24383 -- 
summary(degree(Soil_up_net3)==0) # all of the nodes has connection to others
summary(E(Soil_up_net3)$weight>0) # 12668 positive interactions, 11715 nagative interactions

#2) remove non-significant edges with alpha=0.001

summary(E(Soil_up_net2)$p_value<0.001) # TRUE 14033, FALSE 590417
Soil_up_net4<-delete.edges(Soil_up_net2,which(E(Soil_up_net2)$p_value>=0.001))
Soil_up_net4 # IGRAPH b2b3bc7 UNW- 800 2481 --
Soil_up_net4<-delete.vertices(Soil_up_net4,degree(Soil_up_net4)==0) # IGRAPH 281b75b UNW- 789 2481 --

#3) further simplify the network based on degree
barplot(sort(degree(Soil_up_net4)),main = "Soil_up_net4_barplot")
summary(degree(Soil_up_net4)>15) # 36 nodes with degree larger than 100
Soil_up_net5<-delete.vertices(Soil_up_net4,which(degree(Soil_up_net4)<=5))
Soil_up_net5# IGRAPH ce4509d UNW- 36 536 --
summary(E(Soil_up_net5)$weight>0) # 261 positive and 275 negative

# top 50

Soil_up_net6<-delete.vertices(Soil_up_net4,labels(sort(degree(Soil_up_net4),FALSE)[c(1:(length(V(Soil_up_net4)$name)-50))])) # IGRAPH 56796e5 UNW- 50 937 --
Soil_up_net6<-delete.vertices(Soil_up_net6,degree(Soil_up_net6)==0) #IGRAPH cead187 UNW- 43 71 -- 
summary(E(Soil_up_net6)$weight>0) # logical : 46(negative)     25(positive)

# top100
Soil_up_net7<-delete.vertices(Soil_up_net4,labels(sort(degree(Soil_up_net4),FALSE)[c(1:(length(V(Soil_up_net4)$name)-100))]))
Soil_up_net7<-delete.vertices(Soil_up_net7,degree(Soil_up_net7)==0) #IGRAPH cc57451 UNW- 96 200 -
summary(E(Soil_up_net7)$weight>0) # 1273 (negative)    1228 (positive) 
```


## Plot network

```{r}

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Soil_up_net4.tiff', units="in", width=10, height=8, res=300)

plot(Soil_up_net4,edge.arrow.size=0,edge.width=abs(E(Soil_up_net4)$width/4),vertex.color=V(Soil_up_net4)$vertex_color,vertex.frame.color="#555555",vertex.label=NA,edge.color=c("#ff8989","#9fe899")[1+(E(Soil_up_net4)$weight>0)],vertex.size=degree(Soil_up_net4)/10,layout=layout_with_fr (Soil_up_net4),main="Soil_up_net4") # green edges means positive correlation and red edges mean negative correlation
dev.off()

#plot(Soil_up_net5,edge.arrow.size=0,edge.width=abs(E(Soil_up_net5)$width/3),vertex.color=V(Soil_up_net5)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Soil_up_net5)$label,edge.color=c("#ff8989","#9fe899")[1+(E(Soil_up_net5)$weight>0)],vertex.label.cex=.8,vertex.size=degree(Soil_up_net5)/2,layout=layout_with_dh (Soil_up_net5),main="Soil_up_net5")

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Soil_up_net6.tiff', units="in", width=10, height=8, res=300)
plot(Soil_up_net6,edge.arrow.size=0,edge.width=abs(E(Soil_up_net6)$width/3),vertex.color=V(Soil_up_net6)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Soil_up_net6)$label,vertex.label.color="black",vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,edge.color=c("#ff8989","#9fe899")[1+(E(Soil_up_net6)$weight>0)],vertex.size=degree(Soil_up_net6)/3,layout=layout_with_dh (Soil_up_net6),main="Soil_up_net6")
dev.off()

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Soil_up_net7.tiff', units="in", width=10, height=8, res=300)
plot(Soil_up_net7,edge.arrow.size=0,edge.width=abs(E(Soil_up_net7)$width/3),vertex.color=V(Soil_up_net7)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Soil_up_net7)$label,vertex.label.color="black",vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,edge.color=c("#ff8989","#9fe899")[1+(E(Soil_up_net7)$weight>0)],vertex.size=degree(Soil_up_net7)/3,layout=layout_with_dh (Soil_up_net7),main="Soil_up_net7")
dev.off()
```

#  >>>>> Bulk_1wk_up >>>>>>>

```{r}
## import SparCC and p_value matrix into R 
Bulk_1wk_up_SparCC_cor<-read.csv("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Bulk_1wk_up_OTU_cor_sparcc.csv",header = TRUE,row.names = 1)
Bulk_1wk_up_SparCC_cor<-as.matrix(Bulk_1wk_up_SparCC_cor)
dim(Bulk_1wk_up_SparCC_cor)
Bulk_1wk_up_SparCC_cor[1:5,1:5]

# Read in corresponding p_value dataframe
Bulk_1wk_up_SparCC_pvalue<-read.csv("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Bulk_1wk_up_OTU_two_sided_pvalue.csv",header = TRUE,row.names = 1)
Bulk_1wk_up_SparCC_pvalue<-as.matrix(Bulk_1wk_up_SparCC_pvalue)
dim(Bulk_1wk_up_SparCC_pvalue)
Bulk_1wk_up_SparCC_pvalue[1:5,1:5]
```


## Generate vertices information and add relative abundance information to vertices and color information


```{r}
r_Bulk_1wk_up_phyloseq<-transform_sample_counts(Bulk_1wk_up,function(x) x/sum(x))

Bulk_1wk_up_vertex<-data.frame(OtuID=taxa_names(r_Bulk_1wk_up_phyloseq),tax_table(r_Bulk_1wk_up_phyloseq),size=taxa_sums(r_Bulk_1wk_up_phyloseq)/length(sample_sums(r_Bulk_1wk_up_phyloseq))) # Here size are the average relative abundance of each OTU across all samples within same treatment
dim(Bulk_1wk_up_vertex)
Bulk_1wk_up_vertex[1:5,]
Bulk_1wk_up_vertex<-Bulk_1wk_up_vertex[,c(1,3,7,8)]
head(Bulk_1wk_up_vertex)

Bacteria_phylum_vs_color_checklist<-data.frame(phylum=c('Proteobacteria','Actinobacteria','Bacteria_unclassified','Acidobacteria','Nitrospirae','Chloroflexi','Candidatus_Saccharibacteria','Firmicutes','Bacteroidetes','Verrucomicrobia','candidate_division_WPS-1','Gemmatimonadetes','Planctomycetes','BRC1','Latescibacteria','Chlamydiae','Others','Patescibacteria','Rokubacteria','Cyanobacteria'),colors=c('#ff9082','#98d66d','#fcef5d','#c45cc0','#45b1f9','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#bc8c38','#bcba6d','#91749e','#b2798d','#343837','#b23301','#235931'))
phylum_colors<-as.character(Bacteria_phylum_vs_color_checklist$colors)
phylum_labels<-as.character(Bacteria_phylum_vs_color_checklist$phylum)


tiff("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Phylum_colour_pie.tiff",units = 'in',width=10,height = 8,res=300,compression = 'lzw')
pie(rep (1,17),col=phylum_colors,labels=paste(phylum_labels,phylum_colors,sep=":"),radius = 1.0)
dev.off()

Bulk_1wk_up_vertex$color<-factor(Bulk_1wk_up_vertex$Phylum,levels=phylum_labels,labels=phylum_colors)
dim(Bulk_1wk_up_vertex)
Bulk_1wk_up_vertex[1:10,]
# I want to add label variable to Bulk_1wk_up_vertex
Bulk_1wk_up_vertex<-join (Bulk_1wk_up_vertex,meta,by="OtuID") 

# Originally, I used merge function, but it turned out to be problematic, because after merge it will The rows are by default lexicographically sorted on the common columns. To resolve this problem, I changed merge function to joint function.
dim(Bulk_1wk_up_vertex)
Bulk_1wk_up_vertex[1:5,]
```


## Creat links and have a look at edge adn vertex information

```{r}
# read SparCC correlation matrix and p_value matrix to igraph
Bulk_1wk_up_net<-graph_from_adjacency_matrix(Bulk_1wk_up_SparCC_cor,weighted = TRUE)
Bulk_1wk_up_net  #IGRAPH 8cfb97f DNW- 500 250000 -- 
summary(V(Bulk_1wk_up_net)$name==as.character(Bulk_1wk_up_vertex$OtuID)) 
V(Bulk_1wk_up_net)$name[1:5]
E(Bulk_1wk_up_net)$weight[1:5] # here the weight are just the SparCC coefficience

# Add vertex information to the network
V(Bulk_1wk_up_net)$vertex_color<-as.character(Bulk_1wk_up_vertex$color)
V(Bulk_1wk_up_net)$size<-Bulk_1wk_up_vertex$size
V(Bulk_1wk_up_net)$Phylum<-as.character(Bulk_1wk_up_vertex$Phylum)
V(Bulk_1wk_up_net)$label<-as.character(Bulk_1wk_up_vertex$label)

# Add edge information (p_value) to the network
E(Bulk_1wk_up_net)$p_value<-Bulk_1wk_up_SparCC_pvalue
E(Bulk_1wk_up_net)$p_value[1:10]

## Add edge width to the network
E(Bulk_1wk_up_net)$width<-E(Bulk_1wk_up_net)$weight*10
E(Bulk_1wk_up_net)$width[1:20]
```


------------------------------------------------------------------------------------------------------------------------------------------
In order to simplify the network, including remove of loop and mutual edges. If using simplify function on Bulk_1wk_up_net, this will loose p_value and edge color information (which indicate if the interaction are negarive or positive)

To solve this problem, I will write out the edge and vertices in dataframe format and read into these two dataframe and modify the edge information globally before use ``graph_from_data_frame`` to create network object.
-----------------------------------------------------------------------------------------------------------------------------------------


## write out edges and vertices in dataframe format, and edit on the network

```{r}

# -- generate edge dataframe from whole network---

Bulk_1wk_up_edges<-igraph::as_data_frame(Bulk_1wk_up_net,what = "edges")
#Bulk_1wk_up_edges
dim(Bulk_1wk_up_edges) # 250000 rows and 5 columns
head(Bulk_1wk_up_edges)
# Here the weight indicate the strength of interaction (SparCC value), p_value is the significance of interactions nad width is formulated based on weight * 40.

# -- generate vertices dataframe from whole network---

Bulk_1wk_up_vertices<-igraph::as_data_frame(Bulk_1wk_up_net,what = "vertices")
dim(Bulk_1wk_up_vertices)
head(Bulk_1wk_up_vertices)

# remove edges that from and to the same nodes(selfloops)

Bulk_1wk_up_edges<-Bulk_1wk_up_edges[which(Bulk_1wk_up_edges$from!=Bulk_1wk_up_edges$to),]
dim(Bulk_1wk_up_edges) # 1208900 rows and 5 columns
head(Bulk_1wk_up_edges)

# Creat a new network using the above vertices and edge data frame
Bulk_1wk_up_net1 <- graph_from_data_frame(d=Bulk_1wk_up_edges, vertices=Bulk_1wk_up_vertices, directed=T)#here if I used directed=FALSE, I will lost my p_value information as well as weight information. 
Bulk_1wk_up_net1 # IGRAPH 526a39b DNW- 800 639200 --
E(Bulk_1wk_up_net1)$weight[1:10]

# convert directed network to undirected network

Bulk_1wk_up_net2<-as.undirected(Bulk_1wk_up_net1 ,mode = 'mutual',edge.attr.comb = "mean") # as a result, weight, p_value and width of edges are all averaged for mutual edges.
Bulk_1wk_up_net2 # IGRAPH 1bce66c UNW- 500 124750 - -
E(Bulk_1wk_up_net2)$weight[1:10]
E(Bulk_1wk_up_net2)$p_value[1:10]
E(Bulk_1wk_up_net2)$width[1:10]


## simplify the network from Bulk_1wk_up_net2

#1)  remove non_significant edges with alpha=0.05

summary(E(Bulk_1wk_up_net2)$p_value<0.05) #TRUE 74374 and FALSE 530076   

Bulk_1wk_up_net3<-delete.edges(Bulk_1wk_up_net2,which(E(Bulk_1wk_up_net2)$p_value>=0.05))
Bulk_1wk_up_net3 # IGRAPH f25c247 UNW- 800 24383 -- 
summary(degree(Bulk_1wk_up_net3)==0) # all of the nodes has connection to others
summary(E(Bulk_1wk_up_net3)$weight>0) # 12668 positive interactions, 11715 nagative interactions

#2) remove non-significant edges with alpha=0.001

summary(E(Bulk_1wk_up_net2)$p_value<0.001) # TRUE 14033, FALSE 590417
Bulk_1wk_up_net4<-delete.edges(Bulk_1wk_up_net2,which(E(Bulk_1wk_up_net2)$p_value>=0.001))
Bulk_1wk_up_net4 # IGRAPH b2b3bc7 UNW- 800 2481 --
Bulk_1wk_up_net4<-delete.vertices(Bulk_1wk_up_net4,degree(Bulk_1wk_up_net4)==0) # IGRAPH 281b75b UNW- 789 2481 --

#3) further simplify the network based on degree
barplot(sort(degree(Bulk_1wk_up_net4)),main = "Bulk_1wk_up_net4_barplot")
summary(degree(Bulk_1wk_up_net4)>15) # 36 nodes with degree larger than 100
Bulk_1wk_up_net5<-delete.vertices(Bulk_1wk_up_net4,which(degree(Bulk_1wk_up_net4)<=5))
Bulk_1wk_up_net5# IGRAPH ce4509d UNW- 36 536 --
summary(E(Bulk_1wk_up_net5)$weight>0) # 261 positive and 275 negative

# top 50

Bulk_1wk_up_net6<-delete.vertices(Bulk_1wk_up_net4,labels(sort(degree(Bulk_1wk_up_net4),FALSE)[c(1:(length(V(Bulk_1wk_up_net4)$name)-50))])) # IGRAPH 56796e5 UNW- 50 937 --
Bulk_1wk_up_net6<-delete.vertices(Bulk_1wk_up_net6,degree(Bulk_1wk_up_net6)==0) #IGRAPH cead187 UNW- 43 71 -- 
summary(E(Bulk_1wk_up_net6)$weight>0) # logical : 46(negative)     25(positive)

# top100
Bulk_1wk_up_net7<-delete.vertices(Bulk_1wk_up_net4,labels(sort(degree(Bulk_1wk_up_net4),FALSE)[c(1:(length(V(Bulk_1wk_up_net4)$name)-100))]))
Bulk_1wk_up_net7<-delete.vertices(Bulk_1wk_up_net7,degree(Bulk_1wk_up_net7)==0) #IGRAPH cc57451 UNW- 96 200 -
summary(E(Bulk_1wk_up_net7)$weight>0) # 1273 (negative)    1228 (positive) 
```


## Plot network

```{r}

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Bulk_1wk_up_net4.tiff', units="in", width=10, height=8, res=300)

plot(Bulk_1wk_up_net4,edge.arrow.size=0,edge.width=abs(E(Bulk_1wk_up_net4)$width/4),vertex.color=V(Bulk_1wk_up_net4)$vertex_color,vertex.frame.color="#555555",vertex.label=NA,edge.color=c("#ff8989","#9fe899")[1+(E(Bulk_1wk_up_net4)$weight>0)],vertex.size=degree(Bulk_1wk_up_net4)/10,layout=layout_with_fr (Bulk_1wk_up_net4),main="Bulk_1wk_up_net4") # green edges means positive correlation and red edges mean negative correlation
dev.off()

#plot(Bulk_1wk_up_net5,edge.arrow.size=0,edge.width=abs(E(Bulk_1wk_up_net5)$width/3),vertex.color=V(Bulk_1wk_up_net5)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Bulk_1wk_up_net5)$label,edge.color=c("#ff8989","#9fe899")[1+(E(Bulk_1wk_up_net5)$weight>0)],vertex.label.cex=.8,vertex.size=degree(Bulk_1wk_up_net5)/2,layout=layout_with_dh (Bulk_1wk_up_net5),main="Bulk_1wk_up_net5")

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Bulk_1wk_up_net6.tiff', units="in", width=10, height=8, res=300)
plot(Bulk_1wk_up_net6,edge.arrow.size=0,edge.width=abs(E(Bulk_1wk_up_net6)$width/3),vertex.color=V(Bulk_1wk_up_net6)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Bulk_1wk_up_net6)$label,vertex.label.color="black",vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,edge.color=c("#ff8989","#9fe899")[1+(E(Bulk_1wk_up_net6)$weight>0)],vertex.size=degree(Bulk_1wk_up_net6)/3,layout=layout_with_dh (Bulk_1wk_up_net6),main="Bulk_1wk_up_net6")
dev.off()

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Bulk_1wk_up_net7.tiff', units="in", width=10, height=8, res=300)
plot(Bulk_1wk_up_net7,edge.arrow.size=0,edge.width=abs(E(Bulk_1wk_up_net7)$width/3),vertex.color=V(Bulk_1wk_up_net7)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Bulk_1wk_up_net7)$label,vertex.label.color="black",vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,edge.color=c("#ff8989","#9fe899")[1+(E(Bulk_1wk_up_net7)$weight>0)],vertex.size=degree(Bulk_1wk_up_net7)/3,layout=layout_with_dh (Bulk_1wk_up_net7),main="Bulk_1wk_up_net7")
dev.off()
```

#  >>>>> Bulk_3wk_up >>>>>>>

```{r}
## import SparCC and p_value matrix into R 
Bulk_3wk_up_SparCC_cor<-read.csv("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Bulk_3wk_up_OTU_cor_sparcc.csv",header = TRUE,row.names = 1)
Bulk_3wk_up_SparCC_cor<-as.matrix(Bulk_3wk_up_SparCC_cor)
dim(Bulk_3wk_up_SparCC_cor)
Bulk_3wk_up_SparCC_cor[1:5,1:5]

# Read in corresponding p_value dataframe
Bulk_3wk_up_SparCC_pvalue<-read.csv("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Bulk_3wk_up_OTU_two_sided_pvalue.csv",header = TRUE,row.names = 1)
Bulk_3wk_up_SparCC_pvalue<-as.matrix(Bulk_3wk_up_SparCC_pvalue)
dim(Bulk_3wk_up_SparCC_pvalue)
Bulk_3wk_up_SparCC_pvalue[1:5,1:5]
```


## Generate vertices information and add relative abundance information to vertices and color information


```{r}
r_Bulk_3wk_up_phyloseq<-transform_sample_counts(Bulk_3wk_up,function(x) x/sum(x))

Bulk_3wk_up_vertex<-data.frame(OtuID=taxa_names(r_Bulk_3wk_up_phyloseq),tax_table(r_Bulk_3wk_up_phyloseq),size=taxa_sums(r_Bulk_3wk_up_phyloseq)/length(sample_sums(r_Bulk_3wk_up_phyloseq))) # Here size are the average relative abundance of each OTU across all samples within same treatment
dim(Bulk_3wk_up_vertex)
Bulk_3wk_up_vertex[1:5,]
Bulk_3wk_up_vertex<-Bulk_3wk_up_vertex[,c(1,3,7,8)]
head(Bulk_3wk_up_vertex)

Bacteria_phylum_vs_color_checklist<-data.frame(phylum=c('Proteobacteria','Actinobacteria','Bacteria_unclassified','Acidobacteria','Nitrospirae','Chloroflexi','Candidatus_Saccharibacteria','Firmicutes','Bacteroidetes','Verrucomicrobia','candidate_division_WPS-1','Gemmatimonadetes','Planctomycetes','BRC1','Latescibacteria','Chlamydiae','Others','Patescibacteria','Rokubacteria','Cyanobacteria'),colors=c('#ff9082','#98d66d','#fcef5d','#c45cc0','#45b1f9','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#bc8c38','#bcba6d','#91749e','#b2798d','#343837','#b23301','#235931'))
phylum_colors<-as.character(Bacteria_phylum_vs_color_checklist$colors)
phylum_labels<-as.character(Bacteria_phylum_vs_color_checklist$phylum)


tiff("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Phylum_colour_pie.tiff",units = 'in',width=10,height = 8,res=300,compression = 'lzw')
pie(rep (1,17),col=phylum_colors,labels=paste(phylum_labels,phylum_colors,sep=":"),radius = 1.0)
dev.off()

Bulk_3wk_up_vertex$color<-factor(Bulk_3wk_up_vertex$Phylum,levels=phylum_labels,labels=phylum_colors)
dim(Bulk_3wk_up_vertex)
Bulk_3wk_up_vertex[1:10,]
# I want to add label variable to Bulk_3wk_up_vertex
Bulk_3wk_up_vertex<-join (Bulk_3wk_up_vertex,meta,by="OtuID") 

# Originally, I used merge function, but it turned out to be problematic, because after merge it will The rows are by default lexicographically sorted on the common columns. To resolve this problem, I changed merge function to joint function.
dim(Bulk_3wk_up_vertex)
Bulk_3wk_up_vertex[1:5,]
```


## Creat links and have a look at edge adn vertex information

```{r}
# read SparCC correlation matrix and p_value matrix to igraph
Bulk_3wk_up_net<-graph_from_adjacency_matrix(Bulk_3wk_up_SparCC_cor,weighted = TRUE)
Bulk_3wk_up_net  #IGRAPH 8cfb97f DNW- 500 250000 -- 
summary(V(Bulk_3wk_up_net)$name==as.character(Bulk_3wk_up_vertex$OtuID)) 
V(Bulk_3wk_up_net)$name[1:5]
E(Bulk_3wk_up_net)$weight[1:5] # here the weight are just the SparCC coefficience

# Add vertex information to the network
V(Bulk_3wk_up_net)$vertex_color<-as.character(Bulk_3wk_up_vertex$color)
V(Bulk_3wk_up_net)$size<-Bulk_3wk_up_vertex$size
V(Bulk_3wk_up_net)$Phylum<-as.character(Bulk_3wk_up_vertex$Phylum)
V(Bulk_3wk_up_net)$label<-as.character(Bulk_3wk_up_vertex$label)

# Add edge information (p_value) to the network
E(Bulk_3wk_up_net)$p_value<-Bulk_3wk_up_SparCC_pvalue
E(Bulk_3wk_up_net)$p_value[1:10]

## Add edge width to the network
E(Bulk_3wk_up_net)$width<-E(Bulk_3wk_up_net)$weight*10
E(Bulk_3wk_up_net)$width[1:20]
```


------------------------------------------------------------------------------------------------------------------------------------------
In order to simplify the network, including remove of loop and mutual edges. If using simplify function on Bulk_3wk_up_net, this will loose p_value and edge color information (which indicate if the interaction are negarive or positive)

To solve this problem, I will write out the edge and vertices in dataframe format and read into these two dataframe and modify the edge information globally before use ``graph_from_data_frame`` to create network object.
-----------------------------------------------------------------------------------------------------------------------------------------


## write out edges and vertices in dataframe format, and edit on the network

```{r}

# -- generate edge dataframe from whole network---

Bulk_3wk_up_edges<-igraph::as_data_frame(Bulk_3wk_up_net,what = "edges")
#Bulk_3wk_up_edges
dim(Bulk_3wk_up_edges) # 250000 rows and 5 columns
head(Bulk_3wk_up_edges)
# Here the weight indicate the strength of interaction (SparCC value), p_value is the significance of interactions nad width is formulated based on weight * 40.

# -- generate vertices dataframe from whole network---

Bulk_3wk_up_vertices<-igraph::as_data_frame(Bulk_3wk_up_net,what = "vertices")
dim(Bulk_3wk_up_vertices)
head(Bulk_3wk_up_vertices)

# remove edges that from and to the same nodes(selfloops)

Bulk_3wk_up_edges<-Bulk_3wk_up_edges[which(Bulk_3wk_up_edges$from!=Bulk_3wk_up_edges$to),]
dim(Bulk_3wk_up_edges) # 1208900 rows and 5 columns
head(Bulk_3wk_up_edges)

# Creat a new network using the above vertices and edge data frame
Bulk_3wk_up_net1 <- graph_from_data_frame(d=Bulk_3wk_up_edges, vertices=Bulk_3wk_up_vertices, directed=T)#here if I used directed=FALSE, I will lost my p_value information as well as weight information. 
Bulk_3wk_up_net1 # IGRAPH 526a39b DNW- 800 639200 --
E(Bulk_3wk_up_net1)$weight[1:10]

# convert directed network to undirected network

Bulk_3wk_up_net2<-as.undirected(Bulk_3wk_up_net1 ,mode = 'mutual',edge.attr.comb = "mean") # as a result, weight, p_value and width of edges are all averaged for mutual edges.
Bulk_3wk_up_net2 # IGRAPH 1bce66c UNW- 500 124750 - -
E(Bulk_3wk_up_net2)$weight[1:10]
E(Bulk_3wk_up_net2)$p_value[1:10]
E(Bulk_3wk_up_net2)$width[1:10]


## simplify the network from Bulk_3wk_up_net2

#1)  remove non_significant edges with alpha=0.05

summary(E(Bulk_3wk_up_net2)$p_value<0.05) #TRUE 74374 and FALSE 530076   

Bulk_3wk_up_net3<-delete.edges(Bulk_3wk_up_net2,which(E(Bulk_3wk_up_net2)$p_value>=0.05))
Bulk_3wk_up_net3 # IGRAPH f25c247 UNW- 800 24383 -- 
summary(degree(Bulk_3wk_up_net3)==0) # all of the nodes has connection to others
summary(E(Bulk_3wk_up_net3)$weight>0) # 12668 positive interactions, 11715 nagative interactions

#2) remove non-significant edges with alpha=0.001

summary(E(Bulk_3wk_up_net2)$p_value<0.001) # TRUE 14033, FALSE 590417
Bulk_3wk_up_net4<-delete.edges(Bulk_3wk_up_net2,which(E(Bulk_3wk_up_net2)$p_value>=0.001))
Bulk_3wk_up_net4 # IGRAPH b2b3bc7 UNW- 800 2481 --
Bulk_3wk_up_net4<-delete.vertices(Bulk_3wk_up_net4,degree(Bulk_3wk_up_net4)==0) # IGRAPH 281b75b UNW- 789 2481 --

#3) further simplify the network based on degree
barplot(sort(degree(Bulk_3wk_up_net4)),main = "Bulk_3wk_up_net4_barplot")
summary(degree(Bulk_3wk_up_net4)>15) # 36 nodes with degree larger than 100
Bulk_3wk_up_net5<-delete.vertices(Bulk_3wk_up_net4,which(degree(Bulk_3wk_up_net4)<=5))
Bulk_3wk_up_net5# IGRAPH ce4509d UNW- 36 536 --
summary(E(Bulk_3wk_up_net5)$weight>0) # 261 positive and 275 negative

# top 50

Bulk_3wk_up_net6<-delete.vertices(Bulk_3wk_up_net4,labels(sort(degree(Bulk_3wk_up_net4),FALSE)[c(1:(length(V(Bulk_3wk_up_net4)$name)-50))])) # IGRAPH 56796e5 UNW- 50 937 --
Bulk_3wk_up_net6<-delete.vertices(Bulk_3wk_up_net6,degree(Bulk_3wk_up_net6)==0) #IGRAPH cead187 UNW- 43 71 -- 
summary(E(Bulk_3wk_up_net6)$weight>0) # logical : 46(negative)     25(positive)

# top100
Bulk_3wk_up_net7<-delete.vertices(Bulk_3wk_up_net4,labels(sort(degree(Bulk_3wk_up_net4),FALSE)[c(1:(length(V(Bulk_3wk_up_net4)$name)-100))]))
Bulk_3wk_up_net7<-delete.vertices(Bulk_3wk_up_net7,degree(Bulk_3wk_up_net7)==0) #IGRAPH cc57451 UNW- 96 200 -
summary(E(Bulk_3wk_up_net7)$weight>0) # 1273 (negative)    1228 (positive) 
```


## Plot network

```{r}

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Bulk_3wk_up_net4.tiff', units="in", width=10, height=8, res=300)

plot(Bulk_3wk_up_net4,edge.arrow.size=0,edge.width=abs(E(Bulk_3wk_up_net4)$width/4),vertex.color=V(Bulk_3wk_up_net4)$vertex_color,vertex.frame.color="#555555",vertex.label=NA,edge.color=c("#ff8989","#9fe899")[1+(E(Bulk_3wk_up_net4)$weight>0)],vertex.size=degree(Bulk_3wk_up_net4)/10,layout=layout_with_fr (Bulk_3wk_up_net4),main="Bulk_3wk_up_net4") # green edges means positive correlation and red edges mean negative correlation
dev.off()

#plot(Bulk_3wk_up_net5,edge.arrow.size=0,edge.width=abs(E(Bulk_3wk_up_net5)$width/3),vertex.color=V(Bulk_3wk_up_net5)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Bulk_3wk_up_net5)$label,edge.color=c("#ff8989","#9fe899")[1+(E(Bulk_3wk_up_net5)$weight>0)],vertex.label.cex=.8,vertex.size=degree(Bulk_3wk_up_net5)/2,layout=layout_with_dh (Bulk_3wk_up_net5),main="Bulk_3wk_up_net5")

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Bulk_3wk_up_net6.tiff', units="in", width=10, height=8, res=300)
plot(Bulk_3wk_up_net6,edge.arrow.size=0,edge.width=abs(E(Bulk_3wk_up_net6)$width/3),vertex.color=V(Bulk_3wk_up_net6)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Bulk_3wk_up_net6)$label,vertex.label.color="black",vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,edge.color=c("#ff8989","#9fe899")[1+(E(Bulk_3wk_up_net6)$weight>0)],vertex.size=degree(Bulk_3wk_up_net6)/3,layout=layout_with_dh (Bulk_3wk_up_net6),main="Bulk_3wk_up_net6")
dev.off()

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Bulk_3wk_up_net7.tiff', units="in", width=10, height=8, res=300)
plot(Bulk_3wk_up_net7,edge.arrow.size=0,edge.width=abs(E(Bulk_3wk_up_net7)$width/3),vertex.color=V(Bulk_3wk_up_net7)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Bulk_3wk_up_net7)$label,vertex.label.color="black",vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,edge.color=c("#ff8989","#9fe899")[1+(E(Bulk_3wk_up_net7)$weight>0)],vertex.size=degree(Bulk_3wk_up_net7)/3,layout=layout_with_dh (Bulk_3wk_up_net7),main="Bulk_3wk_up_net7")
dev.off()
```

#  >>>>> Bulk_4wk_up >>>>>>>

```{r}
## import SparCC and p_value matrix into R 
Bulk_4wk_up_SparCC_cor<-read.csv("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Bulk_4wk_up_OTU_cor_sparcc.csv",header = TRUE,row.names = 1)
Bulk_4wk_up_SparCC_cor<-as.matrix(Bulk_4wk_up_SparCC_cor)
dim(Bulk_4wk_up_SparCC_cor)
Bulk_4wk_up_SparCC_cor[1:5,1:5]

# Read in corresponding p_value dataframe
Bulk_4wk_up_SparCC_pvalue<-read.csv("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Bulk_4wk_up_OTU_two_sided_pvalue.csv",header = TRUE,row.names = 1)
Bulk_4wk_up_SparCC_pvalue<-as.matrix(Bulk_4wk_up_SparCC_pvalue)
dim(Bulk_4wk_up_SparCC_pvalue)
Bulk_4wk_up_SparCC_pvalue[1:5,1:5]
```


## Generate vertices information and add relative abundance information to vertices and color information


```{r}
r_Bulk_4wk_up_phyloseq<-transform_sample_counts(Bulk_4wk_up,function(x) x/sum(x))

Bulk_4wk_up_vertex<-data.frame(OtuID=taxa_names(r_Bulk_4wk_up_phyloseq),tax_table(r_Bulk_4wk_up_phyloseq),size=taxa_sums(r_Bulk_4wk_up_phyloseq)/length(sample_sums(r_Bulk_4wk_up_phyloseq))) # Here size are the average relative abundance of each OTU across all samples within same treatment
dim(Bulk_4wk_up_vertex)
Bulk_4wk_up_vertex[1:5,]
Bulk_4wk_up_vertex<-Bulk_4wk_up_vertex[,c(1,3,7,8)]
head(Bulk_4wk_up_vertex)

Bacteria_phylum_vs_color_checklist<-data.frame(phylum=c('Proteobacteria','Actinobacteria','Bacteria_unclassified','Acidobacteria','Nitrospirae','Chloroflexi','Candidatus_Saccharibacteria','Firmicutes','Bacteroidetes','Verrucomicrobia','candidate_division_WPS-1','Gemmatimonadetes','Planctomycetes','BRC1','Latescibacteria','Chlamydiae','Others','Patescibacteria','Rokubacteria','Cyanobacteria'),colors=c('#ff9082','#98d66d','#fcef5d','#c45cc0','#45b1f9','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#bc8c38','#bcba6d','#91749e','#b2798d','#343837','#b23301','#235931'))
phylum_colors<-as.character(Bacteria_phylum_vs_color_checklist$colors)
phylum_labels<-as.character(Bacteria_phylum_vs_color_checklist$phylum)


tiff("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Phylum_colour_pie.tiff",units = 'in',width=10,height = 8,res=300,compression = 'lzw')
pie(rep (1,17),col=phylum_colors,labels=paste(phylum_labels,phylum_colors,sep=":"),radius = 1.0)
dev.off()

Bulk_4wk_up_vertex$color<-factor(Bulk_4wk_up_vertex$Phylum,levels=phylum_labels,labels=phylum_colors)
dim(Bulk_4wk_up_vertex)
Bulk_4wk_up_vertex[1:10,]
# I want to add label variable to Bulk_4wk_up_vertex
Bulk_4wk_up_vertex<-join (Bulk_4wk_up_vertex,meta,by="OtuID") 

# Originally, I used merge function, but it turned out to be problematic, because after merge it will The rows are by default lexicographically sorted on the common columns. To resolve this problem, I changed merge function to joint function.
dim(Bulk_4wk_up_vertex)
Bulk_4wk_up_vertex[1:5,]
```


## Creat links and have a look at edge adn vertex information

```{r}
# read SparCC correlation matrix and p_value matrix to igraph
Bulk_4wk_up_net<-graph_from_adjacency_matrix(Bulk_4wk_up_SparCC_cor,weighted = TRUE)
Bulk_4wk_up_net  #IGRAPH 8cfb97f DNW- 500 250000 -- 
summary(V(Bulk_4wk_up_net)$name==as.character(Bulk_4wk_up_vertex$OtuID)) 
V(Bulk_4wk_up_net)$name[1:5]
E(Bulk_4wk_up_net)$weight[1:5] # here the weight are just the SparCC coefficience

# Add vertex information to the network
V(Bulk_4wk_up_net)$vertex_color<-as.character(Bulk_4wk_up_vertex$color)
V(Bulk_4wk_up_net)$size<-Bulk_4wk_up_vertex$size
V(Bulk_4wk_up_net)$Phylum<-as.character(Bulk_4wk_up_vertex$Phylum)
V(Bulk_4wk_up_net)$label<-as.character(Bulk_4wk_up_vertex$label)

# Add edge information (p_value) to the network
E(Bulk_4wk_up_net)$p_value<-Bulk_4wk_up_SparCC_pvalue
E(Bulk_4wk_up_net)$p_value[1:10]

## Add edge width to the network
E(Bulk_4wk_up_net)$width<-E(Bulk_4wk_up_net)$weight*10
E(Bulk_4wk_up_net)$width[1:20]
```


------------------------------------------------------------------------------------------------------------------------------------------
In order to simplify the network, including remove of loop and mutual edges. If using simplify function on Bulk_4wk_up_net, this will loose p_value and edge color information (which indicate if the interaction are negarive or positive)

To solve this problem, I will write out the edge and vertices in dataframe format and read into these two dataframe and modify the edge information globally before use ``graph_from_data_frame`` to create network object.
-----------------------------------------------------------------------------------------------------------------------------------------


## write out edges and vertices in dataframe format, and edit on the network

```{r}

# -- generate edge dataframe from whole network---

Bulk_4wk_up_edges<-igraph::as_data_frame(Bulk_4wk_up_net,what = "edges")
#Bulk_4wk_up_edges
dim(Bulk_4wk_up_edges) # 250000 rows and 5 columns
head(Bulk_4wk_up_edges)
# Here the weight indicate the strength of interaction (SparCC value), p_value is the significance of interactions nad width is formulated based on weight * 40.

# -- generate vertices dataframe from whole network---

Bulk_4wk_up_vertices<-igraph::as_data_frame(Bulk_4wk_up_net,what = "vertices")
dim(Bulk_4wk_up_vertices)
head(Bulk_4wk_up_vertices)

# remove edges that from and to the same nodes(selfloops)

Bulk_4wk_up_edges<-Bulk_4wk_up_edges[which(Bulk_4wk_up_edges$from!=Bulk_4wk_up_edges$to),]
dim(Bulk_4wk_up_edges) # 1208900 rows and 5 columns
head(Bulk_4wk_up_edges)

# Creat a new network using the above vertices and edge data frame
Bulk_4wk_up_net1 <- graph_from_data_frame(d=Bulk_4wk_up_edges, vertices=Bulk_4wk_up_vertices, directed=T)#here if I used directed=FALSE, I will lost my p_value information as well as weight information. 
Bulk_4wk_up_net1 # IGRAPH 526a39b DNW- 800 639200 --
E(Bulk_4wk_up_net1)$weight[1:10]

# convert directed network to undirected network

Bulk_4wk_up_net2<-as.undirected(Bulk_4wk_up_net1 ,mode = 'mutual',edge.attr.comb = "mean") # as a result, weight, p_value and width of edges are all averaged for mutual edges.
Bulk_4wk_up_net2 # IGRAPH 1bce66c UNW- 500 124750 - -
E(Bulk_4wk_up_net2)$weight[1:10]
E(Bulk_4wk_up_net2)$p_value[1:10]
E(Bulk_4wk_up_net2)$width[1:10]


## simplify the network from Bulk_4wk_up_net2

#1)  remove non_significant edges with alpha=0.05

summary(E(Bulk_4wk_up_net2)$p_value<0.05) #TRUE 74374 and FALSE 530076   

Bulk_4wk_up_net3<-delete.edges(Bulk_4wk_up_net2,which(E(Bulk_4wk_up_net2)$p_value>=0.05))
Bulk_4wk_up_net3 # IGRAPH f25c247 UNW- 800 24383 -- 
summary(degree(Bulk_4wk_up_net3)==0) # all of the nodes has connection to others
summary(E(Bulk_4wk_up_net3)$weight>0) # 12668 positive interactions, 11715 nagative interactions

#2) remove non-significant edges with alpha=0.001

summary(E(Bulk_4wk_up_net2)$p_value<0.001) # TRUE 14033, FALSE 590417
Bulk_4wk_up_net4<-delete.edges(Bulk_4wk_up_net2,which(E(Bulk_4wk_up_net2)$p_value>=0.001))
Bulk_4wk_up_net4 # IGRAPH b2b3bc7 UNW- 800 2481 --
Bulk_4wk_up_net4<-delete.vertices(Bulk_4wk_up_net4,degree(Bulk_4wk_up_net4)==0) # IGRAPH 281b75b UNW- 789 2481 --

#3) further simplify the network based on degree
barplot(sort(degree(Bulk_4wk_up_net4)),main = "Bulk_4wk_up_net4_barplot")
summary(degree(Bulk_4wk_up_net4)>15) # 36 nodes with degree larger than 100
Bulk_4wk_up_net5<-delete.vertices(Bulk_4wk_up_net4,which(degree(Bulk_4wk_up_net4)<=5))
Bulk_4wk_up_net5# IGRAPH ce4509d UNW- 36 536 --
summary(E(Bulk_4wk_up_net5)$weight>0) # 261 positive and 275 negative

# top 50

Bulk_4wk_up_net6<-delete.vertices(Bulk_4wk_up_net4,labels(sort(degree(Bulk_4wk_up_net4),FALSE)[c(1:(length(V(Bulk_4wk_up_net4)$name)-50))])) # IGRAPH 56796e5 UNW- 50 937 --
Bulk_4wk_up_net6<-delete.vertices(Bulk_4wk_up_net6,degree(Bulk_4wk_up_net6)==0) #IGRAPH cead187 UNW- 43 71 -- 
summary(E(Bulk_4wk_up_net6)$weight>0) # logical : 46(negative)     25(positive)

# top100
Bulk_4wk_up_net7<-delete.vertices(Bulk_4wk_up_net4,labels(sort(degree(Bulk_4wk_up_net4),FALSE)[c(1:(length(V(Bulk_4wk_up_net4)$name)-100))]))
Bulk_4wk_up_net7<-delete.vertices(Bulk_4wk_up_net7,degree(Bulk_4wk_up_net7)==0) #IGRAPH cc57451 UNW- 96 200 -
summary(E(Bulk_4wk_up_net7)$weight>0) # 1273 (negative)    1228 (positive) 
```

## Plot network

```{r}

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Bulk_4wk_up_net4.tiff', units="in", width=10, height=8, res=300)

plot(Bulk_4wk_up_net4,edge.arrow.size=0,edge.width=abs(E(Bulk_4wk_up_net4)$width/4),vertex.color=V(Bulk_4wk_up_net4)$vertex_color,vertex.frame.color="#555555",vertex.label=NA,edge.color=c("#ff8989","#9fe899")[1+(E(Bulk_4wk_up_net4)$weight>0)],vertex.size=degree(Bulk_4wk_up_net4)/10,layout=layout_with_fr (Bulk_4wk_up_net4),main="Bulk_4wk_up_net4") # green edges means positive correlation and red edges mean negative correlation
dev.off()

#plot(Bulk_4wk_up_net5,edge.arrow.size=0,edge.width=abs(E(Bulk_4wk_up_net5)$width/3),vertex.color=V(Bulk_4wk_up_net5)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Bulk_4wk_up_net5)$label,edge.color=c("#ff8989","#9fe899")[1+(E(Bulk_4wk_up_net5)$weight>0)],vertex.label.cex=.8,vertex.size=degree(Bulk_4wk_up_net5)/2,layout=layout_with_dh (Bulk_4wk_up_net5),main="Bulk_4wk_up_net5")

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Bulk_4wk_up_net6.tiff', units="in", width=10, height=8, res=300)
plot(Bulk_4wk_up_net6,edge.arrow.size=0,edge.width=abs(E(Bulk_4wk_up_net6)$width/3),vertex.color=V(Bulk_4wk_up_net6)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Bulk_4wk_up_net6)$label,vertex.label.color="black",vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,edge.color=c("#ff8989","#9fe899")[1+(E(Bulk_4wk_up_net6)$weight>0)],vertex.size=degree(Bulk_4wk_up_net6)/3,layout=layout_with_dh (Bulk_4wk_up_net6),main="Bulk_4wk_up_net6")
dev.off()

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Bulk_4wk_up_net7.tiff', units="in", width=10, height=8, res=300)
plot(Bulk_4wk_up_net7,edge.arrow.size=0,edge.width=abs(E(Bulk_4wk_up_net7)$width/3),vertex.color=V(Bulk_4wk_up_net7)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Bulk_4wk_up_net7)$label,vertex.label.color="black",vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,edge.color=c("#ff8989","#9fe899")[1+(E(Bulk_4wk_up_net7)$weight>0)],vertex.size=degree(Bulk_4wk_up_net7)/3,layout=layout_with_dh (Bulk_4wk_up_net7),main="Bulk_4wk_up_net7")
dev.off()
```


#  >>>>> Rhi_1wk_up >>>>>>>

```{r}
## import SparCC and p_value matrix into R 
Rhi_1wk_up_SparCC_cor<-read.csv("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Rhi_1wk_up_OTU_cor_sparcc.csv",header = TRUE,row.names = 1)
Rhi_1wk_up_SparCC_cor<-as.matrix(Rhi_1wk_up_SparCC_cor)
dim(Rhi_1wk_up_SparCC_cor)
Rhi_1wk_up_SparCC_cor[1:5,1:5]

# Read in corresponding p_value dataframe
Rhi_1wk_up_SparCC_pvalue<-read.csv("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Rhi_1wk_up_OTU_two_sided_pvalue.csv",header = TRUE,row.names = 1)
Rhi_1wk_up_SparCC_pvalue<-as.matrix(Rhi_1wk_up_SparCC_pvalue)
dim(Rhi_1wk_up_SparCC_pvalue)
Rhi_1wk_up_SparCC_pvalue[1:5,1:5]
```


## Generate vertices information and add relative abundance information to vertices and color information


```{r}
r_Rhi_1wk_up_phyloseq<-transform_sample_counts(Rhi_1wk_up,function(x) x/sum(x))

Rhi_1wk_up_vertex<-data.frame(OtuID=taxa_names(r_Rhi_1wk_up_phyloseq),tax_table(r_Rhi_1wk_up_phyloseq),size=taxa_sums(r_Rhi_1wk_up_phyloseq)/length(sample_sums(r_Rhi_1wk_up_phyloseq))) # Here size are the average relative abundance of each OTU across all samples within same treatment
dim(Rhi_1wk_up_vertex)
Rhi_1wk_up_vertex[1:5,]
Rhi_1wk_up_vertex<-Rhi_1wk_up_vertex[,c(1,3,7,8)]
head(Rhi_1wk_up_vertex)

Bacteria_phylum_vs_color_checklist<-data.frame(phylum=c('Proteobacteria','Actinobacteria','Bacteria_unclassified','Acidobacteria','Nitrospirae','Chloroflexi','Candidatus_Saccharibacteria','Firmicutes','Bacteroidetes','Verrucomicrobia','candidate_division_WPS-1','Gemmatimonadetes','Planctomycetes','BRC1','Latescibacteria','Chlamydiae','Others','Patescibacteria','Rokubacteria','Cyanobacteria'),colors=c('#ff9082','#98d66d','#fcef5d','#c45cc0','#45b1f9','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#bc8c38','#bcba6d','#91749e','#b2798d','#343837','#b23301','#235931'))
phylum_colors<-as.character(Bacteria_phylum_vs_color_checklist$colors)
phylum_labels<-as.character(Bacteria_phylum_vs_color_checklist$phylum)



tiff("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Phylum_colour_pie.tiff",units = 'in',width=10,height = 8,res=300,compression = 'lzw')
pie(rep (1,17),col=phylum_colors,labels=paste(phylum_labels,phylum_colors,sep=":"),radius = 1.0)
dev.off()

Rhi_1wk_up_vertex$color<-factor(Rhi_1wk_up_vertex$Phylum,levels=phylum_labels,labels=phylum_colors)
dim(Rhi_1wk_up_vertex)
Rhi_1wk_up_vertex[1:10,]
# I want to add label variable to Rhi_1wk_up_vertex
Rhi_1wk_up_vertex<-join (Rhi_1wk_up_vertex,meta,by="OtuID") 

# Originally, I used merge function, but it turned out to be problematic, because after merge it will The rows are by default lexicographically sorted on the common columns. To resolve this problem, I changed merge function to joint function.
dim(Rhi_1wk_up_vertex)
Rhi_1wk_up_vertex[1:5,]
```


## Creat links and have a look at edge adn vertex information

```{r}
# read SparCC correlation matrix and p_value matrix to igraph
Rhi_1wk_up_net<-graph_from_adjacency_matrix(Rhi_1wk_up_SparCC_cor,weighted = TRUE)
Rhi_1wk_up_net  #IGRAPH 8cfb97f DNW- 500 250000 -- 
summary(V(Rhi_1wk_up_net)$name==as.character(Rhi_1wk_up_vertex$OtuID)) 
V(Rhi_1wk_up_net)$name[1:5]
E(Rhi_1wk_up_net)$weight[1:5] # here the weight are just the SparCC coefficience

# Add vertex information to the network
V(Rhi_1wk_up_net)$vertex_color<-as.character(Rhi_1wk_up_vertex$color)
V(Rhi_1wk_up_net)$size<-Rhi_1wk_up_vertex$size
V(Rhi_1wk_up_net)$Phylum<-as.character(Rhi_1wk_up_vertex$Phylum)
V(Rhi_1wk_up_net)$label<-as.character(Rhi_1wk_up_vertex$label)

# Add edge information (p_value) to the network
E(Rhi_1wk_up_net)$p_value<-Rhi_1wk_up_SparCC_pvalue
E(Rhi_1wk_up_net)$p_value[1:10]

## Add edge width to the network
E(Rhi_1wk_up_net)$width<-E(Rhi_1wk_up_net)$weight*10
E(Rhi_1wk_up_net)$width[1:20]
```


------------------------------------------------------------------------------------------------------------------------------------------
In order to simplify the network, including remove of loop and mutual edges. If using simplify function on Rhi_1wk_up_net, this will loose p_value and edge color information (which indicate if the interaction are negarive or positive)

To solve this problem, I will write out the edge and vertices in dataframe format and read into these two dataframe and modify the edge information globally before use ``graph_from_data_frame`` to create network object.
-----------------------------------------------------------------------------------------------------------------------------------------


## write out edges and vertices in dataframe format, and edit on the network

```{r}

# -- generate edge dataframe from whole network---

Rhi_1wk_up_edges<-igraph::as_data_frame(Rhi_1wk_up_net,what = "edges")
#Rhi_1wk_up_edges
dim(Rhi_1wk_up_edges) # 250000 rows and 5 columns
head(Rhi_1wk_up_edges)
# Here the weight indicate the strength of interaction (SparCC value), p_value is the significance of interactions nad width is formulated based on weight * 40.

# -- generate vertices dataframe from whole network---

Rhi_1wk_up_vertices<-igraph::as_data_frame(Rhi_1wk_up_net,what = "vertices")
dim(Rhi_1wk_up_vertices)
head(Rhi_1wk_up_vertices)

# remove edges that from and to the same nodes(selfloops)

Rhi_1wk_up_edges<-Rhi_1wk_up_edges[which(Rhi_1wk_up_edges$from!=Rhi_1wk_up_edges$to),]
dim(Rhi_1wk_up_edges) # 1208900 rows and 5 columns
head(Rhi_1wk_up_edges)

# Creat a new network using the above vertices and edge data frame
Rhi_1wk_up_net1 <- graph_from_data_frame(d=Rhi_1wk_up_edges, vertices=Rhi_1wk_up_vertices, directed=T)#here if I used directed=FALSE, I will lost my p_value information as well as weight information. 
Rhi_1wk_up_net1 # IGRAPH 526a39b DNW- 800 639200 --
E(Rhi_1wk_up_net1)$weight[1:10]

# convert directed network to undirected network

Rhi_1wk_up_net2<-as.undirected(Rhi_1wk_up_net1 ,mode = 'mutual',edge.attr.comb = "mean") # as a result, weight, p_value and width of edges are all averaged for mutual edges.
Rhi_1wk_up_net2 # IGRAPH 1bce66c UNW- 500 124750 - -
E(Rhi_1wk_up_net2)$weight[1:10]
E(Rhi_1wk_up_net2)$p_value[1:10]
E(Rhi_1wk_up_net2)$width[1:10]


## simplify the network from Rhi_1wk_up_net2

#1)  remove non_significant edges with alpha=0.05

summary(E(Rhi_1wk_up_net2)$p_value<0.05) #TRUE 74374 and FALSE 530076   

Rhi_1wk_up_net3<-delete.edges(Rhi_1wk_up_net2,which(E(Rhi_1wk_up_net2)$p_value>=0.05))
Rhi_1wk_up_net3 # IGRAPH f25c247 UNW- 800 24383 -- 
summary(degree(Rhi_1wk_up_net3)==0) # all of the nodes has connection to others
summary(E(Rhi_1wk_up_net3)$weight>0) # 12668 positive interactions, 11715 nagative interactions

#2) remove non-significant edges with alpha=0.001

summary(E(Rhi_1wk_up_net2)$p_value<0.001) # TRUE 14033, FALSE 590417
Rhi_1wk_up_net4<-delete.edges(Rhi_1wk_up_net2,which(E(Rhi_1wk_up_net2)$p_value>=0.001))
Rhi_1wk_up_net4 # IGRAPH b2b3bc7 UNW- 800 2481 --
Rhi_1wk_up_net4<-delete.vertices(Rhi_1wk_up_net4,degree(Rhi_1wk_up_net4)==0) # IGRAPH 281b75b UNW- 789 2481 --

#3) further simplify the network based on degree
barplot(sort(degree(Rhi_1wk_up_net4)),main = "Rhi_1wk_up_net4_barplot")
summary(degree(Rhi_1wk_up_net4)>15) # 36 nodes with degree larger than 100
Rhi_1wk_up_net5<-delete.vertices(Rhi_1wk_up_net4,which(degree(Rhi_1wk_up_net4)<=5))
Rhi_1wk_up_net5# IGRAPH ce4509d UNW- 36 536 --
summary(E(Rhi_1wk_up_net5)$weight>0) # 261 positive and 275 negative

# top 50

Rhi_1wk_up_net6<-delete.vertices(Rhi_1wk_up_net4,labels(sort(degree(Rhi_1wk_up_net4),FALSE)[c(1:(length(V(Rhi_1wk_up_net4)$name)-50))])) # IGRAPH 56796e5 UNW- 50 937 --
Rhi_1wk_up_net6<-delete.vertices(Rhi_1wk_up_net6,degree(Rhi_1wk_up_net6)==0) #IGRAPH cead187 UNW- 43 71 -- 
summary(E(Rhi_1wk_up_net6)$weight>0) # logical : 46(negative)     25(positive)

# top100
Rhi_1wk_up_net7<-delete.vertices(Rhi_1wk_up_net4,labels(sort(degree(Rhi_1wk_up_net4),FALSE)[c(1:(length(V(Rhi_1wk_up_net4)$name)-100))]))
Rhi_1wk_up_net7<-delete.vertices(Rhi_1wk_up_net7,degree(Rhi_1wk_up_net7)==0) #IGRAPH cc57451 UNW- 96 200 -
summary(E(Rhi_1wk_up_net7)$weight>0) # 1273 (negative)    1228 (positive) 
```


## Plot network

```{r}

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Rhi_1wk_up_net4.tiff', units="in", width=10, height=8, res=300)

plot(Rhi_1wk_up_net4,edge.arrow.size=0,edge.width=abs(E(Rhi_1wk_up_net4)$width/4),vertex.color=V(Rhi_1wk_up_net4)$vertex_color,vertex.frame.color="#555555",vertex.label=NA,edge.color=c("#ff8989","#9fe899")[1+(E(Rhi_1wk_up_net4)$weight>0)],vertex.size=degree(Rhi_1wk_up_net4)/10,layout=layout_with_fr (Rhi_1wk_up_net4),main="Rhi_1wk_up_net4") # green edges means positive correlation and red edges mean negative correlation
dev.off()

#plot(Rhi_1wk_up_net5,edge.arrow.size=0,edge.width=abs(E(Rhi_1wk_up_net5)$width/3),vertex.color=V(Rhi_1wk_up_net5)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Rhi_1wk_up_net5)$label,edge.color=c("#ff8989","#9fe899")[1+(E(Rhi_1wk_up_net5)$weight>0)],vertex.label.cex=.8,vertex.size=degree(Rhi_1wk_up_net5)/2,layout=layout_with_dh (Rhi_1wk_up_net5),main="Rhi_1wk_up_net5")

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Rhi_1wk_up_net6.tiff', units="in", width=10, height=8, res=300)
plot(Rhi_1wk_up_net6,edge.arrow.size=0,edge.width=abs(E(Rhi_1wk_up_net6)$width/3),vertex.color=V(Rhi_1wk_up_net6)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Rhi_1wk_up_net6)$label,vertex.label.color="black",vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,edge.color=c("#ff8989","#9fe899")[1+(E(Rhi_1wk_up_net6)$weight>0)],vertex.size=degree(Rhi_1wk_up_net6)/3,layout=layout_with_dh (Rhi_1wk_up_net6),main="Rhi_1wk_up_net6")
dev.off()

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Rhi_1wk_up_net7.tiff', units="in", width=10, height=8, res=300)
plot(Rhi_1wk_up_net7,edge.arrow.size=0,edge.width=abs(E(Rhi_1wk_up_net7)$width/3),vertex.color=V(Rhi_1wk_up_net7)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Rhi_1wk_up_net7)$label,vertex.label.color="black",vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,edge.color=c("#ff8989","#9fe899")[1+(E(Rhi_1wk_up_net7)$weight>0)],vertex.size=degree(Rhi_1wk_up_net7)/3,layout=layout_with_dh (Rhi_1wk_up_net7),main="Rhi_1wk_up_net7")
dev.off()
```

#  >>>>> Rhi_3wk_up >>>>>>>

```{r}
## import SparCC and p_value matrix into R 
Rhi_3wk_up_SparCC_cor<-read.csv("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Rhi_3wk_up_OTU_cor_sparcc.csv",header = TRUE,row.names = 1)
Rhi_3wk_up_SparCC_cor<-as.matrix(Rhi_3wk_up_SparCC_cor)
dim(Rhi_3wk_up_SparCC_cor)
Rhi_3wk_up_SparCC_cor[1:5,1:5]

# Read in corresponding p_value dataframe
Rhi_3wk_up_SparCC_pvalue<-read.csv("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Rhi_3wk_up_OTU_two_sided_pvalue.csv",header = TRUE,row.names = 1)
Rhi_3wk_up_SparCC_pvalue<-as.matrix(Rhi_3wk_up_SparCC_pvalue)
dim(Rhi_3wk_up_SparCC_pvalue)
Rhi_3wk_up_SparCC_pvalue[1:5,1:5]
```


## Generate vertices information and add relative abundance information to vertices and color information


```{r}
r_Rhi_3wk_up_phyloseq<-transform_sample_counts(Rhi_3wk_up,function(x) x/sum(x))

Rhi_3wk_up_vertex<-data.frame(OtuID=taxa_names(r_Rhi_3wk_up_phyloseq),tax_table(r_Rhi_3wk_up_phyloseq),size=taxa_sums(r_Rhi_3wk_up_phyloseq)/length(sample_sums(r_Rhi_3wk_up_phyloseq))) # Here size are the average relative abundance of each OTU across all samples within same treatment
dim(Rhi_3wk_up_vertex)
Rhi_3wk_up_vertex[1:5,]
Rhi_3wk_up_vertex<-Rhi_3wk_up_vertex[,c(1,3,7,8)]
head(Rhi_3wk_up_vertex)

Bacteria_phylum_vs_color_checklist<-data.frame(phylum=c('Proteobacteria','Actinobacteria','Bacteria_unclassified','Acidobacteria','Nitrospirae','Chloroflexi','Candidatus_Saccharibacteria','Firmicutes','Bacteroidetes','Verrucomicrobia','candidate_division_WPS-1','Gemmatimonadetes','Planctomycetes','BRC1','Latescibacteria','Chlamydiae','Others','Patescibacteria','Rokubacteria','Cyanobacteria'),colors=c('#ff9082','#98d66d','#fcef5d','#c45cc0','#45b1f9','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#bc8c38','#bcba6d','#91749e','#b2798d','#343837','#b23301','#235931'))
phylum_colors<-as.character(Bacteria_phylum_vs_color_checklist$colors)
phylum_labels<-as.character(Bacteria_phylum_vs_color_checklist$phylum)


tiff("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Phylum_colour_pie.tiff",units = 'in',width=10,height = 8,res=300,compression = 'lzw')
pie(rep (1,17),col=phylum_colors,labels=paste(phylum_labels,phylum_colors,sep=":"),radius = 1.0)
dev.off()

Rhi_3wk_up_vertex$color<-factor(Rhi_3wk_up_vertex$Phylum,levels=phylum_labels,labels=phylum_colors)
dim(Rhi_3wk_up_vertex)
Rhi_3wk_up_vertex[1:10,]
# I want to add label variable to Rhi_3wk_up_vertex
Rhi_3wk_up_vertex<-join (Rhi_3wk_up_vertex,meta,by="OtuID") 

# Originally, I used merge function, but it turned out to be problematic, because after merge it will The rows are by default lexicographically sorted on the common columns. To resolve this problem, I changed merge function to joint function.
dim(Rhi_3wk_up_vertex)
Rhi_3wk_up_vertex[1:5,]
```


## Creat links and have a look at edge adn vertex information

```{r}
# read SparCC correlation matrix and p_value matrix to igraph
Rhi_3wk_up_net<-graph_from_adjacency_matrix(Rhi_3wk_up_SparCC_cor,weighted = TRUE)
Rhi_3wk_up_net  #IGRAPH 8cfb97f DNW- 500 250000 -- 
summary(V(Rhi_3wk_up_net)$name==as.character(Rhi_3wk_up_vertex$OtuID)) 
V(Rhi_3wk_up_net)$name[1:5]
E(Rhi_3wk_up_net)$weight[1:5] # here the weight are just the SparCC coefficience

# Add vertex information to the network
V(Rhi_3wk_up_net)$vertex_color<-as.character(Rhi_3wk_up_vertex$color)
V(Rhi_3wk_up_net)$size<-Rhi_3wk_up_vertex$size
V(Rhi_3wk_up_net)$Phylum<-as.character(Rhi_3wk_up_vertex$Phylum)
V(Rhi_3wk_up_net)$label<-as.character(Rhi_3wk_up_vertex$label)

# Add edge information (p_value) to the network
E(Rhi_3wk_up_net)$p_value<-Rhi_3wk_up_SparCC_pvalue
E(Rhi_3wk_up_net)$p_value[1:10]

## Add edge width to the network
E(Rhi_3wk_up_net)$width<-E(Rhi_3wk_up_net)$weight*10
E(Rhi_3wk_up_net)$width[1:20]
```


------------------------------------------------------------------------------------------------------------------------------------------
In order to simplify the network, including remove of loop and mutual edges. If using simplify function on Rhi_3wk_up_net, this will loose p_value and edge color information (which indicate if the interaction are negarive or positive)

To solve this problem, I will write out the edge and vertices in dataframe format and read into these two dataframe and modify the edge information globally before use ``graph_from_data_frame`` to create network object.
-----------------------------------------------------------------------------------------------------------------------------------------


## write out edges and vertices in dataframe format, and edit on the network

```{r}

# -- generate edge dataframe from whole network---

Rhi_3wk_up_edges<-igraph::as_data_frame(Rhi_3wk_up_net,what = "edges")
#Rhi_3wk_up_edges
dim(Rhi_3wk_up_edges) # 250000 rows and 5 columns
head(Rhi_3wk_up_edges)
# Here the weight indicate the strength of interaction (SparCC value), p_value is the significance of interactions nad width is formulated based on weight * 40.

# -- generate vertices dataframe from whole network---

Rhi_3wk_up_vertices<-igraph::as_data_frame(Rhi_3wk_up_net,what = "vertices")
dim(Rhi_3wk_up_vertices)
head(Rhi_3wk_up_vertices)

# remove edges that from and to the same nodes(selfloops)

Rhi_3wk_up_edges<-Rhi_3wk_up_edges[which(Rhi_3wk_up_edges$from!=Rhi_3wk_up_edges$to),]
dim(Rhi_3wk_up_edges) # 1208900 rows and 5 columns
head(Rhi_3wk_up_edges)

# Creat a new network using the above vertices and edge data frame
Rhi_3wk_up_net1 <- graph_from_data_frame(d=Rhi_3wk_up_edges, vertices=Rhi_3wk_up_vertices, directed=T)#here if I used directed=FALSE, I will lost my p_value information as well as weight information. 
Rhi_3wk_up_net1 # IGRAPH 526a39b DNW- 800 639200 --
E(Rhi_3wk_up_net1)$weight[1:10]

# convert directed network to undirected network

Rhi_3wk_up_net2<-as.undirected(Rhi_3wk_up_net1 ,mode = 'mutual',edge.attr.comb = "mean") # as a result, weight, p_value and width of edges are all averaged for mutual edges.
Rhi_3wk_up_net2 # IGRAPH 1bce66c UNW- 500 124750 - -
E(Rhi_3wk_up_net2)$weight[1:10]
E(Rhi_3wk_up_net2)$p_value[1:10]
E(Rhi_3wk_up_net2)$width[1:10]


## simplify the network from Rhi_3wk_up_net2

#1)  remove non_significant edges with alpha=0.05

summary(E(Rhi_3wk_up_net2)$p_value<0.05) #TRUE 74374 and FALSE 530076   

Rhi_3wk_up_net3<-delete.edges(Rhi_3wk_up_net2,which(E(Rhi_3wk_up_net2)$p_value>=0.05))
Rhi_3wk_up_net3 # IGRAPH f25c247 UNW- 800 24383 -- 
summary(degree(Rhi_3wk_up_net3)==0) # all of the nodes has connection to others
summary(E(Rhi_3wk_up_net3)$weight>0) # 12668 positive interactions, 11715 nagative interactions

#2) remove non-significant edges with alpha=0.001

summary(E(Rhi_3wk_up_net2)$p_value<0.001) # TRUE 14033, FALSE 590417
Rhi_3wk_up_net4<-delete.edges(Rhi_3wk_up_net2,which(E(Rhi_3wk_up_net2)$p_value>=0.001))
Rhi_3wk_up_net4 # IGRAPH b2b3bc7 UNW- 800 2481 --
Rhi_3wk_up_net4<-delete.vertices(Rhi_3wk_up_net4,degree(Rhi_3wk_up_net4)==0) # IGRAPH 281b75b UNW- 789 2481 --

#3) further simplify the network based on degree
barplot(sort(degree(Rhi_3wk_up_net4)),main = "Rhi_3wk_up_net4_barplot")
summary(degree(Rhi_3wk_up_net4)>15) # 36 nodes with degree larger than 100
Rhi_3wk_up_net5<-delete.vertices(Rhi_3wk_up_net4,which(degree(Rhi_3wk_up_net4)<=5))
Rhi_3wk_up_net5# IGRAPH ce4509d UNW- 36 536 --
summary(E(Rhi_3wk_up_net5)$weight>0) # 261 positive and 275 negative

# top 50

Rhi_3wk_up_net6<-delete.vertices(Rhi_3wk_up_net4,labels(sort(degree(Rhi_3wk_up_net4),FALSE)[c(1:(length(V(Rhi_3wk_up_net4)$name)-50))])) # IGRAPH 56796e5 UNW- 50 937 --
Rhi_3wk_up_net6<-delete.vertices(Rhi_3wk_up_net6,degree(Rhi_3wk_up_net6)==0) #IGRAPH cead187 UNW- 43 71 -- 
summary(E(Rhi_3wk_up_net6)$weight>0) # logical : 46(negative)     25(positive)

# top100
Rhi_3wk_up_net7<-delete.vertices(Rhi_3wk_up_net4,labels(sort(degree(Rhi_3wk_up_net4),FALSE)[c(1:(length(V(Rhi_3wk_up_net4)$name)-100))]))
Rhi_3wk_up_net7<-delete.vertices(Rhi_3wk_up_net7,degree(Rhi_3wk_up_net7)==0) #IGRAPH cc57451 UNW- 96 200 -
summary(E(Rhi_3wk_up_net7)$weight>0) # 1273 (negative)    1228 (positive) 
```


## Plot network

```{r}

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Rhi_3wk_up_net4.tiff', units="in", width=10, height=8, res=300)

plot(Rhi_3wk_up_net4,edge.arrow.size=0,edge.width=abs(E(Rhi_3wk_up_net4)$width/4),vertex.color=V(Rhi_3wk_up_net4)$vertex_color,vertex.frame.color="#555555",vertex.label=NA,edge.color=c("#ff8989","#9fe899")[1+(E(Rhi_3wk_up_net4)$weight>0)],vertex.size=degree(Rhi_3wk_up_net4)/10,layout=layout_with_fr (Rhi_3wk_up_net4),main="Rhi_3wk_up_net4") # green edges means positive correlation and red edges mean negative correlation
dev.off()

#plot(Rhi_3wk_up_net5,edge.arrow.size=0,edge.width=abs(E(Rhi_3wk_up_net5)$width/3),vertex.color=V(Rhi_3wk_up_net5)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Rhi_3wk_up_net5)$label,edge.color=c("#ff8989","#9fe899")[1+(E(Rhi_3wk_up_net5)$weight>0)],vertex.label.cex=.8,vertex.size=degree(Rhi_3wk_up_net5)/2,layout=layout_with_dh (Rhi_3wk_up_net5),main="Rhi_3wk_up_net5")

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Rhi_3wk_up_net6.tiff', units="in", width=10, height=8, res=300)
plot(Rhi_3wk_up_net6,edge.arrow.size=0,edge.width=abs(E(Rhi_3wk_up_net6)$width/3),vertex.color=V(Rhi_3wk_up_net6)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Rhi_3wk_up_net6)$label,vertex.label.color="black",vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,edge.color=c("#ff8989","#9fe899")[1+(E(Rhi_3wk_up_net6)$weight>0)],vertex.size=degree(Rhi_3wk_up_net6)/3,layout=layout_with_dh (Rhi_3wk_up_net6),main="Rhi_3wk_up_net6")
dev.off()

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Rhi_3wk_up_net7.tiff', units="in", width=10, height=8, res=300)
plot(Rhi_3wk_up_net7,edge.arrow.size=0,edge.width=abs(E(Rhi_3wk_up_net7)$width/3),vertex.color=V(Rhi_3wk_up_net7)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Rhi_3wk_up_net7)$label,vertex.label.color="black",vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,edge.color=c("#ff8989","#9fe899")[1+(E(Rhi_3wk_up_net7)$weight>0)],vertex.size=degree(Rhi_3wk_up_net7)/3,layout=layout_with_dh (Rhi_3wk_up_net7),main="Rhi_3wk_up_net7")
dev.off()
```

#  >>>>> Rhi_4wk_up >>>>>>>

```{r}
## import SparCC and p_value matrix into R 
Rhi_4wk_up_SparCC_cor<-read.csv("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Rhi_4wk_up_OTU_cor_sparcc.csv",header = TRUE,row.names = 1)
Rhi_4wk_up_SparCC_cor<-as.matrix(Rhi_4wk_up_SparCC_cor)
dim(Rhi_4wk_up_SparCC_cor)
Rhi_4wk_up_SparCC_cor[1:5,1:5]

# Read in corresponding p_value dataframe
Rhi_4wk_up_SparCC_pvalue<-read.csv("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Rhi_4wk_up_OTU_two_sided_pvalue.csv",header = TRUE,row.names = 1)
Rhi_4wk_up_SparCC_pvalue<-as.matrix(Rhi_4wk_up_SparCC_pvalue)
dim(Rhi_4wk_up_SparCC_pvalue)
Rhi_4wk_up_SparCC_pvalue[1:5,1:5]
```


## Generate vertices information and add relative abundance information to vertices and color information


```{r}
r_Rhi_4wk_up_phyloseq<-transform_sample_counts(Rhi_4wk_up,function(x) x/sum(x))

Rhi_4wk_up_vertex<-data.frame(OtuID=taxa_names(r_Rhi_4wk_up_phyloseq),tax_table(r_Rhi_4wk_up_phyloseq),size=taxa_sums(r_Rhi_4wk_up_phyloseq)/length(sample_sums(r_Rhi_4wk_up_phyloseq))) # Here size are the average relative abundance of each OTU across all samples within same treatment
dim(Rhi_4wk_up_vertex)
Rhi_4wk_up_vertex[1:5,]
Rhi_4wk_up_vertex<-Rhi_4wk_up_vertex[,c(1,3,7,8)]
head(Rhi_4wk_up_vertex)

Bacteria_phylum_vs_color_checklist<-data.frame(phylum=c('Proteobacteria','Actinobacteria','Bacteria_unclassified','Acidobacteria','Nitrospirae','Chloroflexi','Candidatus_Saccharibacteria','Firmicutes','Bacteroidetes','Verrucomicrobia','candidate_division_WPS-1','Gemmatimonadetes','Planctomycetes','BRC1','Latescibacteria','Chlamydiae','Others','Patescibacteria','Rokubacteria','Cyanobacteria'),colors=c('#ff9082','#98d66d','#fcef5d','#c45cc0','#45b1f9','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#bc8c38','#bcba6d','#91749e','#b2798d','#343837','#b23301','#235931'))
phylum_colors<-as.character(Bacteria_phylum_vs_color_checklist$colors)
phylum_labels<-as.character(Bacteria_phylum_vs_color_checklist$phylum)


tiff("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Phylum_colour_pie.tiff",units = 'in',width=10,height = 8,res=300,compression = 'lzw')
pie(rep (1,17),col=phylum_colors,labels=paste(phylum_labels,phylum_colors,sep=":"),radius = 1.0)
dev.off()

Rhi_4wk_up_vertex$color<-factor(Rhi_4wk_up_vertex$Phylum,levels=phylum_labels,labels=phylum_colors)
dim(Rhi_4wk_up_vertex)
Rhi_4wk_up_vertex[1:10,]
# I want to add label variable to Rhi_4wk_up_vertex
Rhi_4wk_up_vertex<-join (Rhi_4wk_up_vertex,meta,by="OtuID") 

# Originally, I used merge function, but it turned out to be problematic, because after merge it will The rows are by default lexicographically sorted on the common columns. To resolve this problem, I changed merge function to joint function.
dim(Rhi_4wk_up_vertex)
Rhi_4wk_up_vertex[1:5,]
```


## Creat links and have a look at edge adn vertex information

```{r}
# read SparCC correlation matrix and p_value matrix to igraph
Rhi_4wk_up_net<-graph_from_adjacency_matrix(Rhi_4wk_up_SparCC_cor,weighted = TRUE)
Rhi_4wk_up_net  #IGRAPH 8cfb97f DNW- 500 250000 -- 
summary(V(Rhi_4wk_up_net)$name==as.character(Rhi_4wk_up_vertex$OtuID)) 
V(Rhi_4wk_up_net)$name[1:5]
E(Rhi_4wk_up_net)$weight[1:5] # here the weight are just the SparCC coefficience

# Add vertex information to the network
V(Rhi_4wk_up_net)$vertex_color<-as.character(Rhi_4wk_up_vertex$color)
V(Rhi_4wk_up_net)$size<-Rhi_4wk_up_vertex$size
V(Rhi_4wk_up_net)$Phylum<-as.character(Rhi_4wk_up_vertex$Phylum)
V(Rhi_4wk_up_net)$label<-as.character(Rhi_4wk_up_vertex$label)

# Add edge information (p_value) to the network
E(Rhi_4wk_up_net)$p_value<-Rhi_4wk_up_SparCC_pvalue
E(Rhi_4wk_up_net)$p_value[1:10]

## Add edge width to the network
E(Rhi_4wk_up_net)$width<-E(Rhi_4wk_up_net)$weight*10
E(Rhi_4wk_up_net)$width[1:20]
```


------------------------------------------------------------------------------------------------------------------------------------------
In order to simplify the network, including remove of loop and mutual edges. If using simplify function on Rhi_4wk_up_net, this will loose p_value and edge color information (which indicate if the interaction are negarive or positive)

To solve this problem, I will write out the edge and vertices in dataframe format and read into these two dataframe and modify the edge information globally before use ``graph_from_data_frame`` to create network object.
-----------------------------------------------------------------------------------------------------------------------------------------


## write out edges and vertices in dataframe format, and edit on the network

```{r}

# -- generate edge dataframe from whole network---

Rhi_4wk_up_edges<-igraph::as_data_frame(Rhi_4wk_up_net,what = "edges")
#Rhi_4wk_up_edges
dim(Rhi_4wk_up_edges) # 250000 rows and 5 columns
head(Rhi_4wk_up_edges)
# Here the weight indicate the strength of interaction (SparCC value), p_value is the significance of interactions nad width is formulated based on weight * 40.

# -- generate vertices dataframe from whole network---

Rhi_4wk_up_vertices<-igraph::as_data_frame(Rhi_4wk_up_net,what = "vertices")
dim(Rhi_4wk_up_vertices)
head(Rhi_4wk_up_vertices)

# remove edges that from and to the same nodes(selfloops)

Rhi_4wk_up_edges<-Rhi_4wk_up_edges[which(Rhi_4wk_up_edges$from!=Rhi_4wk_up_edges$to),]
dim(Rhi_4wk_up_edges) # 1208900 rows and 5 columns
head(Rhi_4wk_up_edges)

# Creat a new network using the above vertices and edge data frame
Rhi_4wk_up_net1 <- graph_from_data_frame(d=Rhi_4wk_up_edges, vertices=Rhi_4wk_up_vertices, directed=T)#here if I used directed=FALSE, I will lost my p_value information as well as weight information. 
Rhi_4wk_up_net1 # IGRAPH 526a39b DNW- 800 639200 --
E(Rhi_4wk_up_net1)$weight[1:10]

# convert directed network to undirected network

Rhi_4wk_up_net2<-as.undirected(Rhi_4wk_up_net1 ,mode = 'mutual',edge.attr.comb = "mean") # as a result, weight, p_value and width of edges are all averaged for mutual edges.
Rhi_4wk_up_net2 # IGRAPH 1bce66c UNW- 500 124750 - -
E(Rhi_4wk_up_net2)$weight[1:10]
E(Rhi_4wk_up_net2)$p_value[1:10]
E(Rhi_4wk_up_net2)$width[1:10]


## simplify the network from Rhi_4wk_up_net2

#1)  remove non_significant edges with alpha=0.05

summary(E(Rhi_4wk_up_net2)$p_value<0.05) #TRUE 74374 and FALSE 530076   

Rhi_4wk_up_net3<-delete.edges(Rhi_4wk_up_net2,which(E(Rhi_4wk_up_net2)$p_value>=0.05))
Rhi_4wk_up_net3 # IGRAPH f25c247 UNW- 800 24383 -- 
summary(degree(Rhi_4wk_up_net3)==0) # all of the nodes has connection to others
summary(E(Rhi_4wk_up_net3)$weight>0) # 12668 positive interactions, 11715 nagative interactions

#2) remove non-significant edges with alpha=0.001

summary(E(Rhi_4wk_up_net2)$p_value<0.001) # TRUE 14033, FALSE 590417
Rhi_4wk_up_net4<-delete.edges(Rhi_4wk_up_net2,which(E(Rhi_4wk_up_net2)$p_value>=0.001))
Rhi_4wk_up_net4 # IGRAPH b2b3bc7 UNW- 800 2481 --
Rhi_4wk_up_net4<-delete.vertices(Rhi_4wk_up_net4,degree(Rhi_4wk_up_net4)==0) # IGRAPH 281b75b UNW- 789 2481 --

#3) further simplify the network based on degree
barplot(sort(degree(Rhi_4wk_up_net4)),main = "Rhi_4wk_up_net4_barplot")
summary(degree(Rhi_4wk_up_net4)>15) # 36 nodes with degree larger than 100
Rhi_4wk_up_net5<-delete.vertices(Rhi_4wk_up_net4,which(degree(Rhi_4wk_up_net4)<=5))
Rhi_4wk_up_net5# IGRAPH ce4509d UNW- 36 536 --
summary(E(Rhi_4wk_up_net5)$weight>0) # 261 positive and 275 negative

# top 50

Rhi_4wk_up_net6<-delete.vertices(Rhi_4wk_up_net4,labels(sort(degree(Rhi_4wk_up_net4),FALSE)[c(1:(length(V(Rhi_4wk_up_net4)$name)-50))])) # IGRAPH 56796e5 UNW- 50 937 --
Rhi_4wk_up_net6<-delete.vertices(Rhi_4wk_up_net6,degree(Rhi_4wk_up_net6)==0) #IGRAPH cead187 UNW- 43 71 -- 
summary(E(Rhi_4wk_up_net6)$weight>0) # logical : 46(negative)     25(positive)

# top100
Rhi_4wk_up_net7<-delete.vertices(Rhi_4wk_up_net4,labels(sort(degree(Rhi_4wk_up_net4),FALSE)[c(1:(length(V(Rhi_4wk_up_net4)$name)-100))]))
Rhi_4wk_up_net7<-delete.vertices(Rhi_4wk_up_net7,degree(Rhi_4wk_up_net7)==0) #IGRAPH cc57451 UNW- 96 200 -
summary(E(Rhi_4wk_up_net7)$weight>0) # 1273 (negative)    1228 (positive) 
```

## Plot network

```{r}

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Rhi_4wk_up_net4.tiff', units="in", width=10, height=8, res=300)

plot(Rhi_4wk_up_net4,edge.arrow.size=0,edge.width=abs(E(Rhi_4wk_up_net4)$width/4),vertex.color=V(Rhi_4wk_up_net4)$vertex_color,vertex.frame.color="#555555",vertex.label=NA,edge.color=c("#ff8989","#9fe899")[1+(E(Rhi_4wk_up_net4)$weight>0)],vertex.size=degree(Rhi_4wk_up_net4)/10,layout=layout_with_fr (Rhi_4wk_up_net4),main="Rhi_4wk_up_net4") # green edges means positive correlation and red edges mean negative correlation
dev.off()

#plot(Rhi_4wk_up_net5,edge.arrow.size=0,edge.width=abs(E(Rhi_4wk_up_net5)$width/3),vertex.color=V(Rhi_4wk_up_net5)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Rhi_4wk_up_net5)$label,edge.color=c("#ff8989","#9fe899")[1+(E(Rhi_4wk_up_net5)$weight>0)],vertex.label.cex=.8,vertex.size=degree(Rhi_4wk_up_net5)/2,layout=layout_with_dh (Rhi_4wk_up_net5),main="Rhi_4wk_up_net5")

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Rhi_4wk_up_net6.tiff', units="in", width=10, height=8, res=300)
plot(Rhi_4wk_up_net6,edge.arrow.size=0,edge.width=abs(E(Rhi_4wk_up_net6)$width/3),vertex.color=V(Rhi_4wk_up_net6)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Rhi_4wk_up_net6)$label,vertex.label.color="black",vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,edge.color=c("#ff8989","#9fe899")[1+(E(Rhi_4wk_up_net6)$weight>0)],vertex.size=degree(Rhi_4wk_up_net6)/3,layout=layout_with_dh (Rhi_4wk_up_net6),main="Rhi_4wk_up_net6")
dev.off()

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Rhi_4wk_up_net7.tiff', units="in", width=10, height=8, res=300)
plot(Rhi_4wk_up_net7,edge.arrow.size=0,edge.width=abs(E(Rhi_4wk_up_net7)$width/3),vertex.color=V(Rhi_4wk_up_net7)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Rhi_4wk_up_net7)$label,vertex.label.color="black",vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,edge.color=c("#ff8989","#9fe899")[1+(E(Rhi_4wk_up_net7)$weight>0)],vertex.size=degree(Rhi_4wk_up_net7)/3,layout=layout_with_dh (Rhi_4wk_up_net7),main="Rhi_4wk_up_net7")
dev.off()
```

#  >>>>> Endo_1wk_up >>>>>>>

```{r}
## import SparCC and p_value matrix into R 
Endo_1wk_up_SparCC_cor<-read.csv("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Endo_1wk_up_OTU_cor_sparcc.csv",header = TRUE,row.names = 1)
Endo_1wk_up_SparCC_cor<-as.matrix(Endo_1wk_up_SparCC_cor)
dim(Endo_1wk_up_SparCC_cor)
Endo_1wk_up_SparCC_cor[1:5,1:5]

# Read in corresponding p_value dataframe
Endo_1wk_up_SparCC_pvalue<-read.csv("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Endo_1wk_up_OTU_two_sided_pvalue.csv",header = TRUE,row.names = 1)
Endo_1wk_up_SparCC_pvalue<-as.matrix(Endo_1wk_up_SparCC_pvalue)
dim(Endo_1wk_up_SparCC_pvalue)
Endo_1wk_up_SparCC_pvalue[1:5,1:5]
```


## Generate vertices information and add relative abundance information to vertices and color information


```{r}
r_Endo_1wk_up_phyloseq<-transform_sample_counts(Endo_1wk_up,function(x) x/sum(x))

Endo_1wk_up_vertex<-data.frame(OtuID=taxa_names(r_Endo_1wk_up_phyloseq),tax_table(r_Endo_1wk_up_phyloseq),size=taxa_sums(r_Endo_1wk_up_phyloseq)/length(sample_sums(r_Endo_1wk_up_phyloseq))) # Here size are the average relative abundance of each OTU across all samples within same treatment
dim(Endo_1wk_up_vertex)
Endo_1wk_up_vertex[1:5,]
Endo_1wk_up_vertex<-Endo_1wk_up_vertex[,c(1,3,7,8)]
head(Endo_1wk_up_vertex)

Bacteria_phylum_vs_color_checklist<-data.frame(phylum=c('Proteobacteria','Actinobacteria','Bacteria_unclassified','Acidobacteria','Nitrospirae','Chloroflexi','Candidatus_Saccharibacteria','Firmicutes','Bacteroidetes','Verrucomicrobia','candidate_division_WPS-1','Gemmatimonadetes','Planctomycetes','BRC1','Latescibacteria','Chlamydiae','Others','Patescibacteria','Rokubacteria','Cyanobacteria'),colors=c('#ff9082','#98d66d','#fcef5d','#c45cc0','#45b1f9','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#bc8c38','#bcba6d','#91749e','#b2798d','#343837','#b23301','#235931'))
phylum_colors<-as.character(Bacteria_phylum_vs_color_checklist$colors)
phylum_labels<-as.character(Bacteria_phylum_vs_color_checklist$phylum)


tiff("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Phylum_colour_pie.tiff",units = 'in',width=10,height = 8,res=300,compression = 'lzw')
pie(rep (1,17),col=phylum_colors,labels=paste(phylum_labels,phylum_colors,sep=":"),radius = 1.0)
dev.off()

Endo_1wk_up_vertex$color<-factor(Endo_1wk_up_vertex$Phylum,levels=phylum_labels,labels=phylum_colors)
dim(Endo_1wk_up_vertex)
Endo_1wk_up_vertex[1:10,]
# I want to add label variable to Endo_1wk_up_vertex
Endo_1wk_up_vertex<-join (Endo_1wk_up_vertex,meta,by="OtuID") 

# Originally, I used merge function, but it turned out to be problematic, because after merge it will The rows are by default lexicographically sorted on the common columns. To resolve this problem, I changed merge function to joint function.
dim(Endo_1wk_up_vertex)
Endo_1wk_up_vertex[1:5,]
```


## Creat links and have a look at edge adn vertex information

```{r}
# read SparCC correlation matrix and p_value matrix to igraph
Endo_1wk_up_net<-graph_from_adjacency_matrix(Endo_1wk_up_SparCC_cor,weighted = TRUE)
Endo_1wk_up_net  #IGRAPH 8cfb97f DNW- 500 250000 -- 
summary(V(Endo_1wk_up_net)$name==as.character(Endo_1wk_up_vertex$OtuID)) 
V(Endo_1wk_up_net)$name[1:5]
E(Endo_1wk_up_net)$weight[1:5] # here the weight are just the SparCC coefficience

# Add vertex information to the network
V(Endo_1wk_up_net)$vertex_color<-as.character(Endo_1wk_up_vertex$color)
V(Endo_1wk_up_net)$size<-Endo_1wk_up_vertex$size
V(Endo_1wk_up_net)$Phylum<-as.character(Endo_1wk_up_vertex$Phylum)
V(Endo_1wk_up_net)$label<-as.character(Endo_1wk_up_vertex$label)

# Add edge information (p_value) to the network
E(Endo_1wk_up_net)$p_value<-Endo_1wk_up_SparCC_pvalue
E(Endo_1wk_up_net)$p_value[1:10]

## Add edge width to the network
E(Endo_1wk_up_net)$width<-E(Endo_1wk_up_net)$weight*10
E(Endo_1wk_up_net)$width[1:20]
```


------------------------------------------------------------------------------------------------------------------------------------------
In order to simplify the network, including remove of loop and mutual edges. If using simplify function on Endo_1wk_up_net, this will loose p_value and edge color information (which indicate if the interaction are negarive or positive)

To solve this problem, I will write out the edge and vertices in dataframe format and read into these two dataframe and modify the edge information globally before use ``graph_from_data_frame`` to create network object.
-----------------------------------------------------------------------------------------------------------------------------------------


## write out edges and vertices in dataframe format, and edit on the network

```{r}

# -- generate edge dataframe from whole network---

Endo_1wk_up_edges<-igraph::as_data_frame(Endo_1wk_up_net,what = "edges")
#Endo_1wk_up_edges
dim(Endo_1wk_up_edges) # 250000 rows and 5 columns
head(Endo_1wk_up_edges)
# Here the weight indicate the strength of interaction (SparCC value), p_value is the significance of interactions nad width is formulated based on weight * 40.

# -- generate vertices dataframe from whole network---

Endo_1wk_up_vertices<-igraph::as_data_frame(Endo_1wk_up_net,what = "vertices")
dim(Endo_1wk_up_vertices)
head(Endo_1wk_up_vertices)

# remove edges that from and to the same nodes(selfloops)

Endo_1wk_up_edges<-Endo_1wk_up_edges[which(Endo_1wk_up_edges$from!=Endo_1wk_up_edges$to),]
dim(Endo_1wk_up_edges) # 1208900 rows and 5 columns
head(Endo_1wk_up_edges)

# Creat a new network using the above vertices and edge data frame
Endo_1wk_up_net1 <- graph_from_data_frame(d=Endo_1wk_up_edges, vertices=Endo_1wk_up_vertices, directed=T)#here if I used directed=FALSE, I will lost my p_value information as well as weight information. 
Endo_1wk_up_net1 # IGRAPH 526a39b DNW- 800 639200 --
E(Endo_1wk_up_net1)$weight[1:10]

# convert directed network to undirected network

Endo_1wk_up_net2<-as.undirected(Endo_1wk_up_net1 ,mode = 'mutual',edge.attr.comb = "mean") # as a result, weight, p_value and width of edges are all averaged for mutual edges.
Endo_1wk_up_net2 # IGRAPH 1bce66c UNW- 500 124750 - -
E(Endo_1wk_up_net2)$weight[1:10]
E(Endo_1wk_up_net2)$p_value[1:10]
E(Endo_1wk_up_net2)$width[1:10]


## simplify the network from Endo_1wk_up_net2

#1)  remove non_significant edges with alpha=0.05

summary(E(Endo_1wk_up_net2)$p_value<0.05) #TRUE 74374 and FALSE 530076   

Endo_1wk_up_net3<-delete.edges(Endo_1wk_up_net2,which(E(Endo_1wk_up_net2)$p_value>=0.05))
Endo_1wk_up_net3 # IGRAPH f25c247 UNW- 800 24383 -- 
summary(degree(Endo_1wk_up_net3)==0) # all of the nodes has connection to others
summary(E(Endo_1wk_up_net3)$weight>0) # 12668 positive interactions, 11715 nagative interactions

#2) remove non-significant edges with alpha=0.001

summary(E(Endo_1wk_up_net2)$p_value<0.001) # TRUE 14033, FALSE 590417
Endo_1wk_up_net4<-delete.edges(Endo_1wk_up_net2,which(E(Endo_1wk_up_net2)$p_value>=0.001))
Endo_1wk_up_net4 # IGRAPH b2b3bc7 UNW- 800 2481 --
Endo_1wk_up_net4<-delete.vertices(Endo_1wk_up_net4,degree(Endo_1wk_up_net4)==0) # IGRAPH 281b75b UNW- 789 2481 --

#3) further simplify the network based on degree
barplot(sort(degree(Endo_1wk_up_net4)),main = "Endo_1wk_up_net4_barplot")
summary(degree(Endo_1wk_up_net4)>15) # 36 nodes with degree larger than 100
Endo_1wk_up_net5<-delete.vertices(Endo_1wk_up_net4,which(degree(Endo_1wk_up_net4)<=5))
Endo_1wk_up_net5# IGRAPH ce4509d UNW- 36 536 --
summary(E(Endo_1wk_up_net5)$weight>0) # 261 positive and 275 negative

# top 50

Endo_1wk_up_net6<-delete.vertices(Endo_1wk_up_net4,labels(sort(degree(Endo_1wk_up_net4),FALSE)[c(1:(length(V(Endo_1wk_up_net4)$name)-50))])) # IGRAPH 56796e5 UNW- 50 937 --
Endo_1wk_up_net6<-delete.vertices(Endo_1wk_up_net6,degree(Endo_1wk_up_net6)==0) #IGRAPH cead187 UNW- 43 71 -- 
summary(E(Endo_1wk_up_net6)$weight>0) # logical : 46(negative)     25(positive)

# top100
Endo_1wk_up_net7<-delete.vertices(Endo_1wk_up_net4,labels(sort(degree(Endo_1wk_up_net4),FALSE)[c(1:(length(V(Endo_1wk_up_net4)$name)-100))]))
Endo_1wk_up_net7<-delete.vertices(Endo_1wk_up_net7,degree(Endo_1wk_up_net7)==0) #IGRAPH cc57451 UNW- 96 200 -
summary(E(Endo_1wk_up_net7)$weight>0) # 1273 (negative)    1228 (positive) 
```


## Plot network

```{r}

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Endo_1wk_up_net4.tiff', units="in", width=10, height=8, res=300)

plot(Endo_1wk_up_net4,edge.arrow.size=0,edge.width=abs(E(Endo_1wk_up_net4)$width/4),vertex.color=V(Endo_1wk_up_net4)$vertex_color,vertex.frame.color="#555555",vertex.label=NA,edge.color=c("#ff8989","#9fe899")[1+(E(Endo_1wk_up_net4)$weight>0)],vertex.size=degree(Endo_1wk_up_net4)/10,layout=layout_with_fr (Endo_1wk_up_net4),main="Endo_1wk_up_net4") # green edges means positive correlation and red edges mean negative correlation
dev.off()

#plot(Endo_1wk_up_net5,edge.arrow.size=0,edge.width=abs(E(Endo_1wk_up_net5)$width/3),vertex.color=V(Endo_1wk_up_net5)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Endo_1wk_up_net5)$label,edge.color=c("#ff8989","#9fe899")[1+(E(Endo_1wk_up_net5)$weight>0)],vertex.label.cex=.8,vertex.size=degree(Endo_1wk_up_net5)/2,layout=layout_with_dh (Endo_1wk_up_net5),main="Endo_1wk_up_net5")

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Endo_1wk_up_net6.tiff', units="in", width=10, height=8, res=300)
plot(Endo_1wk_up_net6,edge.arrow.size=0,edge.width=abs(E(Endo_1wk_up_net6)$width/3),vertex.color=V(Endo_1wk_up_net6)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Endo_1wk_up_net6)$label,vertex.label.color="black",vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,edge.color=c("#ff8989","#9fe899")[1+(E(Endo_1wk_up_net6)$weight>0)],vertex.size=degree(Endo_1wk_up_net6)/3,layout=layout_with_dh (Endo_1wk_up_net6),main="Endo_1wk_up_net6")
dev.off()

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Endo_1wk_up_net7.tiff', units="in", width=10, height=8, res=300)
plot(Endo_1wk_up_net7,edge.arrow.size=0,edge.width=abs(E(Endo_1wk_up_net7)$width/3),vertex.color=V(Endo_1wk_up_net7)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Endo_1wk_up_net7)$label,vertex.label.color="black",vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,edge.color=c("#ff8989","#9fe899")[1+(E(Endo_1wk_up_net7)$weight>0)],vertex.size=degree(Endo_1wk_up_net7)/3,layout=layout_with_dh (Endo_1wk_up_net7),main="Endo_1wk_up_net7")
dev.off()
```

#  >>>>> Endo_3wk_up >>>>>>>

```{r}
## import SparCC and p_value matrix into R 
Endo_3wk_up_SparCC_cor<-read.csv("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Endo_3wk_up_OTU_cor_sparcc.csv",header = TRUE,row.names = 1)
Endo_3wk_up_SparCC_cor<-as.matrix(Endo_3wk_up_SparCC_cor)
dim(Endo_3wk_up_SparCC_cor)
Endo_3wk_up_SparCC_cor[1:5,1:5]

# Read in corresponding p_value dataframe
Endo_3wk_up_SparCC_pvalue<-read.csv("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Endo_3wk_up_OTU_two_sided_pvalue.csv",header = TRUE,row.names = 1)
Endo_3wk_up_SparCC_pvalue<-as.matrix(Endo_3wk_up_SparCC_pvalue)
dim(Endo_3wk_up_SparCC_pvalue)
Endo_3wk_up_SparCC_pvalue[1:5,1:5]
```


## Generate vertices information and add relative abundance information to vertices and color information


```{r}
r_Endo_3wk_up_phyloseq<-transform_sample_counts(Endo_3wk_up,function(x) x/sum(x))

Endo_3wk_up_vertex<-data.frame(OtuID=taxa_names(r_Endo_3wk_up_phyloseq),tax_table(r_Endo_3wk_up_phyloseq),size=taxa_sums(r_Endo_3wk_up_phyloseq)/length(sample_sums(r_Endo_3wk_up_phyloseq))) # Here size are the average relative abundance of each OTU across all samples within same treatment
dim(Endo_3wk_up_vertex)
Endo_3wk_up_vertex[1:5,]
Endo_3wk_up_vertex<-Endo_3wk_up_vertex[,c(1,3,7,8)]
head(Endo_3wk_up_vertex)

Bacteria_phylum_vs_color_checklist<-data.frame(phylum=c('Proteobacteria','Actinobacteria','Bacteria_unclassified','Acidobacteria','Nitrospirae','Chloroflexi','Candidatus_Saccharibacteria','Firmicutes','Bacteroidetes','Verrucomicrobia','candidate_division_WPS-1','Gemmatimonadetes','Planctomycetes','BRC1','Latescibacteria','Chlamydiae','Others','Patescibacteria','Rokubacteria','Cyanobacteria'),colors=c('#ff9082','#98d66d','#fcef5d','#c45cc0','#45b1f9','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#bc8c38','#bcba6d','#91749e','#b2798d','#343837','#b23301','#235931'))
phylum_colors<-as.character(Bacteria_phylum_vs_color_checklist$colors)
phylum_labels<-as.character(Bacteria_phylum_vs_color_checklist$phylum)


tiff("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Phylum_colour_pie.tiff",units = 'in',width=10,height = 8,res=300,compression = 'lzw')
pie(rep (1,17),col=phylum_colors,labels=paste(phylum_labels,phylum_colors,sep=":"),radius = 1.0)
dev.off()

Endo_3wk_up_vertex$color<-factor(Endo_3wk_up_vertex$Phylum,levels=phylum_labels,labels=phylum_colors)
dim(Endo_3wk_up_vertex)
Endo_3wk_up_vertex[1:10,]
# I want to add label variable to Endo_3wk_up_vertex
Endo_3wk_up_vertex<-join (Endo_3wk_up_vertex,meta,by="OtuID") 

# Originally, I used merge function, but it turned out to be problematic, because after merge it will The rows are by default lexicographically sorted on the common columns. To resolve this problem, I changed merge function to joint function.
dim(Endo_3wk_up_vertex)
Endo_3wk_up_vertex[1:5,]
```


## Creat links and have a look at edge adn vertex information

```{r}
# read SparCC correlation matrix and p_value matrix to igraph
Endo_3wk_up_net<-graph_from_adjacency_matrix(Endo_3wk_up_SparCC_cor,weighted = TRUE)
Endo_3wk_up_net  #IGRAPH 8cfb97f DNW- 500 250000 -- 
summary(V(Endo_3wk_up_net)$name==as.character(Endo_3wk_up_vertex$OtuID)) 
V(Endo_3wk_up_net)$name[1:5]
E(Endo_3wk_up_net)$weight[1:5] # here the weight are just the SparCC coefficience

# Add vertex information to the network
V(Endo_3wk_up_net)$vertex_color<-as.character(Endo_3wk_up_vertex$color)
V(Endo_3wk_up_net)$size<-Endo_3wk_up_vertex$size
V(Endo_3wk_up_net)$Phylum<-as.character(Endo_3wk_up_vertex$Phylum)
V(Endo_3wk_up_net)$label<-as.character(Endo_3wk_up_vertex$label)

# Add edge information (p_value) to the network
E(Endo_3wk_up_net)$p_value<-Endo_3wk_up_SparCC_pvalue
E(Endo_3wk_up_net)$p_value[1:10]

## Add edge width to the network
E(Endo_3wk_up_net)$width<-E(Endo_3wk_up_net)$weight*10
E(Endo_3wk_up_net)$width[1:20]
```


------------------------------------------------------------------------------------------------------------------------------------------
In order to simplify the network, including remove of loop and mutual edges. If using simplify function on Endo_3wk_up_net, this will loose p_value and edge color information (which indicate if the interaction are negarive or positive)

To solve this problem, I will write out the edge and vertices in dataframe format and read into these two dataframe and modify the edge information globally before use ``graph_from_data_frame`` to create network object.
-----------------------------------------------------------------------------------------------------------------------------------------


## write out edges and vertices in dataframe format, and edit on the network

```{r}

# -- generate edge dataframe from whole network---

Endo_3wk_up_edges<-igraph::as_data_frame(Endo_3wk_up_net,what = "edges")
#Endo_3wk_up_edges
dim(Endo_3wk_up_edges) # 250000 rows and 5 columns
head(Endo_3wk_up_edges)
# Here the weight indicate the strength of interaction (SparCC value), p_value is the significance of interactions nad width is formulated based on weight * 40.

# -- generate vertices dataframe from whole network---

Endo_3wk_up_vertices<-igraph::as_data_frame(Endo_3wk_up_net,what = "vertices")
dim(Endo_3wk_up_vertices)
head(Endo_3wk_up_vertices)

# remove edges that from and to the same nodes(selfloops)

Endo_3wk_up_edges<-Endo_3wk_up_edges[which(Endo_3wk_up_edges$from!=Endo_3wk_up_edges$to),]
dim(Endo_3wk_up_edges) # 1208900 rows and 5 columns
head(Endo_3wk_up_edges)

# Creat a new network using the above vertices and edge data frame
Endo_3wk_up_net1 <- graph_from_data_frame(d=Endo_3wk_up_edges, vertices=Endo_3wk_up_vertices, directed=T)#here if I used directed=FALSE, I will lost my p_value information as well as weight information. 
Endo_3wk_up_net1 # IGRAPH 526a39b DNW- 800 639200 --
E(Endo_3wk_up_net1)$weight[1:10]

# convert directed network to undirected network

Endo_3wk_up_net2<-as.undirected(Endo_3wk_up_net1 ,mode = 'mutual',edge.attr.comb = "mean") # as a result, weight, p_value and width of edges are all averaged for mutual edges.
Endo_3wk_up_net2 # IGRAPH 1bce66c UNW- 500 124750 - -
E(Endo_3wk_up_net2)$weight[1:10]
E(Endo_3wk_up_net2)$p_value[1:10]
E(Endo_3wk_up_net2)$width[1:10]


## simplify the network from Endo_3wk_up_net2

#1)  remove non_significant edges with alpha=0.05

summary(E(Endo_3wk_up_net2)$p_value<0.05) #TRUE 74374 and FALSE 530076   

Endo_3wk_up_net3<-delete.edges(Endo_3wk_up_net2,which(E(Endo_3wk_up_net2)$p_value>=0.05))
Endo_3wk_up_net3 # IGRAPH f25c247 UNW- 800 24383 -- 
summary(degree(Endo_3wk_up_net3)==0) # all of the nodes has connection to others
summary(E(Endo_3wk_up_net3)$weight>0) # 12668 positive interactions, 11715 nagative interactions

#2) remove non-significant edges with alpha=0.001

summary(E(Endo_3wk_up_net2)$p_value<0.001) # TRUE 14033, FALSE 590417
Endo_3wk_up_net4<-delete.edges(Endo_3wk_up_net2,which(E(Endo_3wk_up_net2)$p_value>=0.001))
Endo_3wk_up_net4 # IGRAPH b2b3bc7 UNW- 800 2481 --
Endo_3wk_up_net4<-delete.vertices(Endo_3wk_up_net4,degree(Endo_3wk_up_net4)==0) # IGRAPH 281b75b UNW- 789 2481 --

#3) further simplify the network based on degree
barplot(sort(degree(Endo_3wk_up_net4)),main = "Endo_3wk_up_net4_barplot")
summary(degree(Endo_3wk_up_net4)>15) # 36 nodes with degree larger than 100
Endo_3wk_up_net5<-delete.vertices(Endo_3wk_up_net4,which(degree(Endo_3wk_up_net4)<=5))
Endo_3wk_up_net5# IGRAPH ce4509d UNW- 36 536 --
summary(E(Endo_3wk_up_net5)$weight>0) # 261 positive and 275 negative

# top 50

Endo_3wk_up_net6<-delete.vertices(Endo_3wk_up_net4,labels(sort(degree(Endo_3wk_up_net4),FALSE)[c(1:(length(V(Endo_3wk_up_net4)$name)-50))])) # IGRAPH 56796e5 UNW- 50 937 --
Endo_3wk_up_net6<-delete.vertices(Endo_3wk_up_net6,degree(Endo_3wk_up_net6)==0) #IGRAPH cead187 UNW- 43 71 -- 
summary(E(Endo_3wk_up_net6)$weight>0) # logical : 46(negative)     25(positive)

# top100
Endo_3wk_up_net7<-delete.vertices(Endo_3wk_up_net4,labels(sort(degree(Endo_3wk_up_net4),FALSE)[c(1:(length(V(Endo_3wk_up_net4)$name)-100))]))
Endo_3wk_up_net7<-delete.vertices(Endo_3wk_up_net7,degree(Endo_3wk_up_net7)==0) #IGRAPH cc57451 UNW- 96 200 -
summary(E(Endo_3wk_up_net7)$weight>0) # 1273 (negative)    1228 (positive) 
```


## Plot network

```{r}

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Endo_3wk_up_net4.tiff', units="in", width=10, height=8, res=300)

plot(Endo_3wk_up_net4,edge.arrow.size=0,edge.width=abs(E(Endo_3wk_up_net4)$width/4),vertex.color=V(Endo_3wk_up_net4)$vertex_color,vertex.frame.color="#555555",vertex.label=NA,edge.color=c("#ff8989","#9fe899")[1+(E(Endo_3wk_up_net4)$weight>0)],vertex.size=degree(Endo_3wk_up_net4)/10,layout=layout_with_fr (Endo_3wk_up_net4),main="Endo_3wk_up_net4") # green edges means positive correlation and red edges mean negative correlation
dev.off()

#plot(Endo_3wk_up_net5,edge.arrow.size=0,edge.width=abs(E(Endo_3wk_up_net5)$width/3),vertex.color=V(Endo_3wk_up_net5)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Endo_3wk_up_net5)$label,edge.color=c("#ff8989","#9fe899")[1+(E(Endo_3wk_up_net5)$weight>0)],vertex.label.cex=.8,vertex.size=degree(Endo_3wk_up_net5)/2,layout=layout_with_dh (Endo_3wk_up_net5),main="Endo_3wk_up_net5")

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Endo_3wk_up_net6.tiff', units="in", width=10, height=8, res=300)
plot(Endo_3wk_up_net6,edge.arrow.size=0,edge.width=abs(E(Endo_3wk_up_net6)$width/3),vertex.color=V(Endo_3wk_up_net6)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Endo_3wk_up_net6)$label,vertex.label.color="black",vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,edge.color=c("#ff8989","#9fe899")[1+(E(Endo_3wk_up_net6)$weight>0)],vertex.size=degree(Endo_3wk_up_net6)/3,layout=layout_with_dh (Endo_3wk_up_net6),main="Endo_3wk_up_net6")
dev.off()

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Endo_3wk_up_net7.tiff', units="in", width=10, height=8, res=300)
plot(Endo_3wk_up_net7,edge.arrow.size=0,edge.width=abs(E(Endo_3wk_up_net7)$width/3),vertex.color=V(Endo_3wk_up_net7)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Endo_3wk_up_net7)$label,vertex.label.color="black",vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,edge.color=c("#ff8989","#9fe899")[1+(E(Endo_3wk_up_net7)$weight>0)],vertex.size=degree(Endo_3wk_up_net7)/3,layout=layout_with_dh (Endo_3wk_up_net7),main="Endo_3wk_up_net7")
dev.off()
```

#  >>>>> Endo_4wk_up >>>>>>>

```{r}
## import SparCC and p_value matrix into R 
Endo_4wk_up_SparCC_cor<-read.csv("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Endo_4wk_up_OTU_cor_sparcc.csv",header = TRUE,row.names = 1)
Endo_4wk_up_SparCC_cor<-as.matrix(Endo_4wk_up_SparCC_cor)
dim(Endo_4wk_up_SparCC_cor)
Endo_4wk_up_SparCC_cor[1:5,1:5]

# Read in corresponding p_value dataframe
Endo_4wk_up_SparCC_pvalue<-read.csv("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Endo_4wk_up_OTU_two_sided_pvalue.csv",header = TRUE,row.names = 1)
Endo_4wk_up_SparCC_pvalue<-as.matrix(Endo_4wk_up_SparCC_pvalue)
dim(Endo_4wk_up_SparCC_pvalue)
Endo_4wk_up_SparCC_pvalue[1:5,1:5]
```


## Generate vertices information and add relative abundance information to vertices and color information


```{r}
r_Endo_4wk_up_phyloseq<-transform_sample_counts(Endo_4wk_up,function(x) x/sum(x))

Endo_4wk_up_vertex<-data.frame(OtuID=taxa_names(r_Endo_4wk_up_phyloseq),tax_table(r_Endo_4wk_up_phyloseq),size=taxa_sums(r_Endo_4wk_up_phyloseq)/length(sample_sums(r_Endo_4wk_up_phyloseq))) # Here size are the average relative abundance of each OTU across all samples within same treatment
dim(Endo_4wk_up_vertex)
Endo_4wk_up_vertex[1:5,]
Endo_4wk_up_vertex<-Endo_4wk_up_vertex[,c(1,3,7,8)]
head(Endo_4wk_up_vertex)

Bacteria_phylum_vs_color_checklist<-data.frame(phylum=c('Proteobacteria','Actinobacteria','Bacteria_unclassified','Acidobacteria','Nitrospirae','Chloroflexi','Candidatus_Saccharibacteria','Firmicutes','Bacteroidetes','Verrucomicrobia','candidate_division_WPS-1','Gemmatimonadetes','Planctomycetes','BRC1','Latescibacteria','Chlamydiae','Others','Patescibacteria','Rokubacteria','Cyanobacteria'),colors=c('#ff9082','#98d66d','#fcef5d','#c45cc0','#45b1f9','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#bc8c38','#bcba6d','#91749e','#b2798d','#343837','#b23301','#235931'))
phylum_colors<-as.character(Bacteria_phylum_vs_color_checklist$colors)
phylum_labels<-as.character(Bacteria_phylum_vs_color_checklist$phylum)


tiff("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Phylum_colour_pie.tiff",units = 'in',width=10,height = 8,res=300,compression = 'lzw')
pie(rep (1,17),col=phylum_colors,labels=paste(phylum_labels,phylum_colors,sep=":"),radius = 1.0)
dev.off()

Endo_4wk_up_vertex$color<-factor(Endo_4wk_up_vertex$Phylum,levels=phylum_labels,labels=phylum_colors)
dim(Endo_4wk_up_vertex)
Endo_4wk_up_vertex[1:10,]
# I want to add label variable to Endo_4wk_up_vertex
Endo_4wk_up_vertex<-join (Endo_4wk_up_vertex,meta,by="OtuID") 

# Originally, I used merge function, but it turned out to be problematic, because after merge it will The rows are by default lexicographically sorted on the common columns. To resolve this problem, I changed merge function to joint function.
dim(Endo_4wk_up_vertex)
Endo_4wk_up_vertex[1:5,]
```


## Creat links and have a look at edge adn vertex information

```{r}
# read SparCC correlation matrix and p_value matrix to igraph
Endo_4wk_up_net<-graph_from_adjacency_matrix(Endo_4wk_up_SparCC_cor,weighted = TRUE)
Endo_4wk_up_net  #IGRAPH 8cfb97f DNW- 500 250000 -- 
summary(V(Endo_4wk_up_net)$name==as.character(Endo_4wk_up_vertex$OtuID)) 
V(Endo_4wk_up_net)$name[1:5]
E(Endo_4wk_up_net)$weight[1:5] # here the weight are just the SparCC coefficience

# Add vertex information to the network
V(Endo_4wk_up_net)$vertex_color<-as.character(Endo_4wk_up_vertex$color)
V(Endo_4wk_up_net)$size<-Endo_4wk_up_vertex$size
V(Endo_4wk_up_net)$Phylum<-as.character(Endo_4wk_up_vertex$Phylum)
V(Endo_4wk_up_net)$label<-as.character(Endo_4wk_up_vertex$label)

# Add edge information (p_value) to the network
E(Endo_4wk_up_net)$p_value<-Endo_4wk_up_SparCC_pvalue
E(Endo_4wk_up_net)$p_value[1:10]

## Add edge width to the network
E(Endo_4wk_up_net)$width<-E(Endo_4wk_up_net)$weight*10
E(Endo_4wk_up_net)$width[1:20]
```


------------------------------------------------------------------------------------------------------------------------------------------
In order to simplify the network, including remove of loop and mutual edges. If using simplify function on Endo_4wk_up_net, this will loose p_value and edge color information (which indicate if the interaction are negarive or positive)

To solve this problem, I will write out the edge and vertices in dataframe format and read into these two dataframe and modify the edge information globally before use ``graph_from_data_frame`` to create network object.
-----------------------------------------------------------------------------------------------------------------------------------------


## write out edges and vertices in dataframe format, and edit on the network

```{r}

# -- generate edge dataframe from whole network---

Endo_4wk_up_edges<-igraph::as_data_frame(Endo_4wk_up_net,what = "edges")
#Endo_4wk_up_edges
dim(Endo_4wk_up_edges) # 250000 rows and 5 columns
head(Endo_4wk_up_edges)
# Here the weight indicate the strength of interaction (SparCC value), p_value is the significance of interactions nad width is formulated based on weight * 40.

# -- generate vertices dataframe from whole network---

Endo_4wk_up_vertices<-igraph::as_data_frame(Endo_4wk_up_net,what = "vertices")
dim(Endo_4wk_up_vertices)
head(Endo_4wk_up_vertices)

# remove edges that from and to the same nodes(selfloops)

Endo_4wk_up_edges<-Endo_4wk_up_edges[which(Endo_4wk_up_edges$from!=Endo_4wk_up_edges$to),]
dim(Endo_4wk_up_edges) # 1208900 rows and 5 columns
head(Endo_4wk_up_edges)

# Creat a new network using the above vertices and edge data frame
Endo_4wk_up_net1 <- graph_from_data_frame(d=Endo_4wk_up_edges, vertices=Endo_4wk_up_vertices, directed=T)#here if I used directed=FALSE, I will lost my p_value information as well as weight information. 
Endo_4wk_up_net1 # IGRAPH 526a39b DNW- 800 639200 --
E(Endo_4wk_up_net1)$weight[1:10]

# convert directed network to undirected network

Endo_4wk_up_net2<-as.undirected(Endo_4wk_up_net1 ,mode = 'mutual',edge.attr.comb = "mean") # as a result, weight, p_value and width of edges are all averaged for mutual edges.
Endo_4wk_up_net2 # IGRAPH 1bce66c UNW- 500 124750 - -
E(Endo_4wk_up_net2)$weight[1:10]
E(Endo_4wk_up_net2)$p_value[1:10]
E(Endo_4wk_up_net2)$width[1:10]


## simplify the network from Endo_4wk_up_net2

#1)  remove non_significant edges with alpha=0.05

summary(E(Endo_4wk_up_net2)$p_value<0.05) #TRUE 74374 and FALSE 530076   

Endo_4wk_up_net3<-delete.edges(Endo_4wk_up_net2,which(E(Endo_4wk_up_net2)$p_value>=0.05))
Endo_4wk_up_net3 # IGRAPH f25c247 UNW- 800 24383 -- 
summary(degree(Endo_4wk_up_net3)==0) # all of the nodes has connection to others
summary(E(Endo_4wk_up_net3)$weight>0) # 12668 positive interactions, 11715 nagative interactions

#2) remove non-significant edges with alpha=0.001

summary(E(Endo_4wk_up_net2)$p_value<0.001) # TRUE 14033, FALSE 590417
Endo_4wk_up_net4<-delete.edges(Endo_4wk_up_net2,which(E(Endo_4wk_up_net2)$p_value>=0.001))
Endo_4wk_up_net4 # IGRAPH b2b3bc7 UNW- 800 2481 --
Endo_4wk_up_net4<-delete.vertices(Endo_4wk_up_net4,degree(Endo_4wk_up_net4)==0) # IGRAPH 281b75b UNW- 789 2481 --

#3) further simplify the network based on degree
barplot(sort(degree(Endo_4wk_up_net4)),main = "Endo_4wk_up_net4_barplot")
summary(degree(Endo_4wk_up_net4)>15) # 36 nodes with degree larger than 100
Endo_4wk_up_net5<-delete.vertices(Endo_4wk_up_net4,which(degree(Endo_4wk_up_net4)<=5))
Endo_4wk_up_net5# IGRAPH ce4509d UNW- 36 536 --
summary(E(Endo_4wk_up_net5)$weight>0) # 261 positive and 275 negative

# top 50

Endo_4wk_up_net6<-delete.vertices(Endo_4wk_up_net4,labels(sort(degree(Endo_4wk_up_net4),FALSE)[c(1:(length(V(Endo_4wk_up_net4)$name)-50))])) # IGRAPH 56796e5 UNW- 50 937 --
Endo_4wk_up_net6<-delete.vertices(Endo_4wk_up_net6,degree(Endo_4wk_up_net6)==0) #IGRAPH cead187 UNW- 43 71 -- 
summary(E(Endo_4wk_up_net6)$weight>0) # logical : 46(negative)     25(positive)

# top100
Endo_4wk_up_net7<-delete.vertices(Endo_4wk_up_net4,labels(sort(degree(Endo_4wk_up_net4),FALSE)[c(1:(length(V(Endo_4wk_up_net4)$name)-100))]))
Endo_4wk_up_net7<-delete.vertices(Endo_4wk_up_net7,degree(Endo_4wk_up_net7)==0) #IGRAPH cc57451 UNW- 96 200 -
summary(E(Endo_4wk_up_net7)$weight>0) # 1273 (negative)    1228 (positive) 
```

## Plot network

```{r}

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Endo_4wk_up_net4.tiff', units="in", width=10, height=8, res=300)

plot(Endo_4wk_up_net4,edge.arrow.size=0,edge.width=abs(E(Endo_4wk_up_net4)$width/4),vertex.color=V(Endo_4wk_up_net4)$vertex_color,vertex.frame.color="#555555",vertex.label=NA,edge.color=c("#ff8989","#9fe899")[1+(E(Endo_4wk_up_net4)$weight>0)],vertex.size=degree(Endo_4wk_up_net4)/10,layout=layout_with_fr (Endo_4wk_up_net4),main="Endo_4wk_up_net4") # green edges means positive correlation and red edges mean negative correlation
dev.off()

#plot(Endo_4wk_up_net5,edge.arrow.size=0,edge.width=abs(E(Endo_4wk_up_net5)$width/3),vertex.color=V(Endo_4wk_up_net5)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Endo_4wk_up_net5)$label,edge.color=c("#ff8989","#9fe899")[1+(E(Endo_4wk_up_net5)$weight>0)],vertex.label.cex=.8,vertex.size=degree(Endo_4wk_up_net5)/2,layout=layout_with_dh (Endo_4wk_up_net5),main="Endo_4wk_up_net5")

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Endo_4wk_up_net6.tiff', units="in", width=10, height=8, res=300)
plot(Endo_4wk_up_net6,edge.arrow.size=0,edge.width=abs(E(Endo_4wk_up_net6)$width/3),vertex.color=V(Endo_4wk_up_net6)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Endo_4wk_up_net6)$label,vertex.label.color="black",vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,edge.color=c("#ff8989","#9fe899")[1+(E(Endo_4wk_up_net6)$weight>0)],vertex.size=degree(Endo_4wk_up_net6)/3,layout=layout_with_dh (Endo_4wk_up_net6),main="Endo_4wk_up_net6")
dev.off()

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Endo_4wk_up_net7.tiff', units="in", width=10, height=8, res=300)
plot(Endo_4wk_up_net7,edge.arrow.size=0,edge.width=abs(E(Endo_4wk_up_net7)$width/3),vertex.color=V(Endo_4wk_up_net7)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Endo_4wk_up_net7)$label,vertex.label.color="black",vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,edge.color=c("#ff8989","#9fe899")[1+(E(Endo_4wk_up_net7)$weight>0)],vertex.size=degree(Endo_4wk_up_net7)/3,layout=layout_with_dh (Endo_4wk_up_net7),main="Endo_4wk_up_net7")
dev.off()


```


# >>>>> Seed_up >>>>>>>

```{r}
## import SparCC and p_value matrix into R 
Seed_up_SparCC_cor<-read.csv("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Seed_up_OTU_cor_sparcc.csv",header = TRUE,row.names = 1)
Seed_up_SparCC_cor<-as.matrix(Seed_up_SparCC_cor)
dim(Seed_up_SparCC_cor)
Seed_up_SparCC_cor[1:5,1:5]

# Read in corresponding p_value dataframe
Seed_up_SparCC_pvalue<-read.csv("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/Seed_up_OTU_two_sided_pvalue.csv",header = TRUE,row.names = 1)
Seed_up_SparCC_pvalue<-as.matrix(Seed_up_SparCC_pvalue)
dim(Seed_up_SparCC_pvalue)
Seed_up_SparCC_pvalue[1:5,1:5]
```


## Generate vertices information and add relative abundance information to vertices and color information


```{r}
r_Seed_up_phyloseq<-transform_sample_counts(Seed_up,function(x) x/sum(x))

Seed_up_vertex<-data.frame(OtuID=taxa_names(r_Seed_up_phyloseq),tax_table(r_Seed_up_phyloseq),size=taxa_sums(r_Seed_up_phyloseq)/length(sample_sums(r_Seed_up_phyloseq))) # Here size are the average relative abundance of each OTU across all samples within same treatment
dim(Seed_up_vertex)
Seed_up_vertex[1:5,]
Seed_up_vertex<-Seed_up_vertex[,c(1,3,7,8)]
head(Seed_up_vertex)

Seed_up_vertex$color<-factor(Seed_up_vertex$Phylum,levels=phylum_labels,labels=phylum_colors)
dim(Seed_up_vertex)
Seed_up_vertex[1:10,]
# I want to add label variable to Seed_up_vertex
Seed_up_vertex<-join (Seed_up_vertex,meta,by="OtuID") 

# Originally, I used merge function, but it turned out to be problematic, because after merge it will The rows are by default lexicographically sorted on the common columns. To resolve this problem, I changed merge function to joint function.
dim(Seed_up_vertex)
Seed_up_vertex[1:5,]
```


## Creat links and have a look at edge adn vertex information

```{r}
# read SparCC correlation matrix and p_value matrix to igraph
Seed_up_net<-graph_from_adjacency_matrix(Seed_up_SparCC_cor,weighted = TRUE)
Seed_up_net  #IGRAPH 8cfb97f DNW- 500 250000 -- 
summary(V(Seed_up_net)$name==as.character(Seed_up_vertex$OtuID)) 
V(Seed_up_net)$name[1:5]
E(Seed_up_net)$weight[1:5] # here the weight are just the SparCC coefficience

# Add vertex information to the network
V(Seed_up_net)$vertex_color<-as.character(Seed_up_vertex$color)
V(Seed_up_net)$size<-Seed_up_vertex$size
V(Seed_up_net)$Phylum<-as.character(Seed_up_vertex$Phylum)
V(Seed_up_net)$label<-as.character(Seed_up_vertex$label)

# Add edge information (p_value) to the network
E(Seed_up_net)$p_value<-Seed_up_SparCC_pvalue
E(Seed_up_net)$p_value[1:10]

## Add edge width to the network
E(Seed_up_net)$width<-E(Seed_up_net)$weight*10
E(Seed_up_net)$width[1:20]
```


------------------------------------------------------------------------------------------------------------------------------------------
In order to simplify the network, including remove of loop and mutual edges. If using simplify function on Seed_up_net, this will loose p_value and edge color information (which indicate if the interaction are negarive or positive)

To solve this problem, I will write out the edge and vertices in dataframe format and read into these two dataframe and modify the edge information globally before use ``graph_from_data_frame`` to create network object.
-----------------------------------------------------------------------------------------------------------------------------------------


## write out edges and vertices in dataframe format, and edit on the network

```{r}

# -- generate edge dataframe from whole network---

Seed_up_edges<-igraph::as_data_frame(Seed_up_net,what = "edges")
#Seed_up_edges
dim(Seed_up_edges) # 250000 rows and 5 columns
head(Seed_up_edges)
# Here the weight indicate the strength of interaction (SparCC value), p_value is the significance of interactions nad width is formulated based on weight * 40.

# -- generate vertices dataframe from whole network---

Seed_up_vertices<-igraph::as_data_frame(Seed_up_net,what = "vertices")
dim(Seed_up_vertices)
head(Seed_up_vertices)

# remove edges that from and to the same nodes(selfloops)

Seed_up_edges<-Seed_up_edges[which(Seed_up_edges$from!=Seed_up_edges$to),]
dim(Seed_up_edges) # 1208900 rows and 5 columns
head(Seed_up_edges)

# Creat a new network using the above vertices and edge data frame
Seed_up_net1 <- graph_from_data_frame(d=Seed_up_edges, vertices=Seed_up_vertices, directed=T)#here if I used directed=FALSE, I will lost my p_value information as well as weight information. 
Seed_up_net1 # IGRAPH 526a39b DNW- 800 639200 --
E(Seed_up_net1)$weight[1:10]

# convert directed network to undirected network

Seed_up_net2<-as.undirected(Seed_up_net1 ,mode = 'mutual',edge.attr.comb = "mean") # as a result, weight, p_value and width of edges are all averaged for mutual edges.
Seed_up_net2 # IGRAPH 1bce66c UNW- 500 124750 - -
E(Seed_up_net2)$weight[1:10]
E(Seed_up_net2)$p_value[1:10]
E(Seed_up_net2)$width[1:10]


## simplify the network from Seed_up_net2

#1)  remove non_significant edges with alpha=0.05

summary(E(Seed_up_net2)$p_value<0.05) #TRUE 74374 and FALSE 530076   

Seed_up_net3<-delete.edges(Seed_up_net2,which(E(Seed_up_net2)$p_value>=0.05))
Seed_up_net3 # IGRAPH f25c247 UNW- 800 24383 -- 
summary(degree(Seed_up_net3)==0) # all of the nodes has connection to others
summary(E(Seed_up_net3)$weight>0) # 12668 positive interactions, 11715 nagative interactions

#2) remove non-significant edges with alpha=0.001

summary(E(Seed_up_net2)$p_value<0.001) # TRUE 14033, FALSE 590417
Seed_up_net4<-delete.edges(Seed_up_net2,which(E(Seed_up_net2)$p_value>=0.001))
Seed_up_net4 # IGRAPH b2b3bc7 UNW- 800 2481 --
Seed_up_net4<-delete.vertices(Seed_up_net4,degree(Seed_up_net4)==0) # IGRAPH 281b75b UNW- 789 2481 --

#3) further simplify the network based on degree
barplot(sort(degree(Seed_up_net4)),main = "Seed_up_net4_barplot")
summary(degree(Seed_up_net4)>15) # 36 nodes with degree larger than 100
Seed_up_net5<-delete.vertices(Seed_up_net4,which(degree(Seed_up_net4)<=5))
Seed_up_net5# IGRAPH ce4509d UNW- 36 536 --
summary(E(Seed_up_net5)$weight>0) # 261 positive and 275 negative

# top 50

Seed_up_net6<-Seed_up_net5
summary(E(Seed_up_net6)$weight>0) # logical : 46(negative)     25(positive)

```


## Plot network

```{r}

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Seed_up_net4.tiff', units="in", width=10, height=8, res=300)

plot(Seed_up_net4,edge.arrow.size=0,edge.width=abs(E(Seed_up_net4)$width/4),vertex.color=V(Seed_up_net4)$vertex_color,vertex.frame.color="#555555",vertex.label=NA,edge.color=c("#ff8989","#9fe899")[1+(E(Seed_up_net4)$weight>0)],vertex.size=degree(Seed_up_net4)/10,layout=layout_with_fr (Seed_up_net4),main="Seed_up_net4") # green edges means positive correlation and red edges mean negative correlation
dev.off()

#plot(Seed_up_net5,edge.arrow.size=0,edge.width=abs(E(Seed_up_net5)$width/3),vertex.color=V(Seed_up_net5)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Seed_up_net5)$label,edge.color=c("#ff8989","#9fe899")[1+(E(Seed_up_net5)$weight>0)],vertex.label.cex=.8,vertex.size=degree(Seed_up_net5)/2,layout=layout_with_dh (Seed_up_net5),main="Seed_up_net5")

tiff('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Seed_up_net6.tiff', units="in", width=10, height=8, res=300)
plot(Seed_up_net6,edge.arrow.size=0,edge.width=abs(E(Seed_up_net6)$width/3),vertex.color=V(Seed_up_net6)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Seed_up_net6)$label,vertex.label.color="black",vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,edge.color=c("#ff8989","#9fe899")[1+(E(Seed_up_net6)$weight>0)],vertex.size=degree(Seed_up_net6)/3,layout=layout_with_dh (Seed_up_net6),main="Seed_up_net6")
dev.off()

```

## United big net6

```{r}
Treatment_checklist<-data.frame(Treatment_colors=c('#9b9696','#ffee32','#b6cc0e','#ed5567','#07aeba','#3a44ff','#f936f6','#723434','#316022','#05fff2','#b105fc'),Treatment_labels=c('Soil','Seed','Bulk_1wk','Bulk_3wk','Bulk_4wk','Rhi_1wk','Rhi_3wk','Rhi_4wk','Endo_1wk','Endo_3wk','Endo_4wk'))
Treatment_colors<-as.character(Treatment_checklist$Treatment_colors)
Treatment_labels<-as.character(Treatment_checklist$Treatment_labels)
tiff("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Treatment_colour_pie.tiff",units = 'in',width=10,height = 8,res=300,compression = 'lzw')
pie(rep (1,11),col=Treatment_colors,labels=paste(Treatment_labels,Treatment_colors,sep=":"),radius = 1.0)
dev.off()

# list all subnetwork
Soil_up_net6
Seed_up_net6
Bulk_1wk_up_net6
Bulk_3wk_up_net6
Bulk_4wk_up_net6
Rhi_1wk_up_net6
Rhi_1wk_up_net6
Rhi_1wk_up_net6
Endo_1wk_up_net6
Endo_3wk_up_net6
Endo_4wk_up_net6

## Define edge color for the big network
E(Soil_up_net6)$color<-'#9b9696'
E(Seed_up_net6)$color<-'#ffee32'
E(Bulk_1wk_up_net6)$color<-'#b6cc0e'
E(Bulk_3wk_up_net6)$color<-'#ed5567'
E(Bulk_4wk_up_net6)$color<-'#07aeba'
E(Rhi_1wk_up_net6)$color<-'#3a44ff'
E(Rhi_3wk_up_net6)$color<-'#f936f6'
E(Rhi_4wk_up_net6)$color<-'#723434'
E(Endo_1wk_up_net6)$color<-'#316022'
E(Endo_3wk_up_net6)$color<-'#05fff2'
E(Endo_4wk_up_net6)$color<-'#b105fc'


big_net6<-union(Soil_up_net6,Seed_up_net6,Bulk_1wk_up_net6,Bulk_3wk_up_net6,Bulk_4wk_up_net6,Rhi_1wk_up_net6,Rhi_3wk_up_net6,Rhi_4wk_up_net6,Endo_1wk_up_net6,Endo_3wk_up_net6,Endo_4wk_up_net6,byname=TRUE)

str(edge_attr(big_net6))

big_net6_edge_col<-data.frame(Soil_col=edge_attr(big_net6,"color_1"),Seed_col=edge_attr(big_net6,"color_2"),Bulk_1wk_col=edge_attr(big_net6,"color_3"),Bulk_3wk_col=edge_attr(big_net6,"color_4"),Bulk_4wk_col=edge_attr(big_net6,"color_5"),Rhi_1wk_col=edge_attr(big_net6,"color_6"),Rhi_3wk_col=edge_attr(big_net6,"color_7"),Rhi_4wk_col=edge_attr(big_net6,"color_8"),Endo_1wk_col=edge_attr(big_net6,"color_9"),Endo_3wk_col=edge_attr(big_net6,"color_10"),Endo_4wk_col=edge_attr(big_net6,"color_11"))

#write.csv(big_net6_edge_col,file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/big_net6_edge_col.csv")

big_net6_edge_col2<-read.csv(file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/big_net6_edge_col_edit.csv",header = TRUE,row.names = 1)
head(big_net6_edge_col2)
E(big_net6)$color<-as.character (big_net6_edge_col2$color)

big_net6_vertex<-data.frame(OtuID=V(big_net6)$name,random=rep("whatever",312))
big_net6_vertex<-join(big_net6_vertex,meta,by="OtuID")
head(big_net6_vertex)
summary(big_net6_vertex$Phylum)
sort(summary(big_net6_vertex$Phylum),decreasing = TRUE)
big_net6_vertex$color<-factor(big_net6_vertex$Phylum,levels = phylum_labels,labels = phylum_colors)
head(big_net6_vertex)
identical( as.factor(V(big_net6)$name),big_net6_vertex$OtuID)
V(big_net6)$color<-as.character( big_net6_vertex$color)  
V(big_net6)$name<-as.character(V(big_net6)$name)
network_nodes_checkup<-read.csv("network_nodes_checkup.csv",header = TRUE) # This csv file is created from taxonomy file by saving it in csv format.
net6_degree_df<-data.frame(degree=sort(degree(big_net6),TRUE))
net6_degree_df<-data.frame(net6_degree_df,OTU=row.names(net6_degree_df))
big_net6_degree_sorted_OTUs<-merge(net6_degree_df,network_nodes_checkup,by="OTU")
big_net6_degree_sorted_OTUs<-big_net6_degree_sorted_OTUs[order(big_net6_degree_sorted_OTUs$degree,decreasing = TRUE),]
big_net6_degree_sorted_OTUs_up<-data.frame(big_net6_degree_sorted_OTUs,OtuID=big_net6_degree_sorted_OTUs$OTU)
big_net6_degree_sorted_OTUs_up2<-join(big_net6_degree_sorted_OTUs_up,meta,by="OtuID")

net6_phylum_distribution<-big_net6_degree_sorted_OTUs_up2 %>%
  plyr::count("Phylum")
net6_phylum_distribution_sort<-net6_phylum_distribution[order(net6_phylum_distribution$freq,decreasing = TRUE),]
net6_phylum_distribution_sort


tiff ('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/big_net6_fr_nontreatment_labels_no_tax_labels.tiff', units="in", width=20, height=15, res=300)
plot(big_net6,layout=layout_with_fr(big_net6),vertex.size=2,vertex.label=V(big_net6)$name,vertex.label.cex=.3,vertex.color=V(big_net6)$color,edge.color=E(big_net6)$color,edge.width=2)
#legend(x=-1.3,y=1,c("Soil","SEED","Bulk_1wk","Bulk_3wk","Bulk_4wk","Rhi_1wk","Rhi_3wk","Rhi_4wk","Endo_1wk","Endo_3wk","Endo_4wk"),lty=1,lwd = 4, col = c('#9b9696','#ffee32','#b6cc0e','#ed5567','#07aeba','#3a44ff','#f936f6','#723434','#316022','#05fff2','#b105fc'),text.col=c('#9b9696','#ffee32','#b6cc0e','#ed5567','#07aeba','#3a44ff','#f936f6','#723434','#316022','#05fff2','#b105fc'),pt.cex = 2,ncol = 1,cex = 3,bty = "n")
#legend(x=-1.0,y=-0.95,c('Proteobacteria','Actinobacteria','Bacteria_unclassified','Acidobacteria','Nitrospirae','Chloroflexi','Candidatus_Saccharibacteria','Firmicutes','Bacteroidetes','Verrucomicrobia','candidate_division_WPS-1','Gemmatimonadetes','Planctomycetes','BRC1','Latescibacteria','Chlamydiae','Others','Patescibacteria','Rokubacteria','Cyanobacteria'),pch = 19, col = c('#ff9082','#98d66d','#fcef5d','#c45cc0','#45b1f9','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#bc8c38','#bcba6d','#91749e','#b2798d','#343837','#b23301','#235931'),ncol = 5,cex = 2,bty = "o")
dev.off()

tiff ('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/big_net6_fr.tiff', units="in", width=30, height=15, res=300)
plot(big_net6,layout=layout_with_fr(big_net6),vertex.size=2,vertex.label=V(big_net6)$name,vertex.label.cex=.5,vertex.color=V(big_net6)$color,edge.color=E(big_net6)$color,edge.width=2)
legend(x=-2,y=1,c("Soil","SEED","Bulk_1wk","Bulk_3wk","Bulk_4wk","Rhi_1wk","Rhi_3wk","Rhi_4wk","Endo_1wk","Endo_3wk","Endo_4wk"),lty=1,lwd = 4, col = c('#9b9696','#ffee32','#b6cc0e','#ed5567','#07aeba','#3a44ff','#f936f6','#723434','#316022','#05fff2','#b105fc'),text.col=c('#9b9696','#ffee32','#b6cc0e','#ed5567','#07aeba','#3a44ff','#f936f6','#723434','#316022','#05fff2','#b105fc'),pt.cex = 2,ncol = 1,cex = 3,bty = "n")
legend(x=-2,y=-0.85,c('Proteobacteria','Actinobacteria','Bacteria_unclassified','Acidobacteria','Nitrospirae','Chloroflexi','Candidatus_Saccharibacteria','Firmicutes','Bacteroidetes','Verrucomicrobia','candidate_division_WPS-1','Gemmatimonadetes','Planctomycetes','BRC1','Latescibacteria','Chlamydiae','Others','Patescibacteria','Rokubacteria','Cyanobacteria'),pch = 19, col = c('#ff9082','#98d66d','#fcef5d','#c45cc0','#45b1f9','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#bc8c38','#bcba6d','#91749e','#b2798d','#343837','#b23301','#235931'),ncol = 5,cex = 2,bty = "o")
dev.off()


tiff ('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/big_net6_fr_update_node_size.tiff', units="in", width=30, height=20, res=200)
plot(big_net6,layout=layout_with_fr(big_net6),vertex.size=degree(big_net6)/16,vertex.label=V(big_net6)$name,vertex.label.cex=.5,vertex.color=V(big_net6)$color,edge.color=E(big_net6)$color,edge.width=2)
legend(x=-1.3,y=1,c("Soil","SEED","Bulk_1wk","Bulk_3wk","Bulk_4wk","Rhi_1wk","Rhi_3wk","Rhi_4wk","Endo_1wk","Endo_3wk","Endo_4wk"),lty=1,lwd = 4, col = c('#9b9696','#ffee32','#b6cc0e','#ed5567','#07aeba','#3a44ff','#f936f6','#723434','#316022','#05fff2','#b105fc'),text.col=c('#9b9696','#ffee32','#b6cc0e','#ed5567','#07aeba','#3a44ff','#f936f6','#723434','#316022','#05fff2','#b105fc'),pt.cex = 2,ncol = 1,cex = 3,bty = "n")
legend(x=-1.5,y=-0.95,c('Proteobacteria','Actinobacteria','Bacteria_unclassified','Acidobacteria','Nitrospirae','Chloroflexi','Candidatus_Saccharibacteria','Firmicutes','Bacteroidetes','Verrucomicrobia','candidate_division_WPS-1','Gemmatimonadetes','Planctomycetes','BRC1','Latescibacteria','Chlamydiae','Others','Patescibacteria','Rokubacteria','Cyanobacteria'),pch = 19, col = c('#ff9082','#98d66d','#fcef5d','#c45cc0','#45b1f9','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#bc8c38','#bcba6d','#91749e','#b2798d','#343837','#b23301','#235931'),ncol = 5,cex = 2,bty = "o")
dev.off()

tiff ('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/big_net6_fr_update_node_size_no_labels.tiff', units="in", width=20, height=20, res=200)
plot(big_net6,layout=layout_with_fr(big_net6),vertex.size=degree(big_net6)/16,vertex.label=V(big_net6)$name,vertex.label.cex=.5,vertex.color=V(big_net6)$color,edge.color=E(big_net6)$color,edge.width=2)
#legend(x=-1.3,y=1,c("Soil","SEED","Bulk_1wk","Bulk_3wk","Bulk_4wk","Rhi_1wk","Rhi_3wk","Rhi_4wk","Endo_1wk","Endo_3wk","Endo_4wk"),lty=1,lwd = 4, col = c('#9b9696','#ffee32','#b6cc0e','#ed5567','#07aeba','#3a44ff','#f936f6','#723434','#316022','#05fff2','#b105fc'),text.col=c('#9b9696','#ffee32','#b6cc0e','#ed5567','#07aeba','#3a44ff','#f936f6','#723434','#316022','#05fff2','#b105fc'),pt.cex = 2,ncol = 1,cex = 3,bty = "n")
#legend(x=-1.0,y=-0.95,c('Proteobacteria','Actinobacteria','Bacteria_unclassified','Acidobacteria','Nitrospirae','Chloroflexi','Candidatus_Saccharibacteria','Firmicutes','Bacteroidetes','Verrucomicrobia','candidate_division_WPS-1','Gemmatimonadetes','Planctomycetes','BRC1','Latescibacteria','Chlamydiae','Others','Patescibacteria','Rokubacteria','Cyanobacteria'),pch = 19, col = c('#ff9082','#98d66d','#fcef5d','#c45cc0','#45b1f9','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#bc8c38','#bcba6d','#91749e','#b2798d','#343837','#b23301','#235931'),ncol = 5,cex = 2,bty = "o")
dev.off()

```


## United big net4

```{r}
Treatment_checklist<-data.frame(Treatment_colors=c('#9b9696','#ffee32','#b6cc0e','#ed5567','#07aeba','#3a44ff','#f936f6','#723434','#316022','#05fff2','#b105fc'),Treatment_labels=c('Soil','SEED','Bulk_1wk','Bulk_3wk','Bulk_4wk','Rhi_1wk','Rhi_3wk','Rhi_4wk','Endo_1wk','Endo_3wk','Endo_4wk'))
Treatment_colors<-as.character(Treatment_checklist$Treatment_colors)
Treatment_labels<-as.character(Treatment_checklist$Treatment_labels)
tiff("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/Treatment_colour_pie.tiff",units = 'in',width=10,height = 8,res=300,compression = 'lzw')
pie(rep (1,11),col=Treatment_colors,labels=paste(Treatment_labels,Treatment_colors,sep=":"),radius = 1.0)
dev.off()

# list all subnetwork
Soil_up_net4
Seed_up_net4
Bulk_1wk_up_net4
Bulk_3wk_up_net4
Bulk_4wk_up_net4
Rhi_1wk_up_net4
Rhi_1wk_up_net4
Rhi_1wk_up_net4
Endo_1wk_up_net4
Endo_3wk_up_net4
Endo_4wk_up_net4

## Define edge color for the big network
E(Soil_up_net4)$color<-'#9b9696'
E(Seed_up_net4)$color<-'#ffee32'
E(Bulk_1wk_up_net4)$color<-'#b6cc0e'
E(Bulk_3wk_up_net4)$color<-'#ed5567'
E(Bulk_4wk_up_net4)$color<-'#07aeba'
E(Rhi_1wk_up_net4)$color<-'#3a44ff'
E(Rhi_3wk_up_net4)$color<-'#f936f6'
E(Rhi_4wk_up_net4)$color<-'#723434'
E(Endo_1wk_up_net4)$color<-'#316022'
E(Endo_3wk_up_net4)$color<-'#05fff2'
E(Endo_4wk_up_net4)$color<-'#b105fc'


big_net4<-union(Soil_up_net4,Seed_up_net4,Bulk_1wk_up_net4,Bulk_3wk_up_net4,Bulk_4wk_up_net4,Rhi_1wk_up_net4,Rhi_3wk_up_net4,Rhi_4wk_up_net4,Endo_1wk_up_net4,Endo_3wk_up_net4,Endo_4wk_up_net4,byname=TRUE)

str(edge_attr(big_net4))

big_net4_edge_col<-data.frame(Soil_col=edge_attr(big_net4,"color_1"),Seed_col=edge_attr(big_net4,"color_2"),Bulk_1wk_col=edge_attr(big_net4,"color_3"),Bulk_3wk_col=edge_attr(big_net4,"color_4"),Bulk_4wk_col=edge_attr(big_net4,"color_5"),Rhi_1wk_col=edge_attr(big_net4,"color_6"),Rhi_3wk_col=edge_attr(big_net4,"color_7"),Rhi_4wk_col=edge_attr(big_net4,"color_8"),Endo_1wk_col=edge_attr(big_net4,"color_9"),Endo_3wk_col=edge_attr(big_net4,"color_10"),Endo_4wk_col=edge_attr(big_net4,"color_11"))

#write.csv(big_net4_edge_col,file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/big_net4_edge_col.csv")

big_net4_edge_col2<-read.csv(file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/big_net4_edge_col_edit.csv",header = TRUE,row.names = 1)
head(big_net4_edge_col2)
E(big_net4)$color<-as.character(big_net4_edge_col2$color)

big_net4_vertex<-data.frame(OtuID=V(big_net4)$name,random=rep("whatever",1281))
big_net4_vertex<-join(big_net4_vertex,meta,by="OtuID")
head(big_net4_vertex)
summary(big_net4_vertex$Phylum)
sort(summary(big_net4_vertex$Phylum),decreasing = TRUE)
big_net4_vertex$color<-factor(big_net4_vertex$Phylum,levels = phylum_labels,labels = phylum_colors)
head(big_net4_vertex)
identical( as.factor(V(big_net4)$name),big_net4_vertex$OtuID)
V(big_net4)$color<-as.character( big_net4_vertex$color)  
V(big_net4)$name<-as.character(V(big_net4)$name)
network_nodes_checkup<-read.csv("network_nodes_checkup.csv",header = TRUE) # This csv file is created from taxonomy file by saving it in csv format.
net4_degree_df<-data.frame(degree=sort(degree(big_net4),TRUE))
net4_degree_df<-data.frame(net4_degree_df,OTU=row.names(net4_degree_df))
big_net4_degree_sorted_OTUs<-merge(net4_degree_df,network_nodes_checkup,by="OTU")
big_net4_degree_sorted_OTUs<-big_net4_degree_sorted_OTUs[order(big_net4_degree_sorted_OTUs$degree,decreasing = TRUE),]
big_net4_degree_sorted_OTUs$Genus[1:25]

big_net4
big_net4_simp<-delete.vertices(big_net4,degree(big_net4)<2)

tiff ('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/big_net4_simp_fr_nontreatment_labels_no_tax_labels.tiff', units="in", width=25, height=20, res=300)
plot(big_net4,layout=layout_with_fr(big_net4_simp),vertex.size=2,vertex.label=V(big_net4)$name,vertex.label.cex=.5,vertex.color=V(big_net4)$color,edge.color=E(big_net4)$color,edge.width=2)
#legend(x=-1.3,y=1,c("Soil","SEED","Bulk_1wk","Bulk_3wk","Bulk_4wk","Rhi_1wk","Rhi_3wk","Rhi_4wk","Endo_1wk","Endo_3wk","Endo_4wk"),lty=1,lwd = 4, col = c('#9b9696','#ffee32','#b6cc0e','#ed5567','#07aeba','#3a44ff','#f936f6','#723434','#316022','#05fff2','#b105fc'),text.col=c('#9b9696','#ffee32','#b6cc0e','#ed5567','#07aeba','#3a44ff','#f936f6','#723434','#316022','#05fff2','#b105fc'),pt.cex = 2,ncol = 1,cex = 3,bty = "n")
#legend(x=-1.0,y=-0.95,c('Proteobacteria','Actinobacteria','Bacteria_unclassified','Acidobacteria','Nitrospirae','Chloroflexi','Candidatus_Saccharibacteria','Firmicutes','Bacteroidetes','Verrucomicrobia','candidate_division_WPS-1','Gemmatimonadetes','Planctomycetes','BRC1','Latescibacteria','Chlamydiae','Others','Patescibacteria','Rokubacteria','Cyanobacteria'),pch = 19, col = c('#ff9082','#98d66d','#fcef5d','#c45cc0','#45b1f9','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#bc8c38','#bcba6d','#91749e','#b2798d','#343837','#b23301','#235931'),ncol = 5,cex = 2,bty = "o")
dev.off()

tiff ('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/big_net4_fr.tiff', units="in", width=25, height=20, res=300)
plot(big_net4,layout=layout_with_fr(big_net4),vertex.size=2,vertex.label=V(big_net4)$name,vertex.label.cex=.5,vertex.color=V(big_net4)$color,edge.color=E(big_net4)$color,edge.width=2)
legend(x=-1.3,y=1,c("Soil","SEED","Bulk_1wk","Bulk_3wk","Bulk_4wk","Rhi_1wk","Rhi_3wk","Rhi_4wk","Endo_1wk","Endo_3wk","Endo_4wk"),lty=1,lwd = 4, col = c('#9b9696','#ffee32','#b6cc0e','#ed5567','#07aeba','#3a44ff','#f936f6','#723434','#316022','#05fff2','#b105fc'),text.col=c('#9b9696','#ffee32','#b6cc0e','#ed5567','#07aeba','#3a44ff','#f936f6','#723434','#316022','#05fff2','#b105fc'),pt.cex = 2,ncol = 1,cex = 3,bty = "n")
legend(x=-1.5,y=-0.95,c('Proteobacteria','Actinobacteria','Bacteria_unclassified','Acidobacteria','Nitrospirae','Chloroflexi','Candidatus_Saccharibacteria','Firmicutes','Bacteroidetes','Verrucomicrobia','candidate_division_WPS-1','Gemmatimonadetes','Planctomycetes','BRC1','Latescibacteria','Chlamydiae','Others','Patescibacteria','Rokubacteria','Cyanobacteria'),pch = 19, col = c('#ff9082','#98d66d','#fcef5d','#c45cc0','#45b1f9','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#bc8c38','#bcba6d','#91749e','#b2798d','#343837','#b23301','#235931'),ncol = 5,cex = 2,bty = "o")
dev.off()


tiff ('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/big_net4_fr_update_node_size.tiff', units="in", width=25, height=20, res=200)
plot(big_net4,layout=layout_with_fr(big_net4),vertex.size=degree(big_net4)/8,vertex.label=V(big_net4)$name,vertex.label.cex=.5,vertex.color=V(big_net4)$color,edge.color=E(big_net4)$color,edge.width=2)
legend(x=-1.3,y=1,c("Soil","SEED","Bulk_1wk","Bulk_3wk","Bulk_4wk","Rhi_1wk","Rhi_3wk","Rhi_4wk","Endo_1wk","Endo_3wk","Endo_4wk"),lty=1,lwd = 4, col = c('#9b9696','#ffee32','#b6cc0e','#ed5567','#07aeba','#3a44ff','#f936f6','#723434','#316022','#05fff2','#b105fc'),text.col=c('#9b9696','#ffee32','#b6cc0e','#ed5567','#07aeba','#3a44ff','#f936f6','#723434','#316022','#05fff2','#b105fc'),pt.cex = 2,ncol = 1,cex = 3,bty = "n")
legend(x=-1.5,y=-0.95,c('Proteobacteria','Actinobacteria','Bacteria_unclassified','Acidobacteria','Nitrospirae','Chloroflexi','Candidatus_Saccharibacteria','Firmicutes','Bacteroidetes','Verrucomicrobia','candidate_division_WPS-1','Gemmatimonadetes','Planctomycetes','BRC1','Latescibacteria','Chlamydiae','Others','Patescibacteria','Rokubacteria','Cyanobacteria'),pch = 19, col = c('#ff9082','#98d66d','#fcef5d','#c45cc0','#45b1f9','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#bc8c38','#bcba6d','#91749e','#b2798d','#343837','#b23301','#235931'),ncol = 5,cex = 2,bty = "o")
dev.off()

tiff ('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Silva_tax_based_results/rarefaction_13021/SparCC/R_Output/big_net4_fr_update_node_size_no_labels.tiff', units="in", width=25, height=20, res=200)
plot(big_net4,layout=layout_in_circle(big_net4),vertex.size=degree(big_net4)/24,vertex.label=V(big_net4)$name,vertex.label.cex=.5,vertex.color=V(big_net4)$color,edge.color=E(big_net4)$color,edge.width=0.5)
#legend(x=-1.3,y=1,c("Soil","SEED","Bulk_1wk","Bulk_3wk","Bulk_4wk","Rhi_1wk","Rhi_3wk","Rhi_4wk","Endo_1wk","Endo_3wk","Endo_4wk"),lty=1,lwd = 4, col = c('#9b9696','#ffee32','#b6cc0e','#ed5567','#07aeba','#3a44ff','#f936f6','#723434','#316022','#05fff2','#b105fc'),text.col=c('#9b9696','#ffee32','#b6cc0e','#ed5567','#07aeba','#3a44ff','#f936f6','#723434','#316022','#05fff2','#b105fc'),pt.cex = 2,ncol = 1,cex = 3,bty = "n")
#legend(x=-1.0,y=-0.95,c('Proteobacteria','Actinobacteria','Bacteria_unclassified','Acidobacteria','Nitrospirae','Chloroflexi','Candidatus_Saccharibacteria','Firmicutes','Bacteroidetes','Verrucomicrobia','candidate_division_WPS-1','Gemmatimonadetes','Planctomycetes','BRC1','Latescibacteria','Chlamydiae','Others','Patescibacteria','Rokubacteria','Cyanobacteria'),pch = 19, col = c('#ff9082','#98d66d','#fcef5d','#c45cc0','#45b1f9','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#bc8c38','#bcba6d','#91749e','#b2798d','#343837','#b23301','#235931'),ncol = 5,cex = 2,bty = "o")
dev.off()

#plot(big_net4,layout=layout_with_gem (big_net4),vertex.size=degree(big_net4)/24,vertex.label=V(big_net4)$name,vertex.label.cex=.5,vertex.color=V(big_net4)$color,edge.color=E(big_net4)$color,edge.width=0.5)
```


## Extract neighbor vertices and calculate positive and negative correlations

```{r}
Bradyrhizobium_neighbors<-ego(Endo_1wk_up_net4,order = 2,mode = 'all',nodes = "Otu000001")
Bradyrhizobium_net_Endo_1wk_net4<-induced_subgraph(Endo_1wk_up_net4,vids = unlist(Bradyrhizobium_neighbors))
plot.igraph(Bradyrhizobium_net_Endo_1wk_net4,layout=layout_with_gem(Bradyrhizobium_net_Endo_1wk_net4))
edge_attr(Bradyrhizobium_net_Endo_1wk_net4)
plot(Bradyrhizobium_net_Endo_1wk_net4,edge.arrow.size=0,edge.width=abs(E(Bradyrhizobium_net_Endo_1wk_net4)$width/3),vertex.color=V(Bradyrhizobium_net_Endo_1wk_net4)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Bradyrhizobium_net_Endo_1wk_net4)$label,vertex.label.color="black",vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,edge.color=c("#ff8989","#9fe899")[1+(E(Bradyrhizobium_net_Endo_1wk_net4)$weight>0)],vertex.size=degree(Bradyrhizobium_net_Endo_1wk_net4)/3,layout=layout_with_dh (Bradyrhizobium_net_Endo_1wk_net4),main="Bradyrhizobium_net_Endo_1wk_net4")


Bradyrhizobium_neighbors<-ego(Endo_3wk_up_net4,order = 2,mode = 'all',nodes = "Otu000001")
Bradyrhizobium_net_Endo_3wk_net4<-induced_subgraph(Endo_3wk_up_net4,vids = unlist(Bradyrhizobium_neighbors))
plot.igraph(Bradyrhizobium_net_Endo_3wk_net4,layout=layout_with_gem(Bradyrhizobium_net_Endo_3wk_net4))
edge_attr(Bradyrhizobium_net_Endo_3wk_net4)
plot(Bradyrhizobium_net_Endo_3wk_net4,edge.arrow.size=0,edge.width=abs(E(Bradyrhizobium_net_Endo_3wk_net4)$width/3),vertex.color=V(Bradyrhizobium_net_Endo_3wk_net4)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Bradyrhizobium_net_Endo_3wk_net4)$label,vertex.label.color="black",vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,edge.color=c("#ff8989","#9fe899")[1+(E(Bradyrhizobium_net_Endo_3wk_net4)$weight>0)],vertex.size=degree(Bradyrhizobium_net_Endo_3wk_net4)/3,layout=layout_with_dh (Bradyrhizobium_net_Endo_3wk_net4),main="Bradyrhizobium_net_Endo_3wk_net4")


Bradyrhizobium_neighbors<-ego(Endo_4wk_up_net4,order = 2,mode = 'all',nodes = "Otu000001")
Bradyrhizobium_net_Endo_4wk_net4<-induced_subgraph(Endo_4wk_up_net4,vids = unlist(Bradyrhizobium_neighbors))
plot.igraph(Bradyrhizobium_net_Endo_4wk_net4,layout=layout_with_gem(Bradyrhizobium_net_Endo_4wk_net4))
edge_attr(Bradyrhizobium_net_Endo_4wk_net4)
plot(Bradyrhizobium_net_Endo_4wk_net4,edge.arrow.size=0,edge.width=abs(E(Bradyrhizobium_net_Endo_4wk_net4)$width/3),vertex.color=V(Bradyrhizobium_net_Endo_4wk_net4)$vertex_color,vertex.frame.color="#555555",vertex.label=V(Bradyrhizobium_net_Endo_4wk_net4)$label,vertex.label.color="black",vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,edge.color=c("#ff8989","#9fe899")[1+(E(Bradyrhizobium_net_Endo_4wk_net4)$weight>0)],vertex.size=degree(Bradyrhizobium_net_Endo_4wk_net4)/3,layout=layout_with_dh (Bradyrhizobium_net_Endo_4wk_net4),main="Bradyrhizobium_net_Endo_4wk_net4")


```

## Calculate topological characters of subnetwork net6

```{r}

## ---- Soil_net6 ------

Soil_up_net6
Soil_up_net6_ratio<-length(which(E(Soil_up_net6)$weight>0))/(length(E(Soil_up_net6)$weight)) #53.19%
Soil_up_net6_density<-edge_density(Soil_up_net6,loops = F) # 29.47%
Soil_up_net6_ave_degree<-sum(degree(Soil_up_net6,mode="all"))/length(V(Soil_up_net6)$name) #14.44
Soil_up_net6_centr_degree<-centr_degree(Soil_up_net6,mode="all",normalized = TRUE)$centralization #0.24
Soil_up_net6_closeness<-sum(closeness(Soil_up_net6,mode = "all",weights=rep(1,361),normalized = TRUE))/length(V(Soil_up_net6)$name) #0.2740913
Soil_up_net6_centr_closeness<-centr_clo(Soil_up_net6,mode = "all",normalized = TRUE)$centralization #0.26
Soil_up_net6_betweeness<-sum(betweenness(Soil_up_net6,directed = FALSE,weight=rep(1,361)))/length(V(Soil_up_net6)$name) #20.3
Soil_up_net6_centr_betweenness<-centr_betw(Soil_up_net6,directed = FALSE)$centralization #0.05
Soil_up_net6_modularity<-modularity(Soil_up_net6,weights = abs(E(Soil_up_net6)$weight),membership(cluster_walktrap(Soil_up_net6,weight=abs(E(Soil_up_net6)$weight),steps=4,modularity=TRUE))) #0.28
Soil_up_net6_distance<-mean_distance(Soil_up_net6,directed = FALSE) #1.83

## ---- Bulk_1wk_up_net6 ------

Bulk_1wk_up_net6
Bulk_1wk_up_net6_ratio<-length(which(E(Bulk_1wk_up_net6)$weight>0))/(length(E(Bulk_1wk_up_net6)$weight)) #53.13%
Bulk_1wk_up_net6_density<-edge_density(Bulk_1wk_up_net6,loops = F) # 40.41%
Bulk_1wk_up_net6_ave_degree<-sum(degree(Bulk_1wk_up_net6,mode="all"))/length(V(Bulk_1wk_up_net6)$name) #19.8
Bulk_1wk_up_net6_centr_degree<-centr_degree(Bulk_1wk_up_net6,mode="all",normalized = TRUE)$centralization #0.41
Bulk_1wk_up_net6_closeness<-sum(closeness(Bulk_1wk_up_net6,mode = "all",weights=rep(1,495),normalized = TRUE))/length(V(Bulk_1wk_up_net6)$name) #0.62
Bulk_1wk_up_net6_centr_closeness<-centr_clo(Bulk_1wk_up_net6,mode = "all",normalized = TRUE)$centralization #0.46
Bulk_1wk_up_net6_betweeness<-sum(betweenness(Bulk_1wk_up_net6,directed = FALSE,weight=rep(1,495)))/length(V(Bulk_1wk_up_net6)$name) #15.26
Bulk_1wk_up_net6_centr_betweenness<-centr_betw(Bulk_1wk_up_net6,directed = FALSE)$centralization #0.05
Bulk_1wk_up_net6_modularity<-modularity(Bulk_1wk_up_net6,weights = abs(E(Bulk_1wk_up_net6)$weight),membership(cluster_walktrap(Bulk_1wk_up_net6,weight=abs(E(Bulk_1wk_up_net6)$weight),steps=4,modularity=TRUE))) #0.14
Bulk_1wk_up_net6_distance<-mean_distance(Bulk_1wk_up_net6,directed = FALSE) #1.62

## ---- Bulk_3wk_up_net6 ------

Bulk_3wk_up_net6
Bulk_3wk_up_net6_ratio<-length(which(E(Bulk_3wk_up_net6)$weight>0))/(length(E(Bulk_3wk_up_net6)$weight)) #53.83%
Bulk_3wk_up_net6_density<-edge_density(Bulk_3wk_up_net6,loops = F) # 36.24%
Bulk_3wk_up_net6_ave_degree<-sum(degree(Bulk_3wk_up_net6,mode="all"))/length(V(Bulk_3wk_up_net6)$name) #17.76
Bulk_3wk_up_net6_centr_degree<-centr_degree(Bulk_3wk_up_net6,mode="all",normalized = TRUE)$centralization #0.23
Bulk_3wk_up_net6_closeness<-sum(closeness(Bulk_3wk_up_net6,mode = "all",weights=rep(1,444),normalized = TRUE))/length(V(Bulk_3wk_up_net6)$name) #0.60
Bulk_3wk_up_net6_centr_closeness<-centr_clo(Bulk_3wk_up_net6,mode = "all",normalized = TRUE)$centralization #0.22
Bulk_3wk_up_net6_betweeness<-sum(betweenness(Bulk_3wk_up_net6,directed = FALSE,weight=rep(1,444)))/length(V(Bulk_3wk_up_net6)$name) #16.38
Bulk_3wk_up_net6_centr_betweenness<-centr_betw(Bulk_3wk_up_net6,directed = FALSE)$centralization #0.02
Bulk_3wk_up_net6_modularity<-modularity(Bulk_3wk_up_net6,weights = abs(E(Bulk_3wk_up_net6)$weight),membership(cluster_walktrap(Bulk_3wk_up_net6,weight=abs(E(Bulk_3wk_up_net6)$weight),steps=4,modularity=TRUE))) #0.15
Bulk_3wk_up_net6_distance<-mean_distance(Bulk_3wk_up_net6,directed = FALSE) #1.67

## ---- Bulk_4wk_up_net6 ------

Bulk_4wk_up_net6
Bulk_4wk_up_net6_ratio<-length(which(E(Bulk_4wk_up_net6)$weight>0))/(length(E(Bulk_4wk_up_net6)$weight)) #52.83%
Bulk_4wk_up_net6_density<-edge_density(Bulk_4wk_up_net6,loops = F) # 38.94%
Bulk_4wk_up_net6_ave_degree<-sum(degree(Bulk_4wk_up_net6,mode="all"))/length(V(Bulk_4wk_up_net6)$name) #19.08
Bulk_4wk_up_net6_centr_degree<-centr_degree(Bulk_4wk_up_net6,mode="all",normalized = TRUE)$centralization #0.20
Bulk_4wk_up_net6_closeness<-sum(closeness(Bulk_4wk_up_net6,mode = "all",weights=rep(1,477),normalized = TRUE))/length(V(Bulk_4wk_up_net6)$name) #0.62
Bulk_4wk_up_net6_centr_closeness<-centr_clo(Bulk_4wk_up_net6,mode = "all",normalized = TRUE)$centralization #0.20
Bulk_4wk_up_net6_betweeness<-sum(betweenness(Bulk_4wk_up_net6,directed = FALSE,weight=rep(1,477)))/length(V(Bulk_4wk_up_net6)$name) #15.64
Bulk_4wk_up_net6_centr_betweenness<-centr_betw(Bulk_4wk_up_net6,directed = FALSE)$centralization #0.02
Bulk_4wk_up_net6_modularity<-modularity(Bulk_4wk_up_net6,weights = abs(E(Bulk_4wk_up_net6)$weight),membership(cluster_walktrap(Bulk_4wk_up_net6,weight=abs(E(Bulk_4wk_up_net6)$weight),steps=4,modularity=TRUE))) #0.17
Bulk_4wk_up_net6_distance<-mean_distance(Bulk_4wk_up_net6,directed = FALSE) #1.64

## ---- Rhi_1wk_up_net6 ------

Rhi_1wk_up_net6
Rhi_1wk_up_net6_ratio<-length(which(E(Rhi_1wk_up_net6)$weight>0))/(length(E(Rhi_1wk_up_net6)$weight)) #58.36%
Rhi_1wk_up_net6_density<-edge_density(Rhi_1wk_up_net6,loops = F) # 26.86%
Rhi_1wk_up_net6_ave_degree<-sum(degree(Rhi_1wk_up_net6,mode="all"))/length(V(Rhi_1wk_up_net6)$name) #13.16
Rhi_1wk_up_net6_centr_degree<-centr_degree(Rhi_1wk_up_net6,mode="all",normalized = TRUE)$centralization #0.28
Rhi_1wk_up_net6_closeness<-sum(closeness(Rhi_1wk_up_net6,mode = "all",weights=rep(1,329),normalized = TRUE))/length(V(Rhi_1wk_up_net6)$name) #0.54
Rhi_1wk_up_net6_centr_closeness<-centr_clo(Rhi_1wk_up_net6,mode = "all",normalized = TRUE)$centralization #0.28
Rhi_1wk_up_net6_betweeness<-sum(betweenness(Rhi_1wk_up_net6,directed = FALSE,weight=rep(1,329)))/length(V(Rhi_1wk_up_net6)$name) #21.2
Rhi_1wk_up_net6_centr_betweenness<-centr_betw(Rhi_1wk_up_net6,directed = FALSE)$centralization #0.06
Rhi_1wk_up_net6_modularity<-modularity(Rhi_1wk_up_net6,weights = abs(E(Rhi_1wk_up_net6)$weight),membership(cluster_walktrap(Rhi_1wk_up_net6,weight=abs(E(Rhi_1wk_up_net6)$weight),steps=4,modularity=TRUE))) #0.18
Rhi_1wk_up_net6_distance<-mean_distance(Rhi_1wk_up_net6,directed = FALSE) #1.87

## ---- Rhi_3wk_up_net6 ------

Rhi_3wk_up_net6
Rhi_3wk_up_net6_ratio<-length(which(E(Rhi_3wk_up_net6)$weight>0))/(length(E(Rhi_3wk_up_net6)$weight)) #61.13%
Rhi_3wk_up_net6_density<-edge_density(Rhi_3wk_up_net6,loops = F) # 41.80%
Rhi_3wk_up_net6_ave_degree<-sum(degree(Rhi_3wk_up_net6,mode="all"))/length(V(Rhi_3wk_up_net6)$name) #20.48
Rhi_3wk_up_net6_centr_degree<-centr_degree(Rhi_3wk_up_net6,mode="all",normalized = TRUE)$centralization #0.28
Rhi_3wk_up_net6_closeness<-sum(closeness(Rhi_3wk_up_net6,mode = "all",weights=rep(1,512),normalized = TRUE))/length(V(Rhi_3wk_up_net6)$name) #0.63
Rhi_3wk_up_net6_centr_closeness<-centr_clo(Rhi_3wk_up_net6,mode = "all",normalized = TRUE)$centralization #0.27
Rhi_3wk_up_net6_betweeness<-sum(betweenness(Rhi_3wk_up_net6,directed = FALSE,weight=rep(1,512)))/length(V(Rhi_3wk_up_net6)$name) #14.56
Rhi_3wk_up_net6_centr_betweenness<-centr_betw(Rhi_3wk_up_net6,directed = FALSE)$centralization #0.03
Rhi_3wk_up_net6_modularity<-modularity(Rhi_3wk_up_net6,weights = abs(E(Rhi_3wk_up_net6)$weight),membership(cluster_walktrap(Rhi_3wk_up_net6,weight=abs(E(Rhi_3wk_up_net6)$weight),steps=4,modularity=TRUE))) #0.19
Rhi_3wk_up_net6_distance<-mean_distance(Rhi_3wk_up_net6,directed = FALSE) #1.59

## ---- Rhi_4wk_up_net6 ------

Rhi_4wk_up_net6
Rhi_4wk_up_net6_ratio<-length(which(E(Rhi_4wk_up_net6)$weight>0))/(length(E(Rhi_4wk_up_net6)$weight)) #51.08%
Rhi_4wk_up_net6_density<-edge_density(Rhi_4wk_up_net6,loops = F) # 49.06%
Rhi_4wk_up_net6_ave_degree<-sum(degree(Rhi_4wk_up_net6,mode="all"))/length(V(Rhi_4wk_up_net6)$name) #24.04
Rhi_4wk_up_net6_centr_degree<-centr_degree(Rhi_4wk_up_net6,mode="all",normalized = TRUE)$centralization #0.33
Rhi_4wk_up_net6_closeness<-sum(closeness(Rhi_4wk_up_net6,mode = "all",weights=rep(1,601),normalized = TRUE))/length(V(Rhi_4wk_up_net6)$name) #0.67
Rhi_4wk_up_net6_centr_closeness<-centr_clo(Rhi_4wk_up_net6,mode = "all",normalized = TRUE)$centralization #0.36
Rhi_4wk_up_net6_betweeness<-sum(betweenness(Rhi_4wk_up_net6,directed = FALSE,weight=rep(1,601)))/length(V(Rhi_4wk_up_net6)$name) #12.48
Rhi_4wk_up_net6_centr_betweenness<-centr_betw(Rhi_4wk_up_net6,directed = FALSE)$centralization #0.03
Rhi_4wk_up_net6_modularity<-modularity(Rhi_4wk_up_net6,weights = abs(E(Rhi_4wk_up_net6)$weight),membership(cluster_walktrap(Rhi_4wk_up_net6,weight=abs(E(Rhi_4wk_up_net6)$weight),steps=4,modularity=TRUE))) #0.13
Rhi_4wk_up_net6_distance<-mean_distance(Rhi_4wk_up_net6,directed = FALSE) #1.51

## ---- Endo_1wk_up_net6 ------

Endo_1wk_up_net6
Endo_1wk_up_net6_ratio<-length(which(E(Endo_1wk_up_net6)$weight>0))/(length(E(Endo_1wk_up_net6)$weight)) #73.22%
Endo_1wk_up_net6_density<-edge_density(Endo_1wk_up_net6,loops = F) # 14.94%
Endo_1wk_up_net6_ave_degree<-sum(degree(Endo_1wk_up_net6,mode="all"))/length(V(Endo_1wk_up_net6)$name) #7.32
Endo_1wk_up_net6_centr_degree<-centr_degree(Endo_1wk_up_net6,mode="all",normalized = TRUE)$centralization #0.16
Endo_1wk_up_net6_closeness<-sum(closeness(Endo_1wk_up_net6,mode = "all",weights=rep(1,183),normalized = TRUE))/length(V(Endo_1wk_up_net6)$name) #0.33
Endo_1wk_up_net6_centr_closeness<-centr_clo(Endo_1wk_up_net6,mode = "all",normalized = TRUE)$centralization #0.24
Endo_1wk_up_net6_betweeness<-sum(betweenness(Endo_1wk_up_net6,directed = FALSE,weight=rep(1,183)))/length(V(Endo_1wk_up_net6)$name) #50.78
Endo_1wk_up_net6_centr_betweenness<-centr_betw(Endo_1wk_up_net6,directed = FALSE)$centralization #0.21
Endo_1wk_up_net6_modularity<-modularity(Endo_1wk_up_net6,weights = abs(E(Endo_1wk_up_net6)$weight),membership(cluster_walktrap(Endo_1wk_up_net6,weight=abs(E(Endo_1wk_up_net6)$weight),steps=4,modularity=TRUE))) #0.50
Endo_1wk_up_net6_distance<-mean_distance(Endo_1wk_up_net6,directed = FALSE) #3.07

## ---- Endo_3wk_up_net6 ------

Endo_3wk_up_net6
Endo_3wk_up_net6_ratio<-length(which(E(Endo_3wk_up_net6)$weight>0))/(length(E(Endo_3wk_up_net6)$weight)) #55.26%
Endo_3wk_up_net6_density<-edge_density(Endo_3wk_up_net6,loops = F) # 24.81%
Endo_3wk_up_net6_ave_degree<-sum(degree(Endo_3wk_up_net6,mode="all"))/length(V(Endo_3wk_up_net6)$name) #12.16
Endo_3wk_up_net6_centr_degree<-centr_degree(Endo_3wk_up_net6,mode="all",normalized = TRUE)$centralization #0.22
Endo_3wk_up_net6_closeness<-sum(closeness(Endo_3wk_up_net6,mode = "all",weights=rep(1,304),normalized = TRUE))/length(V(Endo_3wk_up_net6)$name) #0.53
Endo_3wk_up_net6_centr_closeness<-centr_clo(Endo_3wk_up_net6,mode = "all",normalized = TRUE)$centralization #0.23
Endo_3wk_up_net6_betweeness<-sum(betweenness(Endo_3wk_up_net6,directed = FALSE,weight=rep(1,304)))/length(V(Endo_3wk_up_net6)$name) #21.88
Endo_3wk_up_net6_centr_betweenness<-centr_betw(Endo_3wk_up_net6,directed = FALSE)$centralization #0.04
Endo_3wk_up_net6_modularity<-modularity(Endo_3wk_up_net6,weights = abs(E(Endo_3wk_up_net6)$weight),membership(cluster_walktrap(Endo_3wk_up_net6,weight=abs(E(Endo_3wk_up_net6)$weight),steps=4,modularity=TRUE))) #0.22
Endo_3wk_up_net6_distance<-mean_distance(Endo_3wk_up_net6,directed = FALSE) #1.89

## ---- Endo_4wk_up_net6 ------

Endo_4wk_up_net6
Endo_4wk_up_net6_ratio<-length(which(E(Endo_4wk_up_net6)$weight>0))/(length(E(Endo_4wk_up_net6)$weight)) #50.84%
Endo_4wk_up_net6_density<-edge_density(Endo_4wk_up_net6,loops = F) # 29.06%
Endo_4wk_up_net6_ave_degree<-sum(degree(Endo_4wk_up_net6,mode="all"))/length(V(Endo_4wk_up_net6)$name) #14.24
Endo_4wk_up_net6_centr_degree<-centr_degree(Endo_4wk_up_net6,mode="all",normalized = TRUE)$centralization #0.28
Endo_4wk_up_net6_closeness<-sum(closeness(Endo_4wk_up_net6,mode = "all",weights=rep(1,356),normalized = TRUE))/length(V(Endo_4wk_up_net6)$name) #0.57
Endo_4wk_up_net6_centr_closeness<-centr_clo(Endo_4wk_up_net6,mode = "all",normalized = TRUE)$centralization #0.26
Endo_4wk_up_net6_betweeness<-sum(betweenness(Endo_4wk_up_net6,directed = FALSE,weight=rep(1,356)))/length(V(Endo_4wk_up_net6)$name) #18.52
Endo_4wk_up_net6_centr_betweenness<-centr_betw(Endo_4wk_up_net6,directed = FALSE)$centralization #0.09
Endo_4wk_up_net6_modularity<-modularity(Endo_4wk_up_net6,weights = abs(E(Endo_4wk_up_net6)$weight),membership(cluster_walktrap(Endo_4wk_up_net6,weight=abs(E(Endo_4wk_up_net6)$weight),steps=4,modularity=TRUE))) #0.20
Endo_4wk_up_net6_distance<-mean_distance(Endo_4wk_up_net6,directed = FALSE) #1.76

## ---- Seed_up_net6 ------

Seed_up_net6
Seed_up_net6_ratio<-length(which(E(Seed_up_net6)$weight>0))/(length(E(Seed_up_net6)$weight)) #46.87%
Seed_up_net6_density<-edge_density(Seed_up_net6,loops = F) # 48.48%
Seed_up_net6_ave_degree<-sum(degree(Seed_up_net6,mode="all"))/length(V(Seed_up_net6)$name) #5.33
Seed_up_net6_centr_degree<-centr_degree(Seed_up_net6,mode="all",normalized = TRUE)$centralization #0.24
Seed_up_net6_closeness<-sum(closeness(Seed_up_net6,mode = "all",weights=rep(1,32),normalized = TRUE))/length(V(Seed_up_net6)$name) #0.65
Seed_up_net6_centr_closeness<-centr_clo(Seed_up_net6,mode = "all",normalized = TRUE)$centralization #0.31
Seed_up_net6_betweeness<-sum(betweenness(Seed_up_net6,directed = FALSE,weight=rep(1,32)))/length(V(Seed_up_net6)$name) #3.08
Seed_up_net6_centr_betweenness<-centr_betw(Seed_up_net6,directed = FALSE)$centralization #0.17
Seed_up_net6_modularity<-modularity(Seed_up_net6,weights = abs(E(Seed_up_net6)$weight),membership(cluster_walktrap(Seed_up_net6,weight=abs(E(Seed_up_net6)$weight),steps=4,modularity=TRUE))) #0.15
Seed_up_net6_distance<-mean_distance(Seed_up_net6,directed = FALSE) #1.56
```

## Calculate shared nodes 

```{r}
## before generate networks

Soil_up_names<-data.frame(OTU=labels(sort(taxa_sums(Soil),decreasing = TRUE)[1:500]),Soil_up="Soil_up")
Bulk_1wk_up_names<-data.frame(OTU=labels(sort(taxa_sums(Bulk_1wk),decreasing = TRUE)[1:500]),Bulk_1wk_up="Bulk_1wk_up")
Bulk_3wk_up_names<-data.frame(OTU=labels(sort(taxa_sums(Bulk_3wk),decreasing = TRUE)[1:500]),Bulk_3wk_up="Bulk_3wk_up")
Bulk_4wk_up_names<-data.frame(OTU=labels(sort(taxa_sums(Bulk_4wk),decreasing = TRUE)[1:500]),Bulk_4wk_up="Bulk_4wk_up")

Rhi_1wk_up_names<-data.frame(OTU=labels(sort(taxa_sums(Rhi_1wk),decreasing = TRUE)[1:500]),Rhi_1wk_up="Rhi_1wk_up")
Rhi_3wk_up_names<-data.frame(OTU=labels(sort(taxa_sums(Rhi_3wk),decreasing = TRUE)[1:500]),Rhi_3wk_up="Rhi_3wk_up")
Rhi_4wk_up_names<-data.frame(OTU=labels(sort(taxa_sums(Rhi_4wk),decreasing = TRUE)[1:500]),Rhi_4wk_up="Rhi_4wk_up")


Endo_1wk_up_names<-data.frame(OTU=labels(sort(taxa_sums(Endo_1wk),decreasing = TRUE)[1:500]),Endo_1wk_up="Endo_1wk_up")
Endo_3wk_up_names<-data.frame(OTU=labels(sort(taxa_sums(Endo_3wk),decreasing = TRUE)[1:500]),Endo_3wk_up="Endo_3wk_up")
Endo_4wk_up_names<-data.frame(OTU=labels(sort(taxa_sums(Endo_4wk),decreasing = TRUE)[1:500]),Endo_4wk_up="Endo_4wk_up")

Seed_up_names<-data.frame(OTU=labels((taxa_sums(SEED_top50))),Seed_up="Seed_up")


#--- all shared ------
all_names<-join_all(dfs = list(Soil_up_names,Bulk_1wk_up_names,Bulk_3wk_up_names,Bulk_4wk_up_names,Rhi_1wk_up_names,Rhi_3wk_up_names,Rhi_4wk_up_names,Endo_1wk_up_names,Endo_3wk_up_names,Endo_4wk_up_names,Seed_up_names),by = "OTU",type="full")
dim(all_names)
summary(complete.cases(all_names))[3]
all_names$OTU[complete.cases(all_names)]



# --- Bulk shared  ----

Bulk_names<-join_all(dfs = list(Bulk_1wk_up_names,Bulk_3wk_up_names,Bulk_4wk_up_names),by = "OTU",type="full")
dim(Bulk_names)
summary(complete.cases(Bulk_names))[3]
Bulk_shared_OTUID<-Bulk_names$OTU[complete.cases(Bulk_names)]
Seed_13021_up_tax<-data.frame(tax_table(Seed_13021_up)@.Data)
Seed_13021_up_tax$OTU<-rownames(Seed_13021_up_tax)
Bulk_shared_tax<-Seed_13021_up_tax[match(as.character(Bulk_shared_OTUID),Seed_13021_up_tax$OTU),]
dim(Bulk_shared_tax)
sort(summary(Bulk_shared_tax$Phylum),decreasing = TRUE)
shared_percentage_betweenwk134<-dim(Bulk_shared_tax)[1]/dim(Bulk_names)[1]

# --- Rhi shared  ----

Rhi_names<-join_all(dfs = list(Rhi_1wk_up_names,Rhi_3wk_up_names,Rhi_4wk_up_names),by = "OTU",type="full")
dim(Rhi_names)
summary(complete.cases(Rhi_names))[3]
Rhi_shared_OTUID<-Rhi_names$OTU[complete.cases(Rhi_names)]
Seed_13021_up_tax<-data.frame(tax_table(Seed_13021_up)@.Data)
Seed_13021_up_tax$OTU<-rownames(Seed_13021_up_tax)
Rhi_shared_tax<-Seed_13021_up_tax[match(as.character(Rhi_shared_OTUID),Seed_13021_up_tax$OTU),]
dim(Rhi_shared_tax)
sort(summary(Rhi_shared_tax$Phylum),decreasing = TRUE)
shared_percentage_betweenwk134<-dim(Rhi_shared_tax)[1]/dim(Rhi_names)[1]

# --- Endo shared  ----

Endo_names<-join_all(dfs = list(Endo_1wk_up_names,Endo_3wk_up_names,Endo_4wk_up_names),by = "OTU",type="full")
dim(Endo_names)
summary(complete.cases(Endo_names))[3]
Endo_shared_OTUID<-Endo_names$OTU[complete.cases(Endo_names)]
Seed_13021_up_tax<-data.frame(tax_table(Seed_13021_up)@.Data)
Seed_13021_up_tax$OTU<-rownames(Seed_13021_up_tax)
Endo_shared_tax<-Seed_13021_up_tax[match(as.character(Endo_shared_OTUID),Seed_13021_up_tax$OTU),]
dim(Endo_shared_tax)
sort(summary(Endo_shared_tax$Phylum),decreasing = TRUE)
shared_percentage_betweenwk134<-dim(Endo_shared_tax)[1]/dim(Endo_names)[1]

# --- Net4 shared OTUs ----

Soil_up_net4_names<-data.frame(OTU=vertex_attr(Soil_up_net4)$name,Soil="Soil_up_net4")

Bulk_1wk_net4_names<-data.frame(OTU=vertex_attr(Bulk_1wk_up_net4)$name,Bulk_1wk_up_net4="Bulk_1wk_up_net4")
Bulk_3wk_net4_names<-data.frame(OTU=vertex_attr(Bulk_3wk_up_net4)$name,Bulk_3wk_up_net4="Bulk_3wk_up_net4")
Bulk_4wk_net4_names<-data.frame(OTU=vertex_attr(Bulk_4wk_up_net4)$name,Bulk_4wk_up_net4="Bulk_4wk_up_net4")

Rhi_1wk_net4_names<-data.frame(OTU=vertex_attr(Rhi_1wk_up_net4)$name,Rhi_1wk_up_net4="Rhi_1wk_up_net4")
Rhi_3wk_net4_names<-data.frame(OTU=vertex_attr(Rhi_3wk_up_net4)$name,Rhi_3wk_up_net4="Rhi_3wk_up_net4")
Rhi_4wk_net4_names<-data.frame(OTU=vertex_attr(Rhi_4wk_up_net4)$name,Rhi_4wk_up_net4="Rhi_4wk_up_net4")

Endo_1wk_net4_names<-data.frame(OTU=vertex_attr(Endo_1wk_up_net4)$name,Endo_1wk_up_net4="Endo_1wk_up_net4")
Endo_3wk_net4_names<-data.frame(OTU=vertex_attr(Endo_3wk_up_net4)$name,Endo_3wk_up_net4="Endo_3wk_up_net4")
Endo_4wk_net4_names<-data.frame(OTU=vertex_attr(Endo_4wk_up_net4)$name,Endo_4wk_up_net4="Endo_4wk_up_net4")

Seed_up_net4_names<-data.frame(OTU=vertex_attr(Seed_up_net4)$name,Seed_up_net4="Seed_up_net4")

all_net4_names<-join_all(dfs = list(Soil_up_net4_names,Bulk_1wk_net4_names,Bulk_3wk_net4_names,Bulk_4wk_net4_names,Rhi_1wk_net4_names,Rhi_3wk_net4_names,Rhi_4wk_net4_names,Endo_1wk_net4_names,Endo_3wk_net4_names,Endo_4wk_net4_names,Seed_up_net4_names),by = "OTU",type="full")
dim(all_net4_names)
summary(complete.cases(all_net4_names))[3]
all_net4_names$OTU[complete.cases(all_net4_names)]

## --- all net4 shared ------

all_net4_shared_OTUID<-all_net4_names$OTU[complete.cases(all_net4_names)]
Seed_13021_up_tax<-data.frame(tax_table(Seed_13021_up)@.Data)
Seed_13021_up_tax$OTU<-rownames(Seed_13021_up_tax)
all_net4_shared_tax<-Seed_13021_up_tax[match(as.character(all_net4_shared_OTUID),Seed_13021_up_tax$OTU),]
dim(all_net4_shared_tax)
sort(summary(all_net4_shared_tax$Phylum),decreasing = TRUE)

# --- Bulk net4_shared  ----

Bulk_net4_names<-join_all(dfs = list(Bulk_1wk_net4_names,Bulk_3wk_net4_names,Bulk_4wk_net4_names),by = "OTU",type="full")
dim(Bulk_net4_names)
summary(complete.cases(Bulk_net4_names))[3]
Bulk_net4_shared_OTUID<-Bulk_net4_names$OTU[complete.cases(Bulk_net4_names)]
Seed_13021_up_tax<-data.frame(tax_table(Seed_13021_up)@.Data)
Seed_13021_up_tax$OTU<-rownames(Seed_13021_up_tax)
Bulk_net4_shared_tax<-Seed_13021_up_tax[match(as.character(Bulk_net4_shared_OTUID),Seed_13021_up_tax$OTU),]
dim(Bulk_net4_shared_tax)
sort(summary(Bulk_net4_shared_tax$Phylum),decreasing = TRUE)
shared_percentage_betweenwk134<-dim(Bulk_net4_shared_tax)[1]/dim(Bulk_net4_names)[1]

# --- Rhi net4_shared  ----

Rhi_net4_names<-join_all(dfs = list(Rhi_1wk_net4_names,Rhi_3wk_net4_names,Rhi_4wk_net4_names),by = "OTU",type="full")
dim(Rhi_net4_names)
summary(complete.cases(Rhi_net4_names))[3]
Rhi_net4_shared_OTUID<-Rhi_net4_names$OTU[complete.cases(Rhi_net4_names)]
Seed_13021_up_tax<-data.frame(tax_table(Seed_13021_up)@.Data)
Seed_13021_up_tax$OTU<-rownames(Seed_13021_up_tax)
Rhi_net4_shared_tax<-Seed_13021_up_tax[match(as.character(Rhi_net4_shared_OTUID),Seed_13021_up_tax$OTU),]
dim(Rhi_net4_shared_tax)
sort(summary(Rhi_net4_shared_tax$Phylum),decreasing = TRUE)
shared_percentage_betweenwk134<-dim(Rhi_net4_shared_tax)[1]/dim(Rhi_net4_names)[1]

# --- Endo net4_shared  ----

Endo_net4_names<-join_all(dfs = list(Endo_1wk_net4_names,Endo_3wk_net4_names,Endo_4wk_net4_names),by = "OTU",type="full")
dim(Endo_net4_names)
summary(complete.cases(Endo_net4_names))[3]
Endo_net4_shared_OTUID<-Endo_net4_names$OTU[complete.cases(Endo_net4_names)]
Seed_13021_up_tax<-data.frame(tax_table(Seed_13021_up)@.Data)
Seed_13021_up_tax$OTU<-rownames(Seed_13021_up_tax)
Endo_net4_shared_tax<-Seed_13021_up_tax[match(as.character(Endo_net4_shared_OTUID),Seed_13021_up_tax$OTU),]
dim(Endo_net4_shared_tax)
sort(summary(Endo_net4_shared_tax$Phylum),decreasing = TRUE)
shared_percentage_betweenwk134<-dim(Endo_net4_shared_tax)[1]/dim(Endo_net4_names)[1]

# --- Net6 shared OTUs ---

Soil_up_net6_names<-data.frame(OTU=vertex_attr(Soil_up_net6)$name,Soil="Soil_up_net6")

Bulk_1wk_net6_names<-data.frame(OTU=vertex_attr(Bulk_1wk_up_net6)$name,Bulk_1wk_up_net6="Bulk_1wk_up_net6")
Bulk_3wk_net6_names<-data.frame(OTU=vertex_attr(Bulk_3wk_up_net6)$name,Bulk_3wk_up_net6="Bulk_3wk_up_net6")
Bulk_4wk_net6_names<-data.frame(OTU=vertex_attr(Bulk_4wk_up_net6)$name,Bulk_4wk_up_net6="Bulk_4wk_up_net6")

Rhi_1wk_net6_names<-data.frame(OTU=vertex_attr(Rhi_1wk_up_net6)$name,Rhi_1wk_up_net6="Rhi_1wk_up_net6")
Rhi_3wk_net6_names<-data.frame(OTU=vertex_attr(Rhi_3wk_up_net6)$name,Rhi_3wk_up_net6="Rhi_3wk_up_net6")
Rhi_4wk_net6_names<-data.frame(OTU=vertex_attr(Rhi_4wk_up_net6)$name,Rhi_4wk_up_net6="Rhi_4wk_up_net6")

Endo_1wk_net6_names<-data.frame(OTU=vertex_attr(Endo_1wk_up_net6)$name,Endo_1wk_up_net6="Endo_1wk_up_net6")
Endo_3wk_net6_names<-data.frame(OTU=vertex_attr(Endo_3wk_up_net6)$name,Endo_3wk_up_net6="Endo_3wk_up_net6")
Endo_4wk_net6_names<-data.frame(OTU=vertex_attr(Endo_4wk_up_net6)$name,Endo_4wk_up_net6="Endo_4wk_up_net6")

Seed_up_net6_names<-data.frame(OTU=vertex_attr(Seed_up_net6)$name,Seed_up_net6="Seed_up_net6")

all_net6_names<-join_all(dfs = list(Soil_up_net6_names,Bulk_1wk_net6_names,Bulk_3wk_net6_names,Bulk_4wk_net6_names,Rhi_1wk_net6_names,Rhi_3wk_net6_names,Rhi_4wk_net6_names,Endo_1wk_net6_names,Endo_3wk_net6_names,Endo_4wk_net6_names),by = "OTU",type = "full")
summary(complete.cases(all_net6_names))[3]

## --- all net6 shared ------

all_net6_shared_OTUID<-all_net6_names$OTU[complete.cases(all_net6_names)]
Seed_13021_up_tax<-data.frame(tax_table(Seed_13021_up)@.Data)
Seed_13021_up_tax$OTU<-rownames(Seed_13021_up_tax)
all_net6_shared_tax<-Seed_13021_up_tax[match(as.character(all_net6_shared_OTUID),Seed_13021_up_tax$OTU),]
dim(all_net6_shared_tax)
sort(summary(all_net6_shared_tax$Phylum),decreasing = TRUE)

# --- Bulk net6_shared  ----

Bulk_net6_names<-join_all(dfs = list(Bulk_1wk_net6_names,Bulk_3wk_net6_names,Bulk_4wk_net6_names),by = "OTU",type="full")
dim(Bulk_net6_names)
summary(complete.cases(Bulk_net6_names))[3]
Bulk_net6_shared_OTUID<-Bulk_net6_names$OTU[complete.cases(Bulk_net6_names)]
Seed_13021_up_tax<-data.frame(tax_table(Seed_13021_up)@.Data)
Seed_13021_up_tax$OTU<-rownames(Seed_13021_up_tax)
Bulk_net6_shared_tax<-Seed_13021_up_tax[match(as.character(Bulk_net6_shared_OTUID),Seed_13021_up_tax$OTU),]
dim(Bulk_net6_shared_tax)
sort(summary(Bulk_net6_shared_tax$Phylum),decreasing = TRUE)


# --- Rhi net6_shared  ----

Rhi_net6_names<-join_all(dfs = list(Rhi_1wk_net6_names,Rhi_3wk_net6_names,Rhi_4wk_net6_names),by = "OTU",type="full")
dim(Rhi_net6_names)
summary(complete.cases(Rhi_net6_names))[3]
Rhi_net6_shared_OTUID<-Rhi_net6_names$OTU[complete.cases(Rhi_net6_names)]
Seed_13021_up_tax<-data.frame(tax_table(Seed_13021_up)@.Data)
Seed_13021_up_tax$OTU<-rownames(Seed_13021_up_tax)
Rhi_net6_shared_tax<-Seed_13021_up_tax[match(as.character(Rhi_net6_shared_OTUID),Seed_13021_up_tax$OTU),]
dim(Rhi_net6_shared_tax)
sort(summary(Rhi_net6_shared_tax$Phylum),decreasing = TRUE)

# --- Endo net6_shared  ----

Endo_net6_names<-join_all(dfs = list(Endo_1wk_net6_names,Endo_3wk_net6_names,Endo_4wk_net6_names),by = "OTU",type="full")
dim(Endo_net6_names)
summary(complete.cases(Endo_net6_names))[3]
Endo_net6_shared_OTUID<-Endo_net6_names$OTU[complete.cases(Endo_net6_names)]
Seed_13021_up_tax<-data.frame(tax_table(Seed_13021_up)@.Data)
Seed_13021_up_tax$OTU<-rownames(Seed_13021_up_tax)
Endo_net6_shared_tax<-Seed_13021_up_tax[match(as.character(Endo_net6_shared_OTUID),Seed_13021_up_tax$OTU),]
dim(Endo_net6_shared_tax)
sort(summary(Endo_net6_shared_tax$Phylum),decreasing = TRUE)

```
