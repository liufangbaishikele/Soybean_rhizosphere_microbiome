### SparCC for co-orccurence network

* SparCC input files preparasion

```
## set up working environment and load libraries

```{r}
setwd("/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/SparCC")
library(phyloseq)
library(igraph)
library(plyr)
```

## Read into a whole .shared and cons.taxonomy file to define label factor 
## This step is used to generate labels for each genus, which could be used to merge with four groups of vertex dataframe based on Genus name.

```{r}

# Read in the original strigolactone output OTU table and taxonomy table

Seed_22017<-import_mothur(mothur_shared_file = "/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/rarefaction_22017/Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.subsample.0.03.pick.shared",mothur_constaxonomy_file = "/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/rarefaction_22017/Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.subsample.0.03.cons.taxonomy")
Seed_22017 # 28114 OTUs across 159 samples
colnames(tax_table(Seed_22017))<-c('Kingdom','Phylum','Class','Order','Family','Genus')

# Read in the meta data
Seed_meta<-read.csv('/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/Seed_meta.csv',row.names = 1,header = TRUE)
head(Seed_meta)


Seed_meta_phyloseq<-sample_data(Seed_meta)

#>>>>> Seed_phyloseq_up >>>>> with sample data integrated

Seed_22017_up<-merge_phyloseq(Seed_22017,Seed_meta_phyloseq)
Seed_22017_up
sample_names(Seed_22017_up)[1:5]
Seed_22017_up
identical(sample_names(Seed_22017),sample_names(Seed_22017_up))

# creat label vector

meta<-data.frame(label=1:length(taxa_names(Seed_22017_up)),OtuID=taxa_names(Seed_22017_up),Phylum=tax_table(Seed_22017_up)[,2])
head(meta)
```

#  Subset OTUs to creat input for SparCC analysis
---------------------------------------------------

```{r}
#Subset the samples based on treatment

#---Soil---

Soil<-subset_samples(Seed_22017_up,Treat=="Soil")
Soil
#Extract top 1000 OTUs across 8 samples
Soil_up<-prune_taxa(labels(sort(taxa_sums(Soil),decreasing = TRUE)[1:1000]),Soil)
Soil_up #1000 taxa and 15 samples
min(taxa_sums(Soil_up)) # 50
Soil_up_OTU_count<-data.frame(otu_table(Soil_up))
#write.csv(Soil_up_OTU_count,file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/SparCC/Soil_up_OTU_count.csv")
table(Soil_up_OTU_count==0)[2]/sum(table(Soil_up_OTU_count==0)) #0.05493333

#---Bulk_1wk---

Bulk_1wk<-subset_samples(Seed_22017_up,Treat=="Bulk_1wk")
Bulk_1wk
#Extract top 1000 OTUs across 8 samples
Bulk_1wk_up<-prune_taxa(labels(sort(taxa_sums(Bulk_1wk),decreasing = TRUE)[1:1000]),Bulk_1wk)
Bulk_1wk_up #1000 taxa and 15 samples
min(taxa_sums(Bulk_1wk_up)) # 53
Bulk_1wk_up_OTU_count<-data.frame(otu_table(Bulk_1wk_up))
#write.csv(Bulk_1wk_up_OTU_count,file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/SparCC/Bulk_1wk_up_OTU_count.csv")
table(Bulk_1wk_up_OTU_count==0)[2]/sum(table(Bulk_1wk_up_OTU_count==0)) #0.05546667

#---Bulk_3wk---

Bulk_3wk<-subset_samples(Seed_22017_up,Treat=="Bulk_3wk")
Bulk_3wk
#Extract top 1000 OTUs across 8 samples
Bulk_3wk_up<-prune_taxa(labels(sort(taxa_sums(Bulk_3wk),decreasing = TRUE)[1:1000]),Bulk_3wk)
Bulk_3wk_up #1000 taxa and 15 samples
min(taxa_sums(Bulk_3wk_up)) # 51
Bulk_3wk_up_OTU_count<-data.frame(otu_table(Bulk_3wk_up))
#write.csv(Bulk_3wk_up_OTU_count,file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/SparCC/Bulk_3wk_up_OTU_count.csv")
table(Bulk_3wk_up_OTU_count==0)[2]/sum(table(Bulk_3wk_up_OTU_count==0)) #0.06713333

#---Bulk_4wk---

Bulk_4wk<-subset_samples(Seed_22017_up,Treat=="Bulk_4wk")
Bulk_4wk
#Extract top 1000 OTUs across 8 samples
Bulk_4wk_up<-prune_taxa(labels(sort(taxa_sums(Bulk_4wk),decreasing = TRUE)[1:1000]),Bulk_4wk)
Bulk_4wk_up #1000 taxa and 15 samples
min(taxa_sums(Bulk_4wk_up)) # 49
Bulk_4wk_up_OTU_count<-data.frame(otu_table(Bulk_4wk_up))
#write.csv(Bulk_4wk_up_OTU_count,file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/SparCC/Bulk_4wk_up_OTU_count.csv")
table(Bulk_4wk_up_OTU_count==0)[2]/sum(table(Bulk_4wk_up_OTU_count==0)) #0.06713333

#---Rhi_1wk---

Rhi_1wk<-subset_samples(Seed_22017_up,Treat=="Rhi_1wk")
Rhi_1wk
#Extract top 1000 OTUs across 8 samples
Rhi_1wk_up<-prune_taxa(labels(sort(taxa_sums(Rhi_1wk),decreasing = TRUE)[1:1000]),Rhi_1wk)
Rhi_1wk_up #1000 taxa and 15 samples
min(taxa_sums(Rhi_1wk_up)) # 22
Rhi_1wk_up_OTU_count<-data.frame(otu_table(Rhi_1wk_up))
#write.csv(Rhi_1wk_up_OTU_count,file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/SparCC/Rhi_1wk_up_OTU_count.csv")
table(Rhi_1wk_up_OTU_count==0)[2]/sum(table(Rhi_1wk_up_OTU_count==0)) #0.06713333

# --- Rhi_3wk ---

Rhi_3wk<-subset_samples(Seed_22017_up,Treat=="Rhi_3wk")
Rhi_3wk
#Extract top 1000 OTUs across 8 samples
Rhi_3wk_up<-prune_taxa(labels(sort(taxa_sums(Rhi_3wk),decreasing = TRUE)[1:1000]),Rhi_3wk)
Rhi_3wk_up #1000 taxa and 15 samples
min(taxa_sums(Rhi_3wk_up)) # 49
Rhi_3wk_up_OTU_count<-data.frame(otu_table(Rhi_3wk_up))
#write.csv(Rhi_3wk_up_OTU_count,file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/SparCC/Rhi_3wk_up_OTU_count.csv")
table(Rhi_3wk_up_OTU_count==0)[2]/sum(table(Rhi_3wk_up_OTU_count==0)) #0.06826667 

# --- Rhi_4wk ---

Rhi_4wk<-subset_samples(Seed_22017_up,Treat=="Rhi_4wk")
Rhi_4wk
#Extract top 1000 OTUs across 15 samples
Rhi_4wk_up<-prune_taxa(labels(sort(taxa_sums(Rhi_4wk),decreasing = TRUE)[1:1000]),Rhi_4wk)
Rhi_4wk_up #1000 taxa and 15 samples
min(taxa_sums(Rhi_4wk_up)) # 43
Rhi_4wk_up_OTU_count<-data.frame(otu_table(Rhi_4wk_up))
#write.csv(Rhi_4wk_up_OTU_count,file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/SparCC/Rhi_4wk_up_OTU_count.csv")
table(Rhi_4wk_up_OTU_count==0)[2]/sum(table(Rhi_4wk_up_OTU_count==0)) #0.06826667

# --- Endo_1wk ---

Endo_1wk<-subset_samples(Seed_22017_up,Treat=="Endo_1wk")
Endo_1wk
#Extract top 1000 OTUs across 15 samples
Endo_1wk_up<-prune_taxa(labels(sort(taxa_sums(Endo_1wk),decreasing = TRUE)[1:1000]),Endo_1wk)
Endo_1wk_up #1000 taxa and 15 samples
min(taxa_sums(Endo_1wk_up)) # 22
Endo_1wk_up_OTU_count<-data.frame(otu_table(Endo_1wk_up))
#write.csv(Endo_1wk_up_OTU_count,file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/SparCC/Endo_1wk_up_OTU_count.csv")
table(Endo_1wk_up_OTU_count==0)[2]/sum(table(Endo_1wk_up_OTU_count==0)) #0.2120667

# --- Endo_3wk ---

Endo_3wk<-subset_samples(Seed_22017_up,Treat=="Endo_3wk")
Endo_3wk
#Extract top 1000 OTUs across 15 samples
Endo_3wk_up<-prune_taxa(labels(sort(taxa_sums(Endo_3wk),decreasing = TRUE)[1:1000]),Endo_3wk)
Endo_3wk_up #1000 taxa and 15 samples
min(taxa_sums(Endo_3wk_up)) # 11
Endo_3wk_up_OTU_count<-data.frame(otu_table(Endo_3wk_up))
#write.csv(Endo_3wk_up_OTU_count,file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/SparCC/Endo_3wk_up_OTU_count.csv")
table(Endo_3wk_up_OTU_count==0)[2]/sum(table(Endo_3wk_up_OTU_count==0)) #0.3006

# --- Endo_4wk ---

Endo_4wk<-subset_samples(Seed_22017_up,Treat=="Endo_4wk")
Endo_4wk
#Extract top 1000 OTUs across 15 samples
Endo_4wk_up<-prune_taxa(labels(sort(taxa_sums(Endo_4wk),decreasing = TRUE)[1:1000]),Endo_4wk)
Endo_4wk_up #1000 taxa and 15 samples
min(taxa_sums(Endo_4wk_up)) # 12
Endo_4wk_up_OTU_count<-data.frame(otu_table(Endo_4wk_up))
#write.csv(Endo_4wk_up_OTU_count,file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/SparCC/Endo_4wk_up_OTU_count.csv")
table(Endo_4wk_up_OTU_count==0)[2]/sum(table(Endo_4wk_up_OTU_count==0)) #0.3006

#--- Seed ---
Seed<-subset_samples(Seed_22017_up,Treat=="Seed")
Seed
#Extract top 1000 OTUs across 15 samples
Seed_up<-prune_taxa(labels(sort(taxa_sums(Seed),decreasing = TRUE)[1:258]),Seed)
Seed_up #1000 taxa and 15 samples
min(taxa_sums(Seed_up)) # 11
Seed_up_OTU_count<-data.frame(otu_table(Seed_up))
#write.csv(Seed_up_OTU_count,file="/Users/fangliu/Documents/2018_fugicide_project/Seed_treatment/16S/Seed_16S_Community_analysis/SparCC/Seed_up_OTU_count.csv")
table(Seed_up_OTU_count==0)[2]/sum(table(Seed_up_OTU_count==0)) #0.7012274
```
```
