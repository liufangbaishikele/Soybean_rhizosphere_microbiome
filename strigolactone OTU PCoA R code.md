# Strigolactone project pre-cluster with diff=5 #

OTU table and taxonomy data were imported directly from mothur output
- shared file is the OTU table with label=0.03 means OTU are called at 97% of similarity
- taxonomy file include the OTU taxonomy information. In PCoA plot case, we did not use the taxonomy information
- csv file is the meta data that created using sampleID from otu_table and corresponding treatment information

```
source('http://bioconductor.org/biocLite.R')
biocLite('metagenomeSeq')
library(metagenomeSeq)
library(phyloseq)
library(ggplot2)
library(plyr)
library(vegan)

#----Data load from local computer using phyloseq----#
setwd("G:\\UT\\I Love my project\\1-soybean_16S_data_and_analysis\\aa_R file\\strigolactone project")

strigolactone_phyloseq_origin<-import_mothur(mothur_shared_file = "strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared", mothur_constaxonomy_file ="strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy")
SL_OTU<-otu_table(strigolactone_phyloseq_origin)[,1:48]
SL_taxa<-tax_table(strigolactone_phyloseq_origin)
colnames(SL_taxa)<-c("Kingdom","Phylum","Class","Order","Family","Genus")
#sample_ID<-colnames(SL_OTU)
#write.table(sample_ID,file="temp.csv",sep = ",") # Then I edit the meta file manually and saved
SL_meta_csv<-read.csv("SL_meta.csv",h=TRUE,row.names = 1)
SL_meta<-sample_data(SL_meta_csv)

#---Strigolactone phyloseq experiment_level class-------

#1) 48 samples

#a) original datasets
SL_phyloseq<-phyloseq(SL_OTU,SL_taxa,SL_meta)

#b) original datasets but relative abundance
r_SL_phyloseq<-transform_sample_counts(SL_phyloseq,function(x) x/sum(x))

#c) filtered datasets 
filter_SL_phyloseq<-filter_taxa(SL_phyloseq,function(x) sum(x > 1) > (0.2*length(x)),TRUE)# ntaxa(filter_SL_phyloseqed): 801

#d) filtered & relative abundance
r_filter_SL_phyloseq<-transform_sample_counts(filter_SL_phyloseq,function(x) x/sum(x))

#a)PCoA plot using weighted_Bray-Curtis

#SL_phyloseq_PCoA<-ordinate(SL_phyloseq,method="PCoA",distance="bray")
SL_phyloseq_bray<-vegdist(t(otu_table(SL_phyloseq)),method="bray",binary=FALSE,autotransform=FALSE)
SL_phyloseq_PCoA<-ordinate(SL_phyloseq,method="PCoA",SL_phyloseq_bray)
plot_ordination(SL_phyloseq,SL_phyloseq_PCoA,type="samples",color="Treatment")+
  geom_point(size=3)+facet_wrap(~Treatment)+ggtitle("SL_phyloseq weighted Bray PCoA plot")

plot_ordination(SL_phyloseq,SL_phyloseq_PCoA,type="samples",color="Treatment")+
  geom_point(size=4)+ggtitle("SL_phyloseq weighted Bray PCoA plot")

#a-jaccard) PCoA plot with jaccard distance

SL_phyloseq_jaccard<-vegdist(t(otu_table(SL_phyloseq)),method="jaccard",binary=FALSE,autotransform=FALSE)
SL_phyloseq_PCoA<-ordinate(SL_phyloseq,method="PCoA",SL_phyloseq_jaccard)
plot_ordination(SL_phyloseq,SL_phyloseq_PCoA,type="samples",color="Treatment")+
  geom_point(size=3)+facet_wrap(~Treatment)+ggtitle("SL_phyloseq weighted jaccard PCoA plot")

plot_ordination(SL_phyloseq,SL_phyloseq_PCoA,type="samples",color="Treatment")+
  geom_point(size=4)+ggtitle("SL_phyloseq weighted jaccard PCoA plot")

#a-euclidean) PCoA plot with jaccard distance

SL_phyloseq_euclidean<-vegdist(t(otu_table(SL_phyloseq)),method="euclidean",binary=FALSE,autotransform=FALSE)
SL_phyloseq_PCoA<-ordinate(SL_phyloseq,method="PCoA",SL_phyloseq_euclidean)
plot_ordination(SL_phyloseq,SL_phyloseq_PCoA,type="samples",color="Treatment")+
  geom_point(size=3)+facet_wrap(~Treatment)+ggtitle("SL_phyloseq weighted euclidean PCoA plot")

plot_ordination(SL_phyloseq,SL_phyloseq_PCoA,type="samples",color="Treatment")+
  geom_point(size=4)+ggtitle("SL_phyloseq weighted euclidean PCoA plot")


#b)PCoA plot using weighted_Bray-Curtis

r_SL_phyloseq_bray<-vegdist(t(otu_table(r_SL_phyloseq)),method="bray",binary=FALSE,autotransform=FALSE)
r_SL_phyloseq_PCoA<-ordinate(r_SL_phyloseq,method="PCoA",r_SL_phyloseq_bray)
plot_ordination(r_SL_phyloseq,r_SL_phyloseq_PCoA,type="samples",color="Treatment")+
  geom_point(size=3)+facet_wrap(~Treatment)+ggtitle("r_SL_phyloseq weighted Bray PCoA plot")

plot_ordination(r_SL_phyloseq,r_SL_phyloseq_PCoA,type="samples",color="Treatment")+
  geom_point(size=4)+ggtitle("r_SL_phyloseq weighted Bray PCoA plot")

#c)PCoA plot using weighted_Bray-Curtis

filter_SL_phyloseq_bray<-vegdist(t(otu_table(filter_SL_phyloseq)),method="bray",binary=FALSE,autotransform=FALSE)
filter_SL_phyloseq_PCoA<-ordinate(filter_SL_phyloseq,method="PCoA",filter_SL_phyloseq_bray)
plot_ordination(filter_SL_phyloseq,filter_SL_phyloseq_PCoA,type="samples",color="Treatment")+
  geom_point(size=3)+facet_wrap(~Treatment)+ggtitle("filter_SL_phyloseq weighted Bray PCoA plot")

plot_ordination(filter_SL_phyloseq,filter_SL_phyloseq_PCoA,type="samples",color="Treatment")+
  geom_point(size=4)+ggtitle("filter_SL_phyloseq weighted Bray PCoA plot")


#d)PCoA plot using weighted_Bray-Curtis

r_filter_SL_phyloseq_bray<-vegdist(t(otu_table(r_filter_SL_phyloseq)),method="bray",binary=FALSE,autotransform=FALSE)
r_filter_SL_phyloseq_PCoA<-ordinate(r_filter_SL_phyloseq,method="PCoA",r_filter_SL_phyloseq_bray)
plot_ordination(r_filter_SL_phyloseq,r_filter_SL_phyloseq_PCoA,type="samples",color="Treatment")+
  geom_point(size=3)+facet_wrap(~Treatment)+ggtitle("r_filter_SL_phyloseq weighted Bray PCoA plot")

plot_ordination(r_filter_SL_phyloseq,r_filter_SL_phyloseq_PCoA,type="samples",color="Treatment")+
  geom_point(size=4)+ggtitle("r_filter_SL_phyloseq weighted Bray PCoA plot")


#a) PCoA plot using unweighted_Bray-Curtis
SL_phyloseq_bray_un<-vegdist(t(otu_table(SL_phyloseq)),method="bray",binary=TRUE,autotransform=FALSE,na.rm = FALSE)
# error pop out, haven't figured out
#b)PCoA plot using unweighted_Bray-Curtis
#c)PCoA plot using unweighted_Bray-Curtis
#d)PCoA plot using unweighted_Bray-Curtis

#2) prune dataset to keep 6 replications of each transgenic construct treatment

prune_ID<-rownames(SL_meta)[c(1:6,11:16,21:26,31:36,41:48)]

#e) Keep datasets with samples that barcoded right

prune_SL_phyloseq<-prune_samples(prune_ID,SL_phyloseq) #sample.names(prune_SL_phyloseq)

#f) Keep datasets with samples that barcoded right using relative abundance
r_prune_SL_phyloseq<-transform_sample_counts(prune_SL_phyloseq,function(x) x/sum(x))

#g) filter off OTUs that showed up in less than 20% of the samples
filter_prune_SL_phyloseq<-filter_taxa(prune_SL_phyloseq,function(x) sum(x > 1) > (0.2*length(x)),TRUE)# ntaxa(filter_prune_SL_phyloseqed):935 

#h) filter off OTUs that showed up in less than 20% of the samples, OTU counts transformed to relative abundance
r_filter_prune_SL_phyloseq<-transform_sample_counts(filter_prune_SL_phyloseq,function(x) x/sum(x))


#e) Weighted Bray curtis PCoA plot
prune_SL_phyloseq_bray<-vegdist(t(otu_table(prune_SL_phyloseq)),method="bray",binary=FALSE,autotransform=FALSE)
prune_SL_phyloseq_PCoA<-ordinate(prune_SL_phyloseq,method="PCoA",prune_SL_phyloseq_bray)
plot_ordination(prune_SL_phyloseq,prune_SL_phyloseq_PCoA,type="samples",color="Treatment")+
  geom_point(size=3)+facet_wrap(~Treatment)+ggtitle("prune_SL_phyloseq weighted Bray PCoA plot")

plot_ordination(prune_SL_phyloseq,prune_SL_phyloseq_PCoA,type="samples",color="Treatment")+
  geom_point(size=4)+ggtitle("prune_SL_phyloseq weighted Bray PCoA plot")


#f)Weighted Bray curtis PCoA plot

r_prune_SL_phyloseq_bray<-vegdist(t(otu_table(r_prune_SL_phyloseq)),method="bray",binary=FALSE,autotransform=FALSE)
r_prune_SL_phyloseq_PCoA<-ordinate(r_prune_SL_phyloseq,method="PCoA",r_prune_SL_phyloseq_bray)
plot_ordination(r_prune_SL_phyloseq,r_prune_SL_phyloseq_PCoA,type="samples",color="Treatment")+
  geom_point(size=3)+facet_wrap(~Treatment)+ggtitle("r_prune_SL_phyloseq weighted Bray PCoA plot")

plot_ordination(r_prune_SL_phyloseq,r_prune_SL_phyloseq_PCoA,type="samples",color="Treatment")+
  geom_point(size=4)+ggtitle("r_prune_SL_phyloseq weighted Bray PCoA plot")

#g) weighted Bray Curtis PCoA plot

filter_prune_SL_phyloseq_bray<-vegdist(t(otu_table(filter_prune_SL_phyloseq)),method="bray",binary=FALSE,autotransform=FALSE)
filter_prune_SL_phyloseq_PCoA<-ordinate(filter_prune_SL_phyloseq,method="PCoA",filter_prune_SL_phyloseq_bray)
plot_ordination(filter_prune_SL_phyloseq,filter_prune_SL_phyloseq_PCoA,type="samples",color="Treatment")+
  geom_point(size=3)+facet_wrap(~Treatment)+ggtitle("filter_prune_SL_phyloseq weighted Bray PCoA plot")

plot_ordination(filter_prune_SL_phyloseq,filter_prune_SL_phyloseq_PCoA,type="samples",color="Treatment")+
  geom_point(size=4)+ggtitle("filter_prune_SL_phyloseq weighted Bray PCoA plot")

#h) weighted Bray Curtis PCoA plot

r_filter_prune_SL_phyloseq_bray<-vegdist(t(otu_table(r_filter_prune_SL_phyloseq)),method="bray",binary=FALSE,autotransform=FALSE)
r_filter_prune_SL_phyloseq_PCoA<-ordinate(r_filter_prune_SL_phyloseq,method="PCoA",r_filter_prune_SL_phyloseq_bray)
plot_ordination(r_filter_prune_SL_phyloseq,r_filter_prune_SL_phyloseq_PCoA,type="samples",color="Treatment")+
  geom_point(size=3)+facet_wrap(~Treatment)+ggtitle("r_filter_prune_SL_phyloseq weighted Bray PCoA plot")

plot_ordination(r_filter_prune_SL_phyloseq,r_filter_prune_SL_phyloseq_PCoA,type="samples",color="Treatment")+
  geom_point(size=4)+ggtitle("r_filter_prune_SL_phyloseq weighted Bray PCoA plot")
