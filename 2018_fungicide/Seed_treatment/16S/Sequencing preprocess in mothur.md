##   2018 fungicide seed treatment - 16S rRNA sequencing 


---
In total, there are 164 samples collected from stem canker soybean field planted with soybean seeds treated with different fungicide- Cruiser_maxx, Evergol_Energy as well as non-fungicide control.


**Mothur batch file**

```
make.contigs(inputdir=/staton/projects/soybean_rhizosphere/2018_fungicide/seed_treatment/16S/00_raw_fastq,outputdir=/staton/projects/soybean_rhizosphere/2018_fungicide/seed_treatment/16S/02_mothur,file=Seed.file,oligos=Seed.oligo,processors=30)

summary.seqs(inputdir=/staton/projects/soybean_rhizosphere/2018_fungicide/seed_treatment/16S/02_mothur,fasta=Seed.trim.contigs.fasta,processors=30)

screen.seqs(fasta=Seed.trim.contigs.fasta,group=Seed.contigs.groups,maxambig=0,maxlength=428,processors=30)

summary.seqs(fasta=Seed.trim.contigs.good.fasta,processors=30)

unique.seqs(fasta=Seed.trim.contigs.good.fasta)

count.seqs(name=Seed.trim.contigs.good.names,group=Seed.contigs.good.groups)

summary.seqs(count=Seed.trim.contigs.good.count_table,processors=30)

align.seqs(fasta=Seed.trim.contigs.good.unique.fasta,reference=silva_nr132_V3_V4.align,processors=30)

summary.seqs(fasta=Seed.trim.contigs.good.unique.align,count=Seed.trim.contigs.good.count_table,processors=30)

screen.seqs(fasta=Seed.trim.contigs.good.unique.align,count=Seed.trim.contigs.good.count_table,summary=Seed.trim.contigs.good.unique.summary,start=2,end=17012,maxhomop=8,processors=30)

summary.seqs(fasta=Seed.trim.contigs.good.unique.good.align,count=Seed.trim.contigs.good.good.count_table,processors=30)

filter.seqs(fasta=Seed.trim.contigs.good.unique.good.align,vertical=T,trump=.,processors=30)

unique.seqs(fasta=Seed.trim.contigs.good.unique.good.filter.fasta,count=Seed.trim.contigs.good.good.count_table)

pre.cluster(fasta=Seed.trim.contigs.good.unique.good.filter.unique.fasta,count=Seed.trim.contigs.good.unique.good.filter.count_table,diffs=4,processors=30)

chimera.vsearch(fasta=Seed.trim.contigs.good.unique.good.filter.unique.precluster.fasta,count=Seed.trim.contigs.good.unique.good.filter.unique.precluster.count_table,dereplicate=t,processors=30)

remove.seqs(fasta=Seed.trim.contigs.good.unique.good.filter.unique.precluster.fasta,accnos=Seed.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)

summary.seqs(fasta=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta,count=Seed.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table,processors=30)

classify.seqs(fasta=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta,count=Seed.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table,reference=trainset16_022016.rdp.fasta,taxonomy=trainset16_022016.rdp.tax,cutoff=80,processors=30)

remove.lineage(fasta=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta,count=Seed.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table,taxonomy=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.taxonomy,taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)

summary.tax(taxonomy=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.pick.taxonomy,count=Seed.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)

cluster.split(fasta=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta,count=Seed.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table,taxonomy=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.pick.taxonomy,splitmethod=classify,taxlevel=4,cutoff=0.03,processors=30)

make.shared(list=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list,count=Seed.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table,label=0.03)

classify.otu(list=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list,count=Seed.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table,taxonomy=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.pick.taxonomy,label=0.03)
```

**Summary of the sequencing processing**

1. Before any screen or filter, we have 18891448 raw reads across 164 samples
2. After remove all the bad contigs (e.g., have ambiguous calling during alignment, contigs are too long. Theoretically, the length of the contigs should be 427bp), we have 15759162 reads.  
3. After unique the sequence, we got 3181283 unique sequences out of total sequence of 15759162.
4. Precluster were used for further denoise the unique sequences. Chimera were detected and removed. It resulted in 732481 unique sequences out of 14579303 total sequences.
5. After classify the sequences to taxa, all non-bacteria sequences were removed. After this screen, only 13489302 sequences left.
6. In the end, all of the sequences were clustered to 147667 OTU at 97% similarity.

**Rarefaction and singleton remove** 

-- randomized subsample of reads to the same depth is one strategy used to mitigate the sequencing depth bias between samples.
-- In addition, singletons are tend to be sequencing error instead of real unique OTU,so, before downward community analysis, these singletons were removed from the OTU count table.
-- Here are the commands used to process it:

```
sub.sample(list=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list,taxonomy=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy,count=Seed.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table,size=13021,persample=true)

remove.rare(list=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.subsample.list,count=Seed.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.subsample.count_table,nseqs=1,groups=all,bygroup=f)

make.shared(list=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.subsample.0.03.pick.list,count=Seed.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.subsample.pick.count_table,label=0.03)

classify.otu(list=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.subsample.0.03.pick.list,count=Seed.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.subsample.pick.count_table,taxonomy=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.subsample.taxonomy)

phylotype(taxonomy=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.subsample.taxonomy)

remove.rare(list=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.subsample.tx.list,count=Seed.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.subsample.count_table,nseqs=1,groups=all,bygroup=f)

classify.otu(list=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.subsample.tx.1.pick.list,count=Seed.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.subsample.pick.count_table,taxonomy=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.subsample.taxonomy)
```
-- After the rarefaction and remove singletons, we got all seed samples removed except ``Seed_CM_4`` sample. And, we can see a dramatically high bacteria diversity within this sample. We are not sure, if this is cause by some contamination on the seed. Or is was just biased by the sequencing depth. If it was not contanination, why this sample yield high sequencing depth while the other samples did not. Does this have something to do with the DNA extraction efficiency? 

**Generate the rarefaction plot to visualized the sufficience of sequencing depth**

```
rarefaction.single(shared=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.subsample.0.03.pick.shared,label=0.03,calc=sobs,iters=1000,freq=100)

rarefaction.single(shared=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared,label=0.03,calc=sobs,iters=1000,freq=100)

```

**Generate the Lefse input file using mothur command - make.lefse** - The input file is OTU clustered based on their genus level taxonomy using mothur command ``phylotype``

-- As we ask for different specific questions from the same dataset, we did subsample the samples to a smaller subset to do the differential abundance analysis.

```
make.lefse(shared=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.subsample.tx.1.pick.shared,constaxonomy=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.subsample.tx.1.pick.1.cons.taxonomy,design=Week1_compartment.design,groups=104B_1w_ML-104E_1w_ML-104R_1w_ML-105B_1w_ML-105E_1w_ML-105R_1w_ML-108B_1w_ML-108E_1w_ML-108R_1w_ML-109B_1w_ML-109E_1w_ML-109R_1w_ML-201B_1w_ML-201E_1w_ML-201R_1w_ML-202B_1w_ML-202E_1w_ML-202R_1w_ML-206B_1w_ML-206E_1w_ML-206R_1w_ML-207B_1w_ML-207E_1w_ML-207R_1w_ML-209B_1w_ML-209E_1w_ML-209R_1w_ML-301B_1w_ML-301E_1w_ML-301R_1w_ML-302B_1w_ML-302E_1w_ML-302R_1w_ML-306B_1w_ML-306E_1w_ML-309B_1w_ML-309E_1w_ML-309R_1w_ML-310B_1w_ML-310E_1w_ML-310R_1w_ML-411B_1w_ML-411E_1w_ML-411R_1w_ML,scale=totalgroup)
```








