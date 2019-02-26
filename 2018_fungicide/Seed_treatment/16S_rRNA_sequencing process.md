##   2018 fungicide seed treatment - 16S rRNA sequencing 

----
In total, there are 164 samples collected from stem canker soybean field planted with soybean seeds treated with different fungicide- Cruiser_maxx, Evergol_Energy as well as non-fungicide control.
----

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
2. After remove all the bad contigs (e.g., have ambiguous calling during alignment, contigs are too long). Theoretically, the length of the contigs should be 428bp.  
3. 





