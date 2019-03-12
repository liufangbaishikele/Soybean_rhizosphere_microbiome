## ITS sequence process using mothur

----
The sequence process pipeline of ITS is different from 16S
----

**Mothur batch file**

```
make.contigs(inputdir=/staton/projects/soybean_rhizosphere/2018_fungicide/seed_treatment/ITS/00_raw_fastq,outputdir=/staton/projects/soybean_rhizosphere/2018_fungicide/seed_treatment/ITS/02_mothur,file=seed.file,oligos=seed.oligo,processors=30)

summary.seqs(inputdir=/staton/projects/soybean_rhizosphere/2018_fungicide/seed_treatment/ITS/02_mothur,fasta=seed.trim.contigs.fasta,processors=30)

screen.seqs(fasta=seed.trim.contigs.fasta,group=seed.contigs.groups,maxambig=0,maxlength=462,processors=30)

summary.seqs(fasta=seed.trim.contigs.good.fasta,processors=30)

unique.seqs(fasta=seed.trim.contigs.good.fasta)

count.seqs(name=seed.trim.contigs.good.names,group=seed.contigs.good.groups)

summary.seqs(count=seed.trim.contigs.good.count_table,processors=30)

pre.cluster(fasta=seed.trim.contigs.good.unique.fasta,count=seed.trim.contigs.good.count_table,diffs=3,processors=30)

summary.seqs(fasta=seed.trim.contigs.good.unique.precluster.fasta,count=seed.trim.contigs.good.unique.precluster.count_table)

chimera.vsearch(fasta=seed.trim.contigs.good.unique.precluster.fasta,count=seed.trim.contigs.good.unique.precluster.count_table,dereplicate=t,processors=30)

remove.seqs(fasta=seed.trim.contigs.good.unique.precluster.fasta,accnos=seed.trim.contigs.good.unique.precluster.denovo.vsearch.accnos)

summary.seqs(fasta=seed.trim.contigs.good.unique.precluster.pick.fasta,count=seed.trim.contigs.good.unique.precluster.denovo.vsearch.pick.count_table,processors=30)

classify.seqs(fasta=seed.trim.contigs.good.unique.precluster.pick.fasta,count=seed.trim.contigs.good.unique.precluster.denovo.vsearch.pick.count_table,reference=UNITEv6_sh_97.fasta,taxonomy=UNITEv6_sh_97.tax,cutoff=80,processors=30)

remove.lineage(fasta=seed.trim.contigs.good.unique.precluster.pick.fasta,count=seed.trim.contigs.good.unique.precluster.denovo.vsearch.pick.count_table,taxonomy=seed.trim.contigs.good.unique.precluster.pick.UNITEv6_sh_97.wang.taxonomy,taxon=unknown-k__Metazoa-k__Plantae-k__Rhizaria-k__Chromista)

summary.tax(taxonomy=seed.trim.contigs.good.unique.precluster.pick.UNITEv6_sh_97.wang.pick.taxonomy,count=seed.trim.contigs.good.unique.precluster.denovo.vsearch.pick.pick.count_table)

pairwise.seqs(fasta=seed.trim.contigs.good.unique.precluster.pick.pick.fasta,cutoff=0.1,align=needleman,output=lt,countends=T,calc=onegap,processors=30)

cluster(phylip=seed.trim.contigs.good.unique.precluster.pick.pick.phylip.dist,count=seed.trim.contigs.good.unique.precluster.denovo.vsearch.pick.pick.count_table,method=opti,cutoff=0.03,processors=30)

#cluster.split(fasta=seed.trim.contigs.good.unique.precluster.pick.pick.fasta,count=seed.trim.contigs.good.unique.precluster.denovo.vsearch.pick.pick.count_table,taxonomy=seed.trim.contigs.good.unique.precluster.pick.UNITEv6_sh_97.wang.pick.taxonomy,splitmethod=classify,taxlevel=4,cutoff=0.03,processors=30) Actually, I tried this command but it give me error with the reason that sequences are not at the same length

make.shared(list=seed.trim.contigs.good.unique.precluster.pick.pick.phylip.opti_mcc.list,count=seed.trim.contigs.good.unique.precluster.denovo.vsearch.pick.pick.count_table,label=0.03)

classify.otu(list=seed.trim.contigs.good.unique.precluster.pick.pick.phylip.opti_mcc.list,count=seed.trim.contigs.good.unique.precluster.denovo.vsearch.pick.pick.count_table,taxonomy=seed.trim.contigs.good.unique.precluster.pick.UNITEv6_sh_97.wang.pick.taxonomy,label=0.03)

#Phylotype based clustering

#phylotype(taxonomy=seed.trim.contigs.good.unique.precluster.pick.UNITEv6_sh_97.wang.pick.taxonomy)

#make.shared(list=seed.trim.contigs.good.unique.precluster.pick.UNITEv6_sh_97.wang.pick.tx.list,count=seed.trim.contigs.good.unique.precluster.denovo.vsearch.pick.pick.count_table,label=1)

#classify.otu(list=seed.trim.contigs.good.unique.precluster.pick.UNITEv6_sh_97.wang.pick.tx.list,count=seed.trim.contigs.good.unique.precluster.denovo.vsearch.pick.pick.count_table,taxonomy=seed.trim.contigs.good.unique.precluster.pick.UNITEv6_sh_97.wang.pick.taxonomy,label=1)
```

** Processing notes

1. Before any screening or filtering, we have 14146478 reads across 164 samples.
2. After the first step screening based on ambiguous calling, have 11527802 reads left
3. After unique process, above sequence were represented by 846518 unique sequence
4. pre.cluster were used as a further step of condense unique sequence and eliminate sequencing error. After precluter, unique sequence were reduced to 266831.
5. Vsearch was used to detect chimera. All chimera were removed. Left 11499406 reads left
6. After remove non-fungi sequences, got 3981643 sequences left.
7. Stuck at this step ... ``pairwise.seqs`` took forever.






