```
make.contigs(inputdir=/nics/d/home/fliu21/16S_cultivar_proj/raw_fastq,outputdir=/nics/d/home/fliu21/16S_cultivar_proj/analysis/002_mothur_analysis/,file=cultivar.file,oligos=cultivar.oligo,processors=16)

summary.seqs(inputdir=/nics/d/home/fliu21/16S_cultivar_proj/analysis/002_mothur_analysis/,fasta=cultivar.trim.contigs.fasta,processors=16)

screen.seqs(fasta=cultivar.trim.contigs.fasta,group=cultivar.contigs.groups,maxambig=0,maxlength=428,processors=16)

summary.seqs(fasta=cultivar.trim.contigs.good.fasta,processors=16)

unique.seqs(fasta=cultivar.trim.contigs.good.fasta)

count.seqs(name=cultivar.trim.contigs.good.names,group=cultivar.contigs.good.groups,processors=16)

summary.seqs(count=cultivar.trim.contigs.good.count_table,processors=16)

align.seqs(fasta=cultivar.trim.contigs.good.unique.fasta,reference=silva_V3_4.fasta,processors=16)

summary.seqs(fasta=cultivar.trim.contigs.good.unique.align,count=cultivar.trim.contigs.good.count_table,processors=16)

screen.seqs(fasta=cultivar.trim.contigs.good.unique.align,count=cultivar.trim.contigs.good.count_table,summary=cultivar.trim.contigs.good.unique.summary,start=2,end=17012,maxhomop=8,processors=16)

summary.seqs(fasta=cultivar.trim.contigs.good.unique.good.align,count=cultivar.trim.contigs.good.good.count_table,processors=16)

filter.seqs(fasta=cultivar.trim.contigs.good.unique.good.align,vertical=T,trump=.,processors=16)

#unique.seqs(fasta=cultivar.trim.contigs.good.unique.good.filter.fasta,count=cultivar.trim.contigs.good.good.count_table,processors=8)
#processors is not a valid parameter.

unique.seqs(fasta=cultivar.trim.contigs.good.unique.good.filter.fasta,count=cultivar.trim.contigs.good.good.count_table)

pre.cluster(fasta=cultivar.trim.contigs.good.unique.good.filter.unique.fasta,count=cultivar.trim.contigs.good.unique.good.filter.count_table,diffs=4,processors=16)

chimera.vsearch(fasta=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.fasta,count=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.count_table,dereplicate=t,processors=16)

remove.seqs(fasta=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.fasta,accnos=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)

summary.seqs(fasta=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta,count=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table,processors=16)

classify.seqs(fasta=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta,count=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table,reference=trainset9_032012.pds.fasta,taxonomy=trainset9_032012.pds.tax,cutoff=80,processors=16)

remove.lineage(fasta=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta,count=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table,taxonomy=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy,taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)

summary.tax(taxonomy=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy,count=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)

cluster.split(fasta=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta,count=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table,taxonomy=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy,splitmethod=classify,taxlevel=4,cutoff=0.03,processors=16)

make.shared(list=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list,count=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table,label=0.03)

classify.otu(list=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list,count=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table,taxonomy=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy,label=0.03)
```

