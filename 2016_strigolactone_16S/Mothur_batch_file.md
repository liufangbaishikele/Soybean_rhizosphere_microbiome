Bactch file

```
make.contigs(inputdir=/staton/projects/soybean_rhizosphere/2016_strigolactone/16S_2016_strigolactone/Mothur_analysis/01_raw_fastq,outputdir=/staton/projects/soybean_rhizosphere/2016_strigolactone/16S_2016_strigolactone/Mothur_analysis/02_mothur,file=strigolactone.file,oligos=strigolactone.oligo,processors=20)

summary.seqs(inputdir=/staton/projects/soybean_rhizosphere/2016_strigolactone/16S_2016_strigolactone/Mothur_analysis/02_mothur,fasta=strigolactone.trim.contigs.fasta,processors=20)

        Using 20 processors.

                        Start   End     NBases  Ambigs  Polymer NumSeqs
        Minimum:        1       31      31      0       3       1
        2.5%-tile:      1       402     402     0       4       74631
        25%-tile:       1       403     403     0       4       746309
        Median:         1       422     422     0       5       1492618
        75%-tile:       1       427     427     0       5       2238927
        97.5%-tile:     1       428     428     8       8       2910605
        Maximum:        1       514     514     135     243     2985235
        Mean:   1       416     416     0       5
        # of Seqs:      2985235

screen.seqs(fasta=strigolactone.trim.contigs.fasta,group=strigolactone.contigs.groups,maxambig=0,maxlength=428,processors=20)

summary.seqs(fasta=strigolactone.trim.contigs.good.fasta,processors=20)
        Using 20 processors.

                        Start   End     NBases  Ambigs  Polymer NumSeqs
        Minimum:        1       40      40      0       3       1
        2.5%-tile:      1       402     402     0       4       56680
        25%-tile:       1       403     403     0       4       566792
        Median:         1       422     422     0       5       1133583
        75%-tile:       1       427     427     0       5       1700374
        97.5%-tile:     1       428     428     0       8       2210485
        Maximum:        1       428     428     0       38      2267164
        Mean:   1       415     415     0       5
        # of Seqs:      2267164 

count.seqs(name=strigolactone.trim.contigs.good.names,group=strigolactone.contigs.good.groups)
align.seqs(fasta=strigolactone.trim.contigs.good.unique.fasta,reference=silva_nr132_V3_V4.align,processors=20)

summary.seqs(fasta=strigolactone.trim.contigs.good.unique.align,count=strigolactone.trim.contigs.good.count_table,processors=20)

        Using 20 processors.

                        Start   End     NBases  Ambigs  Polymer NumSeqs
        Minimum:        1       6       2       0       1       1
        2.5%-tile:      2       17012   401     0       4       56680
        25%-tile:       2       17012   402     0       4       566792
        Median:         2       17012   420     0       5       1133583
        75%-tile:       2       17012   426     0       5       1700374
        97.5%-tile:     2       17012   427     0       8       2210485
        Maximum:        17011   17012   428     0       38      2267164
        Mean:   33      17010   414     0       5
        # of unique seqs:       866103
        total # of seqs:        2267164

screen.seqs(fasta=strigolactone.trim.contigs.good.unique.align,count=strigolactone.trim.contigs.good.count_table,summary=strigolactone.trim.contigs.good.unique.summary,start=2,end=17012,maxhomop=8,processors=20)

summary.seqs(fasta=strigolactone.trim.contigs.good.unique.good.align,count=strigolactone.trim.contigs.good.good.count_table,processors=20)

        Using 20 processors.

                        Start   End     NBases  Ambigs  Polymer NumSeqs
        Minimum:        1       17012   371     0       3       1
        2.5%-tile:      2       17012   401     0       4       55469
        25%-tile:       2       17012   402     0       4       554681
        Median:         2       17012   421     0       5       1109362
        75%-tile:       2       17012   426     0       5       1664043
        97.5%-tile:     2       17012   427     0       8       2163255
        Maximum:        2       17012   428     0       8       2218723
        Mean:   1       17012   415     0       5
        # of unique seqs:       825593
        total # of seqs:        2218723

filter.seqs(fasta=strigolactone.trim.contigs.good.unique.good.align,vertical=T,trump=.,processors=20)

pre.cluster(fasta=strigolactone.trim.contigs.good.unique.good.filter.unique.fasta,count=strigolactone.trim.contigs.good.unique.good.filter.count_table,diffs=4,processors=20)

chimera.vsearch(fasta=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.fasta,count=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.count_table,dereplicate=t,processors=20)

remove.seqs(fasta=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.fasta,accnos=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)

summary.seqs(fasta=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta,count=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table,processors=20)

classify.seqs(fasta=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta,count=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table,reference=trainset16_022016.rdp.fasta,taxonomy=trainset16_022016.rdp.tax,cutoff=80,processors=20)

remove.lineage(fasta=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta,count=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table,taxonomy=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.taxonomy,taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)

summary.tax(taxonomy=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.pick.taxonomy,count=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)

cluster.split(fasta=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta,count=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table,taxonomy=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.pick.taxonomy,splitmethod=classify,taxlevel=4,cutoff=0.03,processors=20)

make.shared(list=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list,count=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table,label=0.03)

classify.otu(list=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list,count=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table,taxonomy=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.pick.taxonomy,label=0.03)
```
