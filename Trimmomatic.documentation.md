# Demultiplexed and trimmed reads are further trimmed using trimmomatic software using SLIDINGWINDOW:4:15

**trimmomatic job file**

```
#$ -N Trimmomatic_cultivar
#$ -q medium*
#$ -t 1-139
#$ -cwd
#$ -S /bin/bash

R1=`head -$SGE_TASK_ID <(ls /lustre/projects/staton/projects/soybean_strigolactones/16S_raw_fastq/fastq_gunzip/trimmomatic_cultivar_QC/trimmed_cultivar_demultiplex.fas$
R2=`head -$SGE_TASK_ID <(ls /lustre/projects/staton/projects/soybean_strigolactones/16S_raw_fastq/fastq_gunzip/trimmomatic_cultivar_QC/trimmed_cultivar_demultiplex.fas$
BASE=`sed 's/trimmed_\(.*\)_R1.fastq/\1/' R1_list | head -$SGE_TASK_ID | tail -1`

module load trimmomatic/0.36
trimmomatic \
 PE \
 -threads 1 \
 -trimlog .log \
 -phred33 \
  $R1 \
  $R2 \
  $BASE.R1.paired.trimmomatic.fastq \
  $BASE.R1.unpaired.trimmomatic.fastq \
  $BASE.R2.paired.trimmomatic.fastq \
  $BASE.R2.unpaired.trimmomatic.fastq \
 SLIDINGWINDOW:4:15 MINLEN:30 \
 >& $BASE.trim_output.txt

```
- make.contigs (.........)
*After make.contigs, here is the summary, we got lots of ambigous base call that assigned to be N*

## Before screen.seqs


**Strigolactone projecct -- 48 samples, V3_V4 region**
- summary.seqs
```
Using 8 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       3       3       0       1       1
2.5%-tile:      1       176     176     0       4       32573
25%-tile:       1       313     313     0       4       325727
Median:         1       361     361     5       5       651453
75%-tile:       1       404     404     13      6       977179
97.5%-tile:     1       429     429     28      6       1270333
Maximum:        1       560     560     58      295     1302905
Mean:   1       350.867 350.867 8.20351 4.92022
# of Seqs:      1302905

```
 **Popolus project -- 12 samples, V4 region**
 ```
 mothur > summary.seqs(fasta=treatment.trim.contigs.fasta,processors=8)

Using 8 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       217     217     0       3       1
2.5%-tile:      1       253     253     0       4       25467
25%-tile:       1       253     253     0       4       254667
Median:         1       253     253     0       5       509334
75%-tile:       1       253     253     0       5       764001
97.5%-tile:     1       381     381     2       8       993201
Maximum:        1       502     502     71      251     1018667
Mean:   1       259.284 259.284 0.312466        4.91969
# of Seqs:      1018667
 ```
 **Xiaolong project -- 15 sample V3_V4 region**
 
 ```
 mothur > summary.seqs(fasta=column.trim.contigs.fasta,processors=8)

Using 8 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       297     297     0       3       1
2.5%-tile:      1       439     439     0       4       301091
25%-tile:       1       442     442     0       5       3010909
Median:         1       464     464     1       5       6021818
75%-tile:       1       465     465     4       6       9032727
97.5%-tile:     1       466     466     22      7       11742545
Maximum:        1       602     602     285     301     12043635
Mean:   1       456.39  456.39  3.5451  5.28368
# of Seqs:      12043635 ```
 ```
 
- screen.seqs(fasta=strigolactone.trim.contigs.fasta,group=strigolactone.contigs.groups,maxambig=0,maxlength=429)


## After screen.seqs


**Strigolactone project -- 33.7% contigs left**
```
mothur > summary.seqs(fasta=strigolactone.trim.contigs.good.fasta,processors=8)

Using 8 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       10      10      0       2       1
2.5%-tile:      1       264     264     0       4       10988
25%-tile:       1       404     404     0       4       109876
Median:         1       421     421     0       5       219752
75%-tile:       1       429     429     0       5       329627
97.5%-tile:     1       429     429     0       6       428515
Maximum:        1       429     429     0       206     439502
Mean:   1       408.788 408.788 0       4.89575
# of Seqs:      439502
```



**Populus project -- 926176/1018667= 90.9% contigs left**
```
mothur > summary.seqs(fasta=treatment.trim.contigs.good.fasta,processors=8)

Using 8 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       253     253     0       3       1
2.5%-tile:      1       253     253     0       4       23155
25%-tile:       1       253     253     0       4       231545
Median:         1       253     253     0       4       463089
75%-tile:       1       253     253     0       5       694633
97.5%-tile:     1       254     254     0       7       903022
Maximum:        1       292     292     0       32      926176
Mean:   1       253.316 253.316 0       4.63625
# of Seqs:      926176
```
**Xiaolong project -- 5614928/12043635=46.6% contigs left**

```
mothur > summary.seqs(fasta=column.trim.contigs.good.fasta,processors=8)

Using 8 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       301     301     0       3       1
2.5%-tile:      1       440     440     0       4       140374
25%-tile:       1       442     442     0       5       1403733
Median:         1       464     464     0       5       2807465
75%-tile:       1       465     465     0       6       4211197
97.5%-tile:     1       466     466     0       7       5474555
Maximum:        1       466     466     0       74      5614928
Mean:   1       455.839 455.839 0       5.22675
# of Seqs:      5614928

```
- unique.seqs(fasta=strigolactone.trim.contigs.good.fasta)
- count.seqs(name=strigolactone.trim.contigs.good.names,group=strigolactone.contigs.good.groups)
- summary.seqs(count=strigolactone.trim.contigs.good.count_table,processors=8)
```
mothur > summary.seqs(count=strigolactone.trim.contigs.good.count_table,processors=8)
Using strigolactone.trim.contigs.good.unique.fasta as input file for the fasta parameter.

Using 8 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       10      10      0       2       1
2.5%-tile:      1       264     264     0       4       10988
25%-tile:       1       404     404     0       4       109876
Median:         1       421     421     0       5       219752
75%-tile:       1       429     429     0       5       329627
97.5%-tile:     1       429     429     0       6       428515
Maximum:        1       429     429     0       206     439502
Mean:   1       408.788 408.788 0       4.89575
# of unique seqs:       255207
total # of seqs:        439502

```

- align.seqs(fasta=strigolactone.trim.contigs.good.unique.fasta,reference=silva_V3_V4.fasta,processors=8)

```
mothur > summary.seqs(fasta=strigolactone.trim.contigs.good.unique.align,count=strigolactone.trim.contigs.good.count_table,processors=8)

Using 8 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        0       0       0       0       1       1
2.5%-tile:      2       9993    258     0       4       10988
25%-tile:       2       17016   403     0       4       109876
Median:         2       17016   417     0       5       219752
75%-tile:       2       17016   428     0       5       329627
97.5%-tile:     1379    17016   428     0       6       428515
Maximum:        17016   17016   429     0       110     439502
Mean:   108.802 16734.7 407.344 0       4.89088
# of unique seqs:       255207
total # of seqs:        439502
```

- screen.seqs(fasta=strigolactone.trim.contigs.good.unique.align,count=strigolactone.trim.contigs.good.count_table,summary=strigolactone.trim.contigs.good.unique.summary,start=2,end=17016,maxhomop=8,processors=8)
- summary.seqs(fasta=strigolactone.trim.contigs.good.unique.good.align,count=strigolactone.trim.contigs.good.good.count_table,processors=8)

```
mothur > summary.seqs(fasta=strigolactone.trim.contigs.good.unique.good.align,count=strigolactone.trim.contigs.good.good.count_table,processors=8)

Using 8 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        2       17016   359     0       3       1
2.5%-tile:      2       17016   403     0       4       10015
25%-tile:       2       17016   403     0       4       100150
Median:         2       17016   427     0       5       200300
75%-tile:       2       17016   428     0       5       300449
97.5%-tile:     2       17016   428     0       6       390584
Maximum:        2       17016   429     0       8       400598
Mean:   2       17016   417.021 0       4.88047
# of unique seqs:       219152
total # of seqs:        400598
```

- filter.seqs(fasta=strigolactone.trim.contigs.good.unique.good.align,vertical=T,trump=.,processors=8)
- unique.seqs(fasta=strigolactone.trim.contigs.good.unique.good.filter.fasta,count=strigolactone.trim.contigs.good.good.count_table)
- summary.seqs(fasta=strigolactone.trim.contigs.good.unique.good.filter.unique.fasta,count=
- mothur > summary.seqs(fasta=strigolactone.trim.contigs.good.unique.good.filter.unique.fasta,count=strigolactone.trim.contigs.good.unique.good.filter.count_table,processors=8)

```
                          Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       956     359     0       3       1
2.5%-tile:      1       956     403     0       4       10015
25%-tile:       1       956     403     0       4       100150
Median:         1       956     427     0       5       200300
75%-tile:       1       956     428     0       5       300449
97.5%-tile:     1       956     428     0       6       390584
Maximum:        1       956     429     0       8       400598
Mean:   1       956     417.021 0       4.88047
# of unique seqs:       218809
total # of seqs:        400598

```
-  pre.cluster(fasta=strigolactone.trim.contigs.good.unique.good.filter.unique.fasta,count=strigolactone.trim.contigs.good.unique.good.filter.count_table,diffs=4,processors=8)

- summary.seqs(fasta=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.fasta,count=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.count_table,processors=8)

```
Using 8 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       956     359     0       3       1
2.5%-tile:      1       956     403     0       4       10015
25%-tile:       1       956     403     0       4       100150
Median:         1       956     428     0       5       200300
75%-tile:       1       956     428     0       5       300449
97.5%-tile:     1       956     428     0       6       390584
Maximum:        1       956     429     0       8       400598
Mean:   1       956     417.037 0       4.8737
# of unique seqs:       96996
total # of seqs:        400598

```
-  chimera.vsearch(fasta=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.fasta,count=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.count_table,dereplicate=t,processors=8)

- remove.seqs(fasta=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.fasta,accnos=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)

- summary.seqs(fasta=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta,count=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table,processors=8)

```
Using 8 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       956     359     0       4       1
2.5%-tile:      1       956     403     0       4       8957
25%-tile:       1       956     403     0       4       89564
Median:         1       956     428     0       5       179127
75%-tile:       1       956     428     0       5       268690
97.5%-tile:     1       956     428     0       6       349297
Maximum:        1       956     429     0       8       358253
Mean:   1       956     417.259 0       4.8438
# of unique seqs:       61227
total # of seqs:        358253

```
- cluster.split(fasta=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta,count=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table,taxonomy=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy,splitmethod=classify,taxlevel=4,cutoff=0.15,processors=8)

-  make.shared(list=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list,count=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table,label=0.03)

- classify.otu(list=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list,count=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table,taxonomy=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy,label=0.03)

Here we got 11631 OTUs.
The OTU table file: strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared
 
Corresponding taxonomy file of each OTU to genus level: 
strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy
 
Both of the OTU table and taxomoy file were transfered to local computer using FileZilla for further analysis in R


 




