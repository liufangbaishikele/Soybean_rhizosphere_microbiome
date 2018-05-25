#                                      Cultivar project - 16S 3rd run and 4th run



## 3rd_run

-----

#### Transfer raw data from local computer to beacon server using scp
- logout from beacon server
- Run scp(secure copy)
```
scp -r -P 22 /g/UT/I_love_my_project/3_soybean_16S_data_and_analysis_third_run/Fang_40747707 fliu21@login.beacon.nics.utk.edu:/nics/d/home/fliu21/16S_cultivar_proj
```
#### Make a raw_data direcotry inside of ``/nics/d/home/fliu21/16S_cultivar_proj`` and move all folders from Fang_40747707 to raw_data
```
mkdir raw_data
mv /Fang_40747707/* ./raw_data
ls ./raw_data

A1-49258314  B2-49259310  C3-49255318  E1-49264285  G1-49267223
A2-49265307  B3-49266218  D1-49256307  E2-49267221  G2-49260325
A3-49263305  C1-49257326  D2-49231525  F1-49265306  H1-49263304
B1-49261310  C2-49255315  D3-49257324  F2-49267222  H2-49265308
```
#### Make a raw_fastq file inside of 16S_cultivar_proj and extract all of the fastq.gz file from raw_data daugher direcotry to raw_fastq
```
cd /nics/d/home/fliu21/16S_cultivar_proj/raw_data
for file in *; do cd $file; cp * ../../raw_fastq/; cd ../; done
```
#### Uncompress fastq.gz file
```
cd /nics/d/home/fliu21/16S_cultivar_proj/raw_fastq
gunzip *
```
-----

## 4th_run

-------
SAME process like 3rd run, include:
- Transfer data from local computer to beacon server using secure copy (scp)
- Make a direcotry of 4th_run
- Inside of 4th_run directory, create analysis and raw_fastq folder.
- Move all of the .gz file into raw_fastq folder and uncompress files to fastq (gunzip)

## 3rd_and_4th_run_combine

--------
cat file1 file2 > file_combine will be used to combine 3rd and 4th run fastq to a combined fastq file
BEFORE the above command, I need to rename fastq file name from 4th_run, because they are the same as 3rd run

#### Change file name from 4th_run/raw_fastq 
e.g., mv AgB10_S13_L001_R1_001.fastq $(echo AgB10_S13_L001_R1_001.fastq | sed "s/_/-/g") , which will change fastq file from underscore format to AgB10-S13-L001-R1-001.fastq

**Create a for loop to change all fastq file from underscore format to dash format**
```
for file in *.fastq; do echo $file; mv $file $(echo $file | sed 's/_/-/g'); done
ls
AgB10-S13-L001-R1-001.fastq  AgFresh1-S6-L001-R1-001.fastq   ForB10-S18-L001-R1-001.fastq  ForFresh1-S1-L001-R1-001.fastq
AgB10-S13-L001-R2-001.fastq  AgFresh1-S6-L001-R2-001.fastq   ForB10-S18-L001-R2-001.fastq  ForFresh1-S1-L001-R2-001.fastq
AgB11-S14-L001-R1-001.fastq  AgFresh2-S7-L001-R1-001.fastq   ForB11-S19-L001-R1-001.fastq  ForFresh2-S2-L001-R1-001.fastq
AgB11-S14-L001-R2-001.fastq  AgFresh2-S7-L001-R2-001.fastq   ForB11-S19-L001-R2-001.fastq  ForFresh2-S2-L001-R2-001.fastq
AgB12-S15-L001-R1-001.fastq  AgFresh3-S8-L001-R1-001.fastq   ForB12-S20-L001-R1-001.fastq  ForFresh3-S3-L001-R1-001.fastq
AgB12-S15-L001-R2-001.fastq  AgFresh3-S8-L001-R2-001.fastq   ForB12-S20-L001-R2-001.fastq  ForFresh3-S3-L001-R2-001.fastq
AgB8-S11-L001-R1-001.fastq   AgFresh4-S9-L001-R1-001.fastq   ForB8-S16-L001-R1-001.fastq   ForFresh4-S4-L001-R1-001.fastq
AgB8-S11-L001-R2-001.fastq   AgFresh4-S9-L001-R2-001.fastq   ForB8-S16-L001-R2-001.fastq   ForFresh4-S4-L001-R2-001.fastq
AgB9-S12-L001-R1-001.fastq   AgFresh5-S10-L001-R1-001.fastq  ForB9-S17-L001-R1-001.fastq   ForFresh5-S5-L001-R1-001.fastq
AgB9-S12-L001-R2-001.fastq   AgFresh5-S10-L001-R2-001.fastq  ForB9-S17-L001-R2-001.fastq   ForFresh5-S5-L001-R2-001.fastq
```
### Link raw_fastq from 3rd_run and 4th_run to /nics/d/home/fliu21/16S_cultivar_proj/3rd_and_4th_combine

```
ln -s /nics/d/home/fliu21/16S_cultivar_proj/3rd_run/raw_fastq/*.fastq /nics/d/home/fliu21/16S_cultivar_proj/3rd_and_4th_combine
ln -s /nics/d/home/fliu21/16S_cultivar_proj/4th_run/raw_fastq/*.fastq /nics/d/home/fliu21/16S_cultivar_proj/3rd_and_4th_combine
```

**NOTES with ln**

Link is lke a shortcut of the original file or folder. 
- IF you want to link files, but accidently linked the folder. rm folder wil work for this purpose (rm folder/ will not work). 
- IF you remove a file inside of one linked folder, actually the orginal file would be removed. BE CAREFULL, only delete links you built, do not remove files inside.

### Concanate fastq file from 3rd_run and 4th_run to one.

```
cp /nics/d/home/fliu21/16S_cultivar_proj/3rd_run/raw_fastq/treatment_ID /nics/d/home/fliu21/16S_cultivar_proj/3rd_and_4th_combine
nano treatment_ID >>>
AgB10
AgB11
AgB12
AgB8
AgB9
AgFresh1
AgFresh2
AgFresh3
AgFresh4
AgFresh5
ForB10
ForB11
ForB12
ForB8
ForB9
ForFresh1
ForFresh2
ForFresh3
ForFresh4
ForFresh5
<<<
# concatenate two R1 files into one
for file in $(cat treatment_ID); do echo $file; echo $file*R1*.fastq; cat $(echo $file*R1*.fastq) > $file-R1.fastq; done # Here the dash charactor should be changed to underscore as dash is not allowed in sample name in mothur.
for file in $(cat treatment_ID); do echo $file; echo $file*R2*.fastq; cat $(echo $file*R2*.fastq) > $file-R2.fastq; done
```

### Compare read counts among 3rd_run, 4th_run and 3rd_and_4th_run

```
grep -c "^@M" *_001.fastq > 3rd_read_count
grep -c "^@M" *-001.fastq > 4th_read_count
awk -F: '{print $1 "\t" $2}' 3rd_read_count > 3rd-read-count
awk -F: '{print $1 "\t" $2}' 4th_read_count > 4th-read-count
cd /nics/d/home/fliu21/16S_cultivar_proj/3rd_and_4th_combine/combine_raw_fastq
grep -c "^@M" *.fastq > combine_read_count
awk -F: '{print $1 "\t" $2}' combine_read_count > combine-read-count
paste -d "\t" 3rd-read-count 4th-read-count combine-read-count > reads_compare.txt
```



### Customize silva reference to fit with your primer sets

- [Refer to this link](https://github.com/mothur/mothur/issues/235)

- Alternatively, you can use this [mothur blog](http://blog.mothur.org/2016/07/07/Customization-for-your-region/) as a reference

I will demonstrate this process follow the above link:

1) First download a E.coli 16S fasta file from this site-"http://www.ncbi.nlm.nih.gov/nuccore/174375?report=fasta" and transfer this file to beacon server. MAKE SURE that the sequences are in one line.

2) Creat a **pcrTest.oligos** file that contains primer information (In my case, the forward primer is 341F and reverse primer is 805R)
```
forward CCTACGGGNGGCWGCAG
reverse GACTACHVGGGTATCTAATCC
```
3) RUN following command inside of mothur:
```
pcr.seqs(fasta=E_coli_16S.fasta,oligos=pcrTest.oligos)
Output File Names: E_coli_16S.pcr.fasta

align.seqs(fasta=E_coli_16S.pcr.fasta,reference=silva.bacteria.fasta)
Output File Names: E_coli_16S.pcr.align

summary.seqs(fasta=E_coli_16S.pcr.align,processors=8)
Using 8 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        6428    23440   427     0       6       1
2.5%-tile:      6428    23440   427     0       6       1
25%-tile:       6428    23440   427     0       6       1
Median:         6428    23440   427     0       6       1
75%-tile:       6428    23440   427     0       6       1
97.5%-tile:     6428    23440   427     0       6       1
Maximum:        6428    23440   427     0       6       1
Mean:   6428    23440   427     0       6
# of Seqs:      1

Output File Names: E_coli_16S.pcr.summary
```

4) Customize silva reference based on start and end location summary

```
mothur > pcr.seqs(fasta=silva.bacteria.fasta,start=6428,end=23440,keepdots=F,processors=8)
output File Names:silva.bacteria.pcr.fasta
```
5) I changed this silva_bacteria.pcr.fasta to silva_V3_4.fasta to make more sense to me.
----

### Now ready for mothur pipeline

**Mothur install on beacon server via bioconda**
- First install anaconda
- Using anaconda cloud to install mothur software (version 1.39.5). This version of mothur support chimera.vsearch and Opit-distance-based clustering.

**Install anaconda**
```
wget https://repo.continuum.io/archive/Anaconda3-4.4.0-Linux-x86_64.sh
bash Anaconda3-4.4.0-Linux-x86_64.sh
```


1) Build cultivar_proj.file for make.contigs command
```
cd /lustre/medusa/fliu21/16S_cultivar_proj/raw_fastq
ls *R1.fastq > R1_list
ls *R2.fastq > R2_list
sed 's/\([^.*]\)_.*/\1/' R1_list > treatment_ID
paste --delimiters="\t" treatment_ID R1_list R2_list > trial.files
```
2) Make.contigs and summerize the contigs length distribution
inputdir and outputdir were used to specify input directory and output directory for make.contigs command.
```
make.contigs(inputdir=/lustre/medusa/fliu21/16S_cultivar_proj/raw_data,outputdir=/lustre/medusa/fliu21/16S_cultivar_proj/analysis,file=trial.files,processors=8)
summary.seqs(inputdir=/lustre/medusa/fliu21/16S_cultivar_proj/analysis,fasta=trial.trim.contigs.fasta,processors=8)

mothur > summary.seqs(fasta=trial.trim.contigs.fasta,processors=8)

Using 8 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       35      35      0       2       1
2.5%-tile:      1       276     276     0       4       219492
25%-tile:       1       440     440     0       5       2194920
Median:         1       452     452     0       5       4389840
75%-tile:       1       464     464     4       6       6584760
97.5%-tile:     1       471     471     20      8       8560188
Maximum:        1       552     552     301     275     8779679
Mean:   1       443.331 443.331 3.24964 5.86621
# of Seqs:      8779679
```

Based on the above summary results, I decide to use **maxlength=471 and maxambig=0** (Theoritically, the length should be 805-341+1=465)

3) screen.seqs to screen off sequences with any ambiguous bases and the maximum length is set to 471bp.

```
screen.seqs(fasta=trial.trim.contigs.fasta,group=trials.contigs.groups,maxambig=0,maxlength=471,processors=8)
summary.seqs(fasta=trial.trim.contigs.good.fasta,processors=8)
mothur > summary.seqs(fasta=trial.trim.contigs.good.fasta,processors=8)

Using 8 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       35      35      0       2       1
2.5%-tile:      1       287     287     0       4       119050
25%-tile:       1       440     440     0       5       1190491
Median:         1       452     452     0       5       2380982
75%-tile:       1       464     464     0       6       3571472
97.5%-tile:     1       466     466     0       8       4642913
Maximum:        1       471     471     0       275     4761962
Mean:   1       444.206 444.206 0       5.40731
# of Seqs:      4761962
```
After this screen process, 54.2% (4761962/8779679) sequences left

4) Unique the contigs, count sequences and summary unique sequences

```
unique.seqs(fasta=trial.trim.contigs.good.fasta,processors=8)
count.seqs(name=trial.trim.contigs.good.names,group=trial.contigs.good.groups,processors=8)
summary.seqs(count=trial.trim.contigs.good.count_table,processors=8)
Using 8 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       35      35      0       2       1
2.5%-tile:      1       287     287     0       4       119050
25%-tile:       1       440     440     0       5       1190491
Median:         1       452     452     0       5       2380982
75%-tile:       1       464     464     0       6       3571472
97.5%-tile:     1       466     466     0       8       4642913
Maximum:        1       471     471     0       275     4761962
Mean:   1       444.206 444.206 0       5.40731
# of unique seqs:       2793447
total # of seqs:        4761962

```
After unique, the sequences decreased from 4761962 to 2793447

5) Align the unique sequence to customized reference - silva_V3_V4.fasta AND summarize the alignment.

```
align.seqs(fasta=trial.trim.contigs.good.unique.fasta,reference=silva_V3_V4.fasta,processors=8)
summary.seqs(fasta=trial.trim.contigs.good.unique.align,count=trial.trim.contigs.good.count_table,processors=8)
Using 8 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        0       0       0       0       1       1
2.5%-tile:      2       17012   11      0       3       119050
25%-tile:       2       17012   401     0       5       1190491
Median:         2       17012   413     0       5       2380982
75%-tile:       2       17012   425     0       5       3571472
97.5%-tile:     16039   17012   427     0       8       4642913
Maximum:        17012   17012   450     0       37      4761962
Mean:   478.969 16877.4 398.418 0       5.01075
# of unique seqs:       2793447
total # of seqs:        4761962
```
5) Screen alignment to make sure they start and end at the same location.AND summary the screened alignments.

```

summary.seqs(fasta=trial.trim.contigs.good.unique.good.align,count=trial.trim.contigs.good.good.count_table,processors=8                    )

Using 8 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       17012   367     0       3       1
2.5%-tile:      2       17012   401     0       4       109978
25%-tile:       2       17012   401     0       5       1099773
Median:         2       17012   413     0       5       2199546
75%-tile:       2       17012   425     0       5       3299318
97.5%-tile:     2       17012   427     0       8       4289113
Maximum:        2       17012   450     0       8       4399090
Mean:   1.99999 17012   413.176 0       5.08345
# of unique seqs:       2535598
total # of seqs:        4399090

```
6) filter and unique

```
filter.seqs(fasta=trial.trim.contigs.good.unique.good.align,vertical=T,trump=.,processors=8)
unique.seqs(fasta=trial.trim.contigs.good.unique.good.filter.fasta,count=trial.trim.contigs.good.good.count_table)
```
7) pre.cluster to further denoidse sequencing error. Here our target gene is 465bp, so set this diffs=4 (i.e., we allow one basepair of ambiguitybetween unique sequences).

```
pre.cluster(fasta=trial.trim.contigs.good.unique.good.filter.unique.fasta,count=trial.trim.contigs.good.unique.good.filter.count_table,diffs=4,processors=8)
```
8) Chimera detection using chimera.vsearch. When I typped in chimera.vsearch directly inside mothur, it give me an error: mothurvsearch file not exist, and mothur requires the vsearch executable.

Using Anaconda cloud to install vsearch package first.

```
conda install -c bioconda vsearch=2.4.3
```
After this vsearch package was installed via Anaconda cloud, I can run chimera.vsearch command now.

```
chimera.vsearch(fasta=trial.trim.contigs.good.unique.good.filter.unique.precluster.fasta,count=trial.trim.contigs.good.unique.good.filter.unique.precluster.count_table,dereplicate=t,processors=8)

remove.seqs(fasta=trial.trim.contigs.good.unique.good.filter.unique.precluster.fasta,accnos=trial.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)

# Now if we sumarize the sequences, we will see that 4399090-3980420 sequences are detected as chimera and removed from the unique sequence. We could also see that after precluster process, total unique sequence number decreased from 2535598 to 326914.

summary.seqs(fasta=trial.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta,count=trial.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table,processors=8)

Using 8 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       1274    376     0       3       1
2.5%-tile:      1       1274    401     0       4       99511
25%-tile:       1       1274    401     0       5       995106
Median:         1       1274    413     0       5       1990211
75%-tile:       1       1274    425     0       5       2985316
97.5%-tile:     1       1274    427     0       8       3880910
Maximum:        2       1274    449     0       8       3980420
Mean:   1       1274    413.21  0       5.07729
# of unique seqs:       326914
total # of seqs:        3980420
```
9) Classify.seqs and remove.lineage - Now we are going to classify those good sequences to taxon based on provided trainset fasta and trainset taxonomy file. Here the cutoff parameter is used to specify the minimum bootstrp value for each taxon assignment. Wang classification method is used for assign taxon.This method also runs a bootstrapping algorithmn to find the confidence limit of the assignment by randomly choosing with replacement 1/8 of the kmers in the query and then finding the taxonomy. To remove unwanted taxons, which include Chloroplast,mitochondria,unknown,archeae,and eukaryota, here we use remove.lineage command. After this make a summary.

```
classify.seqs(fasta=trial.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta,count=trial.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table,reference=trainset9_032012.pds.fasta,taxonomy=trainset9_032012.pds.tax,cutoff=80,processors=8)
remove.lineage(fasta=trial.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta,count=trial.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table,taxonomy=trial.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy,taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
summary.tax(taxonomy=trial.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy,count=trial.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)
```
10) Now we are ready for otu clustering, yah!!
Inside of mothur, there are three ways of OTU clustering, which include

- distance matrix based clustering (this more like de-novo assembly with a similarity cutoff of OTU at 97%, to make this process more fast, mothur comes with a taxonomy based de-novo clustering. This means, first it will split sequences to smaller bins based on their taxonomy information and clustering within each bin)
- phylotype based clustering (Cluster sequences based on the taxonomy information to genus level, i.e., each genus is a individual OTU)
- phylogenetic tree based OTU clustering (This is computing consuming, but more informative for downward UniFrac or weighted-UniFrac based PCoA analysis)

In my case, I only have 20 samples, it will real fast. I tried all of those three methods out of curiosity

- Opti distance based clustering. Got 95991 OTUs.
```
cluster.split(fasta=trial.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=trial.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=trial.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, splitmethod=classify, taxlevel=4, cutoff=0.03,processors=16)

 make.shared(list=trial.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list,count=trial.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table,label=0.03)
 
classify.otu(list=trial.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list,count=trial.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table,taxonomy=trial.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy,label=0.03)
```
- Phylotype (taxonomy) based clustering. Got 704 OTUs.
```
phylotype(taxonomy=trial.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy)

make.shared(list=trial.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.list,count=trial.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table,label=1) #Here label=1 means OTU will be assigned base on genus level taxon

classify.otu(list=trial.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.list,count=trial.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table,taxonomy=trial.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy,label=1)
```
- Generate phylogenetic tree

**Interactive mode**
```
dist.seqs(fasta=trial.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta,output=lt,processors=16)
# This take 13h to finish and the generated distance file is 357G. This caused big problem of downward clearcut process. So I added a cutoff parameter.
dist.seqs(fasta=trial.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta,output=lt,cutoff=0.1,processors=16)
# here cutoff=0.1 means that any sequences that has distance larger than 0.1, we will not cluster them into OTU. We tell mothur that any distances that larger than 0.1 will not be stored. BUT this took more than 24h, which is out of the job time I requested, don't know why.
```
**Submit job - batch mode**

- First, creat a batch file, which include mothur command used to build phylogenetic tree
Creat a batch file named **phylogenetic.batch**: 

```
dist.seqs(fasta=trial.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta,output=lt,processors=16)
clearcut(phylip=trial.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.phylip.dist)
```
- Second, creat a job file named **phylogenetic.sge**

```
#PBS -S /bin/bash
#PBS -A UT-TNEDU0031
#PBS -l nodes=2,walltime=24:00:00
#PBS -N Phylogenetic.mothur

cd $PBS_O_WORKDIR
mothur phylogenetic.batch
```
The job finished dist.seqs process, but got an error:
```
Clearcut: Incorrect number of distance values on row, M04398_48_000000000-B86V3_1_1102_14991_20478. Expected 81530, and found 81529.
Clearcut: Syntax error in distance matrix at offset 22858849901.
```

After the whole pipeline in mothur, I got 3973971 reads left, which were clustered to OTUs based on opti distance matrix.

Glance of **OTU table**

```
label  Group      numOtus  Otu00001  Otu00002   ......
0.03   AgB10      95991    5         1              .
0.03   AgB11      95991    1         1              .
0.03   AgB12      95991    0         1              .
0.03   AgB8       95991    1         0              .
0.03   AgB9       95991    0         2              .
0.03   AgFresh1   95991    2         2              .
0.03   AgFresh2   95991    2         2              .
0.03   AgFresh3   95991    0         0              .
0.03   AgFresh4   95991    2         0              .
0.03   AgFresh5   95991    3         0              .
0.03   ForB10     95991    15243     5996           .
0.03   ForB11     95991    7554      3469           .
0.03   ForB12     95991    16785     6623           .
0.03   ForB8      95991    4617      1953           .
0.03   ForB9      95991    14496     4789           .
0.03   ForFresh1  95991    7205      4581           .
0.03   ForFresh2  95991    7342      4626           .
0.03   ForFresh3  95991    1552      1037           .
0.03   ForFresh4  95991    9766      6788           .
0.03   ForFresh5  95991    7328      5227           .
```
Glance of taxonomy table

```
OTU             Size    Taxonomy
Otu00001        91904   Bacteria(100);"Verrucomicrobia"(100);Spartobacteria(100);Spartobacteria_order_incertae_sedis(100);...
Otu00002        45098   Bacteria(100);"Acidobacteria"(100);Acidobacteria_Gp2(100);Acidobacteria_Gp2_order_incertae_sedis(100);...
Otu00003        39670   Bacteria(100);"Proteobacteria"(100);Alphaproteobacteria(100);Rhizobiales(100);Rhizobiales_unclassified(100);...
Otu00004        38405   Bacteria(100);"Acidobacteria"(100);Acidobacteria_Gp2(100);Acidobacteria_Gp2_order_incertae_sedis(100);...
Otu00005        35507   Bacteria(100);"Proteobacteria"(100);Alphaproteobacteria(100);Rhizobiales(100);Bradyrhizobiaceae(100);...
...
```

## NOW, OTU count table and its corresponding taxonomy file were transfered to local computer for downward analysis in R

