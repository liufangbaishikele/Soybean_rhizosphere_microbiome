


# Soybean rhizosphere microbiome-16S read analysis  #

*First MiSeq run*

- The first run yields 5,312,110 reads, with 3,223,228 belongs to cultivar project and 1,309,647 belongs to strigolactone project. On average, 71.6% of the quality score is higher than 30.  After further demultiplexing, we got 4,541,875 reads left.
                                                        
  

##                                                    Sequence decompress

```
cd /lustre/projects/staton/projects/soybean_strigolactones/16S_raw_fastq/soybean_strigolactone_16S_03_02_2017
for file in *; do echo $file; cd $file; cp * ../../fastq_gunzip/; cd ..; done
cd /lustre/projects/staton/projects/soybean_strigolactones/16S_raw_fastq/fastq_gunzip
gunzip *

```

When submit sample sheet, the format looks like this one:

Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
SLESCCol11_SLESCrar34,SAMPLE1,16Spool4,,A001,CGTCGGTAA,,
SLESCCol12_SLESCrar35,SAMPLE2,16Spool4,,A002,GTGTCCAAT,,
Tips to make file arrangement easier

I7_index_ID need to be reverse complemented from Lundberg material. Otherwise, these I7_index could not be recognized during demultiplexing.
Sample_ID shoud be the combined sample name that barcoded with same reverse barcode. In this form: sampletreat#_TGA_sampletreat#_ACT This sample ID will be the directory name of each R1 and R2 read files.
Sample_Name will be part of R1 and R2 file name.


##                                                       Demultiplex 


### Library prep protocol refer to ### 

>Lundberg, D.S., Yourstone, S., Mieczkowski, P., Jones, C.D. and Dangl, J.L., 2013. Practical innovations for high-throughput amplicon sequencing. Nature methods, 10(10), pp.999-1002.

### Barcode1 ###
>F-Bc1_Fs2,4,6 seqeunce is *TGA....TCACTCCTACGGGNGGCWGCAG*

### Barcode2 ###
>F-Bc2_Fs1,3,5 sequence is *ACT....TCACTCCTACGGGNGGCWGCAG*

### demultiplex script and jobs ###

**1) Read1 demultiplex scripts**


1) Read1 barcoded with ACT
```
# "$1" is read1.fastq, "$2" is sampleID that barcoded with TGA, "$3" is sampleID that barcoded with ACT
grep -B1 -A2 "ACT....TCACTCCTACGGG[ATGC]GGC[AT]GCAG" "$1" | grep -v "^--$" > "$3"_R1.fastq
grep "^@M04398" "$3"_R1.fastq > "$3"_header1.txt
awk '{print$1}' "$3"_header1.txt > "$3"_header2.txt
```
2) Read1 barcoded with TGA

```
# "$1" is read1.fastq, "$2" is sampleID that barcoded with TGA, "$3" is sampleID that barcoded with TGA
grep -B1 -A2 "TGA....TCACTCCTACGGG[ATGC]GGC[AT]GCAG" "$1" | grep -v "^--$" > "$2"_R1.fastq
grep "^@M04398" "$2"_R1.fastq > "$2"_header1.txt
awk '{print$1}' "$2"_header1.txt > "$2"_header2.txt

```
Read1 demultiplex job file

```
#$ -N demultiplex
#$ -cwd
#$ -q medium*

sh  demultiplex_update.sh  SAMPLE1_S1_L001_R1_001.fastq    Gene1_01      AgCV1_01
sh  demultiplex_update.sh  SAMPLE2_S2_L001_R1_001.fastq    Gene10_01     AgCV1_02
sh  demultiplex_update.sh  SAMPLE3_S3_L001_R1_001.fastq    Gene14_01     AgCV1_03
sh  demultiplex_update.sh  SAMPLE4_S4_L001_R1_001.fastq    CT_01         AgCV1_04
sh  demultiplex_update.sh  SAMPLE5_S5_L001_R1_001.fastq    HypIII_01     AgCV1_05
sh  demultiplex_update.sh  SAMPLE6_S6_L001_R1_001.fastq    Ag_B_01       AgCV1_06
sh  demultiplex_update.sh  SAMPLE7_S7_L001_R1_001.fastq    For_B_01      AgCV1_07
sh  demultiplex_update.sh  SAMPLE8_S8_L001_R1_001.fastq    AgCV1_09      AgCV1_08
sh  demultiplex_update.sh  SAMPLE9_S9_L001_R1_001.fastq    Gene1_02      AgCV2_01
sh  demultiplex_update.sh  SAMPLE10_S10_L001_R1_001.fastq  Gene10_02     AgCV2_02

          ...                          ...                    ...           ...
          ...                          ...                    ...           ...
          ...                          ...                    ...           ...
          
sh  demultiplex_update.sh  SAMPLE95_S95_L001_R1_001.fastq  For_B_12     ForCV6_07
sh  demultiplex_update.sh  SAMPLE96_S96_L001_R1_001.fastq  Blank        ForCV6_08
    
```
Here are some tricky stuffs I got, because dos and unix character are different.
Every file named by the fifth column of job file, got output with extra character, e.g., AgCV1_01.fastq became AgCV1_01?.fastq
But those file names could not be recognized by bash command, mv does not work. So, here came with the solution 

```
for file in *; do mv "$file" "$(echo $file | sed s'/\?r//g')"; done
```

To avoid this error (told by Miriam):

1) check job file using head | od -c to check if there are any \r characters

2) transfer file from dos to unix form by 

```dos2unix filename
```

**2) Read2 demultiplex**

Using unpublished lab code - *extract_seq_list_from_fasta_file.py*

Job file

```
#$ -N Read2_extract
#$ -cwd
#$ -q medium*

module load python/3.5.1
module load biopython/1.65


python  extract_seq_list_from_fasta_file.py     AgCV1_1_header2.txt     SAMPLE1_S1_L001_R2_001.fastq
python  extract_seq_list_from_fasta_file.py     AgCV1_2_header2.txt     SAMPLE2_S2_L001_R2_001.fastq
python  extract_seq_list_from_fasta_file.py     AgCV1_3_header2.txt     SAMPLE3_S3_L001_R2_001.fastq
python  extract_seq_list_from_fasta_file.py     AgCV1_4_header2.txt     SAMPLE4_S4_L001_R2_001.fastq
python  extract_seq_list_from_fasta_file.py     AgCV1_5_header2.txt     SAMPLE5_S5_L001_R2_001.fastq
python  extract_seq_list_from_fasta_file.py     AgCV1_6_header2.txt     SAMPLE6_S6_L001_R2_001.fastq
python  extract_seq_list_from_fasta_file.py     AgCV1_7_header2.txt     SAMPLE7_S7_L001_R2_001.fastq
python  extract_seq_list_from_fasta_file.py     AgCV1_8_header2.txt     SAMPLE8_S8_L001_R2_001.fastq
python  extract_seq_list_from_fasta_file.py     AgCV2_1_header2.txt     SAMPLE9_S9_L001_R2_001.fastq
                  ...                                  ...                      ...
                  ...                                  ...                      ...
                  ...                                  ...                      ...

python  extract_seq_list_from_fasta_file.py     ForCV6_3_header2.txt    SAMPLE91_S91_L001_R2_001.fastq
python  extract_seq_list_from_fasta_file.py     ForCV6_4_header2.txt    SAMPLE92_S92_L001_R2_001.fastq
python  extract_seq_list_from_fasta_file.py     ForCV6_5_header2.txt    SAMPLE93_S93_L001_R2_001.fastq
python  extract_seq_list_from_fasta_file.py     ForCV6_6_header2.txt    SAMPLE94_S94_L001_R2_001.fastq
python  extract_seq_list_from_fasta_file.py     ForCV6_7_header2.txt    SAMPLE95_S95_L001_R2_001.fastq
python  extract_seq_list_from_fasta_file.py     ForCV6_8_header2.txt    SAMPLE96_S96_L001_R2_001.fastq               

```
output are named like this 

*SAMPLE1_S1_L001_R2_001.fastq.filtered*

So, file name were then renamed from corresponding filter file to Read2 file

e.g., *mv SAMPLE1_S1_L001_R2_001.fastq.filtered   AgCV1_01_R2.fastq*

At the end, all of the demultiplexed read1 and corresponding read2 file were checked to make sure they have same number of lines after extraction

## Trimming primer, frameshift and molecule tag sequence using cutadapt2

1) Read1 trimming

```
conda install cutadapt -y
for line in $(cat R1_ID_for_cutadapt)
do
        echo $line
        cutadapt -g ACTCCTACGGGNGGCWGCAG --overlap=20 -o trimmed_$line  $line
done

```
2) Read2 trimming

```
for line in $(cat R2_ID_for_cutadapt)
do
        echo $line
        cutadapt -g GGACTACHVGGGTWTCTAAT --overlap=20 -o trimmed_$line  $line
done

```

## Read file arrangement ##

All of the demultiplexed and trimmed reads from cultivar project are saved in

``/lustre/projects/staton/projects/soybean_strigolactones/16S_raw_fastq/fastq_gunzip/trimmed_cultivar_demultiplex.fastq``


While those from strigolactone project are saved in 

``/lustre/projects/staton/projects/soybean_strigolactones/16S_raw_fastq/fastq_gunzip/trimmed_strigolactone_demultiplex.fastq``

**Now read files are ready for subsequent analysis using mothur**


##                                                  Mothur--sequence processing- Cultivar project##

Sequence preprocess refer to [MiSeq Sop] (https://www.mothur.org/wiki/MiSeq_SOP)

###                                                      Procedure ###

- 01_make.file

``ls *_R1.fastq | column -t >R1_list``

``ls *_R2.fastq | column -t >R2_list``

``sed 's/trimmed_\(.*\)_R1.fastq/\1/' R1_list > sample_ID``

``paste --delimiters= " " sample_ID R1_list R2_list | column -t > cultivar.file``

** All of file name generated later begin with cultivar**
** mv all fastq file and cultivar.file into trimmed_raw_read directory** 

- 02_make.contigs and summary

``make.contigs(inputdir=/lustre/projects/staton/projects/soybean_strigolactones/16S_raw_fastq/fastq_gunzip/trimmed_cultivar_demultiplex.fastq/trimmed_raw_read,outputdir=/lustre/projects/staton/projects/soybean_strigolactones/16S_raw_fastq/fastq_gunzip/trimmed_cultivar_demultiplex.fastq/mothur,file=cultivar.file,processors=8)``

``summary.seqs(fasta=cultivar.trim.contigs.fasta,processors=8)`` 
```
After make.contigs:
 Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       267     267     0       3       1
2.5%-tile:      1       403     403     0       4       80657
25%-tile:       1       409     409     0       4       806565
Median:         1       429     429     1       5       1613130
75%-tile:       1       429     429     4       6       2419694
97.5%-tile:     1       430     430     16      8       3145602
Maximum:        1       569     569     268     300     3226258
Mean:   1       422.28  422.28  2.95994 5.10718

```
**number of Seqs: 3226258**

Most of the sequence with length of 429 and 430,with 429bp (806-338+1-20-20)being our expected length. There are lots of sequences with ambiguous base more than 4bp, this indicate that a long target region may easy to introduce sequence error even with V3 300X2 pair_end sequencing

Therefore, in the following screen process, I changed ``maxambig`` parameter from ``0`` to ``3``, meanwhile, ``maxlength`` changed to ``450``

- 03_screen.seqs
```
screen.seqs(fasta=cultivar.trim.contigs.fasta,group=cultivar.contigs.groups,summary=cultivar.trim.contigs.summary,maxambig=3,maxlength=450,processors=8)
```
```
summary.seqs(fasta=cultivar.trim.contigs.good.fasta,processors=8)
```
```
Using 8 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       270     270     0       3       1
2.5%-tile:      1       404     404     0       4       58233
25%-tile:       1       406     406     0       4       582329
Median:         1       429     429     0       5       1164658
75%-tile:       1       429     429     1       6       1746987
97.5%-tile:     1       430     430     3       6       2271083
Maximum:        1       450     450     3       29      2329315
Mean:   1       421.332 421.332 0.838009        4.9679
# of Seqs:      2329315
```
**number of Seqs:2329315**

- 04_unique.seqs

```
unique.seqs(fasta=cultivar.trim.contigs.good.fasta)
```
This reduced seqs from 2329315 to 1665541

- 05_count.seqs.

Error: when doing count, your group file contains more than 1 sequence name M04398_37_000000000-AVWAN_1_1116_17500_4330. Go back to the contigs.fasta file, report and name file, found that ForCV5_05 and Fresh_Ag_05 got same seqID.And this sequence is barcoded with ACT, should belong to ForCV5_05. So I go ahead and deleted corresponding seqs from corresponding files.

``count.seqs(name=cultivar.trim.contigs.good.names,group=cultivar.contigs.good.groups,processors=8)``

During this count process, it generate a count table with column being the count of each unique seq in each sample and rows beling unique sequence ID

- custom 16S reference-- start and end location determination

Refer to [Git hub] (https://github.com/mothur/mothur/issues/235) 
In my case, the start location is 6428, and end location is 23444

- pcr.seqs
```
pcr.seqs(fasta=silva.bacteria.fasta, start=6428, end=23444, keepdots=F, processors=8)
```
After the pcr.seqs process, it generated a file named silva.bacteria.pcr.fasta and this is my customized reference. To make this easy to remember and make more sense, Ichanged the file name to silva_V3_V4.fasta

- 06_align.seqs

 ```
 align.seqs(fasta=cultivar.trim.contigs.good.unique.fasta,reference=silva_V3_V4.fasta,processors=8)
 ```
- 07_summary.seqs

```
summary.seqs(fasta=cultivar.trim.contigs.good.unique.align,count=cultivar.trim.contigs.good.count_table,processors=8)

Using 8 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       8       3       0       2       1
2.5%-tile:      2       17016   403     0       4       58233
25%-tile:       2       17016   405     0       4       582329
Median:         2       17016   428     0       5       1164658
75%-tile:       2       17016   428     1       6       1746987
97.5%-tile:     2       17016   429     3       6       2271083
Maximum:        6701    17016   449     3       29      2329315
Mean:   22.7809 17014.7 420.165 0.838006        4.96751
# of unique seqs:       1665541
total # of seqs:        2329315
```
As we could see from the above summary results, most of the alignments start at 2 and end at 17016. The total unique alignment  number is 1665541, with total alignment number being 2329315.

- 08_screen.seqs  
At this point, we want to keep sequences that all started at 2 and ended at 17016. Meantime, we set maximum polymer number to be 8 as another way of quality control. 

```
screen.seqs(fasta=cultivar.trim.contigs.good.unique.align,count=cultivar.trim.contigs.good.count_table,summary=cultivar.trim.contigs.good.unique.summary,start=2,end=17016,maxhomop=8,processors=8)

summary.seqs(fasta=cultivar.trim.contigs.good.unique.good.align,count=cultivar.trim.contigs.good.good.count_table,processors=8)

Using 8 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       17016   369     0       3       1
2.5%-tile:      2       17016   403     0       4       57192
25%-tile:       2       17016   407     0       4       571916
Median:         2       17016   428     0       5       1143831
75%-tile:       2       17016   428     1       6       1715746
97.5%-tile:     2       17016   429     3       6       2230469
Maximum:        2       17016   449     3       8       2287660
Mean:   1.99999 17016   420.418 0.838114        4.95788
# of unique seqs:       1628156
total # of seqs:        2287660
```
- 09_filter.seqs  

As we could see from the above summary, although the number of base pairs between 369 and 449, but the length of the whole length of alignment is 17015. This is caused by the format of our reference, there are lots of dots and dashes inside of each alignment. Below is a short example of the alignment.

```
>M04398_37_000000000-AVWAN_1_1101_8581_2023
.G---G-G------G-A-A--TA-TT--GG-A-C-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AA-T-G-G--GC-----------------------------------------------------------------------------------------------------GC-A----A-----------------------------------------------------------------------------------------------------G-C--C-T-G-A-T-C-C-A---GC-C-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------A--T-GCC-G-C-G-T------------------------------------------------------------------------------------------------------------------------------------------------------------G-A-G--T--GA------------------------------------T-G--A--A--G-G-CC---------------------------------------------------------------------------------------TT-AG--------------------------------------------------------------------------------------G-G-T-
```
In order to make the downward distance_based OTU clustering easier and faster, we want to remove those useless dashes column if they are vertically the same among all of the alignment. vertical=T is the parameter define this. trump=. means that in one column, if one of alignment has ".", this column will be remove.
Here is the command we use to do this filter process.

``filter.seqs(fasta=cultivar.trim.contigs.good.unique.good.align,vertical=T,trump=.,processors=8)``

**Here is a brief summary after filter process. As it shown, 15838 columns are removed**

Length of filtered alignment: 1178
Number of columns removed: 15838
Length of the original alignment: 17016
Number of sequences used to construct filter: 1628156

- 10_unique.seqs (after filterr may make few sequences non unique)
``
unique.seqs(fasta=cultivar.trim.contigs.good.unique.good.filter.fasta,count=cultivar.trim.contigs.good.good.count_table)
1628156 1625778
``
Then number of unique alignment decreased from 1628156 to 1625778.

- 11_pre.cluster
If we think about the amplicon sequencing technique, the most common consideration is PCR error. As we know, error could also introduced during sequencing process. To denoising those errors and make the subsequent OTU calling more accurate, we use pre.cluster to further denoise the data. Usually, we allow one base pair of ambiguity per 100bp, this means that two sequence will to clustered to one if they have 1bp of difference along 100bp length. In our case, the length of contigs are around 429, so we set this diffs=5. In fact, we tried using a stringent denoising set up using diffs=2, which caused severe problems during subsequent OTU clustering, with the generated distance file larger than 700GB.  

```
pre.cluster(fasta=cultivar.trim.contigs.good.unique.good.filter.unique.fasta,count=cultivar.trim.contigs.good.unique.good.filter.count_table,diffs=2,processors=8)
```
- 12_chimera.vsearch and remove.seqs

Another quality control strategy is to detect the chimera and remove them. Chimera is caused by false alignment due to sequence error or uncorrect alignment.

```
chimera.vsearch(fasta=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.fasta,count=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.count_table,dereplicate=t,processors=8)

remove.seqs(fasta=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.fasta,accnos=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)
```

- 13_classify.seqs

At this point, we used several strategy to control the sequence quality, including length, polymer #, chimera. But we still need another step to remove unwanted sequences that belong to Chloroplast, Mitochondria, unknown, Archaea and Eukaryota. So before we could do this, we first need to classify the sequences to taxon.
```
classify.seqs(fasta=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta,count=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table,reference=trainset9_032012.pds.fasta,taxonomy=trainset9_032012.pds.tax,cutoff=80,processors=8)
```
- 14_remove.lineage
Here we are ready to remove all of the unwanted taxon
```
remove.lineage(fasta=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta,count=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table,taxonomy=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy,taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
```
- 15_summary.tax

Make a summary of the taxon.
```
summary.tax(taxonomy=current,count=current)
Using cultivar.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table as input file for the count parameter.
Using cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy as input file for the taxonomy parameter.
```
Here is what the taxonomy.summary file looks like:
And Phylotype based OTU clustering just used this taxonoy information to cluster sequences into OTUs. Taxlevel ranges from 1 to 6, represent Kingdom, phylum, class,order,family and genus.
```
taxlevel  rankID         taxon                                     daughterlevels  total    AgCV1_01  AgCV1_02  AgCV1_03
0         0              Root                                      1               2053643  63212     11012     16321
1         0.1            Bacteria                                  23              2053643  63212     11012     16321
2         0.1.1          "Acidobacteria"                           22              166485   942       307       1452
3         0.1.1.1        "Acidobacteria"_unclassified              1               333      3         2         4
4         0.1.1.1.1      "Acidobacteria"_unclassified              1               333      3         2         4
5         0.1.1.1.1.1    "Acidobacteria"_unclassified              1               333      3         2         4
6         0.1.1.1.1.1.1  unclassified                              0               333      3         2         4
3         0.1.1.2        Acidobacteria_Gp1                         1               39741    149       1         21
4         0.1.1.2.1      Acidobacteria_Gp1_order_incertae_sedis    1               39741    149       1         21
5         0.1.1.2.1.1    Acidobacteria_Gp1_family_incertae_sedis   1               39741    149       1         21
6         0.1.1.2.1.1.1  Gp1                                       0               39741    149       1         21
3         0.1.1.3        Acidobacteria_Gp10                        1               1490     6         6         17
4         0.1.1.3.1      Acidobacteria_Gp10_order_incertae_sedis   1               1490     6         6         17
5         0.1.1.3.1.1    Acidobacteria_Gp10_family_incertae_sedis  1               1490     6         6         17
6         0.1.1.3.1.1.1  Gp10                                      0               1490     6         6         17
3         0.1.1.4        Acidobacteria_Gp11                        1               228      0         1         1              
```

**Finally, we are ready to cluster our sequences to OTUs**

In fact, we have three options to process this OTU clustering

1) Distance matrixed based (This is medium computing consuming. If the sequence well preprocessed in forward procedure, suggest using this because this will give far more information than phylotype based clustering)
2) Phylotype based (This is the least computing consuming )
3) Phylogenetic based (This is computing aggresive, if several samples, easy to handle; if more samples it will take forever)

In my case, I used both distance based OTU clustering and phylotype based clustering.
 
 
**1) Average neribour distance based OTU clustering** 



- cluster.split
```
cluster.split(fasta=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, splitmethod=classify, taxlevel=4, cutoff=0.20,processors=10)
```
This process will generate a list file: cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list
*If you got mothur version higher than 1.39.0, you will be able to use opticlust method. Otherwise, average-neighbor, nearest_nerghbor and furthest_neighbor method*

- make.shared

We will set the similarity cutoff at 97% , i.e., label=0.03

```
make.shared(list=culticar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list, count=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)
```
This command will generate the OTU table that we are going to use. This is the format:

```
label   Group   numOtus Otu000001       Otu000002       Otu000003       Otu000004      
0.03    AgCV1_01        250316  0       978     16450   0       0       0
0.03    AgCV1_02        250316  97      1812    37      0       0       0
0.03    AgCV1_03        250316  102     1817    0       0       0       0
0.03    AgCV1_04        250316  61      798     36      0       0       0
0.03    AgCV1_05        250316  1727    568     0       48      0       0
0.03    AgCV1_06        250316  0       330     62      0       0       1
0.03    AgCV1_07        250316  0       188     100     74      0       27
0.03    AgCV1_08        250316  200     198     76      0       0       10
0.03    AgCV1_09        250316  0       4625    0       0       758     431
0.03    AgCV1_10        250316  3447    3684    0       0       0       0
0.03    AgCV2_01        250316  630     380     253     297     0       43
0.03    AgCV2_02        250316  170     63      48      52      0       45
0.03    AgCV2_03        250316  144     98      72      0       0       49
0.03    AgCV2_04        250316  140     88      61      0       0       63
0.03    AgCV2_05        250316  6431    146     592     0       0       0
0.03    AgCV2_06        250316  436     2126    101     93      0       0
0.03    AgCV2_07        250316  0       818     22      0       61      0

```
- classify.otu

classify all of the OTUs to their taxonoly using classify.otu
```
classify.otu(list=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list, count=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, label=0.03)
```
Gernated taxonomy file look like this:

```
    OTU          Size    Taxonomy
Otu000001       178415  Bacteria(100);"Proteobacteria"(100);Gammaproteobacteria(100);Pseudomonadales(100);Pseudomonadaceae(100);Pseudomonas(100);
Otu000002       111061  Bacteria(100);"Proteobacteria"(100);Alphaproteobacteria(100);Rhizobiales(100);Rhizobiaceae(100);Rhizobium(94);
Otu000003       56934   Bacteria(100);"Proteobacteria"(100);Gammaproteobacteria(100);Pseudomonadales(100);Pseudomonadaceae(100);Pseudomonas(100);
Otu000004       41375   Bacteria(100);"Proteobacteria"(100);Gammaproteobacteria(100);"Enterobacteriales"(100);Enterobacteriaceae(100);
Otu000005       39473   Bacteria(100);"Proteobacteria"(100);Gammaproteobacteria(100);Pseudomonadales(100);Pseudomonadaceae(100);Pseudomonas(100);
Otu000006       39064   Bacteria(100);"Proteobacteria"(100);Betaproteobacteria(100);Burkholderiales(100);Burkholderiaceae(100);Burkholderia(100);
Otu000007       34317   Bacteria(100);"Proteobacteria"(100);Gammaproteobacteria(100);"Enterobacteriales"(100);Enterobacteriaceae(100);
Otu000008       30362   Bacteria(100);"Proteobacteria"(100);Gammaproteobacteria(100);Pseudomonadales(100);Pseudomonadaceae(100);Pseudomonas(100);
Otu000009       30118   Bacteria(100);"Proteobacteria"(100);Gammaproteobacteria(100);Xanthomonadales(100);Xanthomonadaceae(100);Stenotrophomonas(100);
Otu000010       26774   Bacteria(100);"Proteobacteria"(100);Gammaproteobacteria(100);Pseudomonadales(100);Pseudomonadaceae(100);Pseudomonas(100);
Otu000011       25575   Bacteria(100);"Proteobacteria"(100);Gammaproteobacteria(100);Xanthomonadales(100);Xanthomonadaceae(100);Stenotrophomonas(100);
Otu000012       23072   Bacteria(100);"Proteobacteria"(100);Gammaproteobacteria(100);Xanthomonadales(100);Xanthomonadaceae(100);Rhodanobacter(100);
Otu000013       23013   Bacteria(100);"Proteobacteria"(100);Alphaproteobacteria(100);Sphingomonadales(100);Sphingomonadaceae(100);Novosphingobium(98);
Otu000014       21996   Bacteria(100);"Proteobacteria"(100);Gammaproteobacteria(100);"Enterobacteriales"(100);Enterobacteriaceae(100);
```


**2) Phylotype based OTU clustering**



- phylotype
``phylotype(taxonomy=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy)``

- make.shared
``make.shared(list=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.list,count=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table,label=1)``

- classify.otu
``classify.otu(list=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.list,count=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table,taxonomy=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy,label=1)``


### Summary


- **At the end of this whole pipeline, we got 250316 OTUs that belongs to 2053643 reads. The original reads number is 3,232,228. As a summary, after the whole pipeline, we got 63.5% reads left. AND the first screen.seqs is the main screen process, which removed about 902913 reads (27.9%)**

- **Using phylotype based OTU clustering we got about 661 OTUs, the total # of reads at the end of the pipeline is the same as distance-based OTU clustering**
***



## Stringent screening ##



** Changing ``maxambig`` parameter from ``3`` to ``0``**
** ``Pre.cluster`` still use ``diffs=5``**
#### At the end of the whole pipeline, we got **23908 OTUs** and  **total number of reads generated is 1033385 (32.06%)**





##                                          Microbial community analysis in R                                              ##



### File transfer

After OTU clustering, generated ``.shared`` and ``.con.taxonomy`` file were transfered to local computer using FileZilla. R software is used for subsequent statitic analysis. 





















