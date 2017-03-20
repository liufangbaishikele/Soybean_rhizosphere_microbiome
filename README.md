# Soybean rhizosphere microbiome-soil type and cultivar impracts
  First run data

### sequence decompress

```
cd /lustre/projects/staton/projects/soybean_strigolactones/16S_raw_fastq/soybean_strigolactone_16S_03_02_2017
for file in *; do echo $file; cd $file; cp * ../../fastq_gunzip/; cd ..; done
cd /lustre/projects/staton/projects/soybean_strigolactones/16S_raw_fastq/fastq_gunzip
gunzip *
```
### demultiplex 

**Library prep protocol refer to** 

>Lundberg, D.S., Yourstone, S., Mieczkowski, P., Jones, C.D. and Dangl, J.L., 2013. Practical innovations for high-throughput amplicon sequencing. Nature methods, 10(10), pp.999-1002.

**Barcode1**
>F-Bc1_Fs2,4,6 seqeunce is *TGA....TCACTCCTACGGG.GGC.GCAG*

**Barcode2**
>F-Bc2_Fs1,3,5 sequence is *ACT....TCACTCCTACGGG.GGC.GCAG*

**demultiplex script**

**Read1 demultiplex**


1) Read1 barcoded with ACT
```
# "$1" is read1.fastq, "$2" is sampleID that barcoded with TGA, "$3" is sampleID that barcoded with ACT, "$4" read2.fastq
grep -B1 -A2 "ACT....TCACTCCTACGGG.GGC.GCAG" "$1" | grep -v "^--$" > "$3"_R1.fastq
grep "^@M04398" "$3"_R1.fastq > "$3"_header1.txt
awk '{print$1}' "$3"_header1.txt > "$3"_header2.txt
```
2) Read1 barcoded with TGA

```
# "$1" is read1.fastq, "$2" is sampleID that barcoded with TGA, "$3" is sampleID that barcoded with ACT, "$4" read2.fastq
grep -B1 -A2 "TGA....TCACTCCTACGGG.GGC.GCAG" "$1" | grep -v "^--$" > "$2"_R1.fastq
grep "^@M04398" "$2"_R1.fastq > "$2"_header1.txt
awk '{print$1}' "$2"_header1.txt > "$2"_header2.txt

```


*Now read files are ready for subsequent analysis using mothur*

## Mothur sequence processing ##

Sequence preprocess refer to [MiSeq Sop] (https://www.mothur.org/wiki/MiSeq_SOP)

### Modification include ###
- custom 16S reference region, start and end location determination
  Start location is 6428 end location is 23444

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
number of Seqs:      **3226258**


- screen.seqs ------







