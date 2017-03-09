# Soybean rhizosphere microbiome-soil type and cultivar impracts
Preliminary analysis for 16S seq data from Ag&amp;Forest cultivar project 


### sequence decompress

```
cd /lustre/projects/staton/projects/soybean_strigolactones/16S_raw_fastq/soybean_strigolactone_16S_03_02_2017
for file in *; do echo $file; cd $file; cp * ../../fastq_gunzip/; cd ..; done

```
### demultiplex each of the 96 fastq file based on TGA and ACT barcode

**Library prep protocol refer to** 

>Lundberg, D.S., Yourstone, S., Mieczkowski, P., Jones, C.D. and Dangl, J.L., 2013. Practical innovations for high-throughput amplicon sequencing. Nature methods, 10(10), pp.999-1002.

Barcode1------>F-Bc1_2,4,6 seqeunce is *TGA....TCACTCCTACGGG.GGC.GCAG*

Barcode2------>F-Bc2_1,3,5 sequence is *ACT....TCACTCCTACGGG.GGC.GCAG*
