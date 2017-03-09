# Soybean rhizosphere microbiome-soil type and cultivar impracts
Preliminary analysis for 16S seq data from Ag&amp;Forest cultivar project 
## sequence decompress

```
cd /lustre/projects/staton/projects/soybean_strigolactones/16S_raw_fastq/soybean_strigolactone_16S_03_02_2017
for file in *; do echo $file; cd $file; cp * ../../fastq_gunzip/; cd ..; done

```
## demultiplex each of the 96 fastq file based on TGA and ACT barcode

>F-Bc1_2,4,6 seqeunce is *TGA....TCACTCCTACGGG.GGC.GCAG*

>F-Bc2_1,3,5 sequence is *ACT....TCACTCCTACGGG.GGC.GCAG*
