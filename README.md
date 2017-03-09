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

**Barcode1**
>F-Bc1_Fs2,4,6 seqeunce is *TGA....TCACTCCTACGGG.GGC.GCAG*

**Barcode2**
>F-Bc2_Fs1,3,5 sequence is *ACT....TCACTCCTACGGG.GGC.GCAG*

**demultiplex script**

"$1" is read1.fastq, "$2" is sampleID that barcoded with TGA, "$3" is sampleID that barcoded with ACT, "$4" read2.fastq

```

grep -B1 -A2 "ACT....TCACTCCTACGGG.GGC.GCAG" "$1" | grep -v "^--$" > "$3"_R1.fastq
grep "^@M04398" "$3"_R1.fastq > "$3"_header1.txt
awk '{print$1}' "$3"_header1.txt > "$3"_header2.txt

for line in $(cat "$3"_header1.txt)
do
       echo $line | grep "^@" >> "$3"_header2.txt
done


for header in $(cat "$3"_header2.txt)
do
        grep -A3 "$header" "$4" >> "$3"_R2.fastq

done

grep -B1 -A2 "TGA....TCACTCCTACGGG.GGC.GCAG" "$1" | grep -v "^--$" > "$2"_R1.fastq
grep "^@M04398" "$2"_R1.fastq > "$2"_header1.txt
awk '{print$1}' "$2"_header1.txt > "$2"_header2.txt

for line in $(cat "$2"_header1.txt)
do
       echo $line | grep "^@" >> "$2"_header2.txt
done

for header in $(cat "$2"_header2.txt)
do
        grep -A3 "$header" "$4" >> "$2"_R2.fastq
done

```






