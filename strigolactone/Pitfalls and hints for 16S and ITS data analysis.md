## Pitfalls and hints for 16S and ITS data analysis

**16S_and_ITS**

* Always download the most latest reference database for alignment and classification

**16S only**

* After make contigs, check the summary results and decide for the screen parameters.
* For classification, it has rdp and pds training set, with the latter training set has 119 mitochondrial 16S rRNA gene sequences as members of the Rickettsiales and four 18S rRNA gene sequences as members of the Eukarya. For Rhizosphere and bulk soil samples, it did not influence much. For seed and endo, I also tested the difference of output. It turned out, rpd performed better than pds.
* rarefy list file and generate shared file and taxonomy file OR do rarefaction in R. **Need confirmation in R**
* If no big difference, always rarefy list file and then generate downward shared, taxonomy file.

**ITS only**

* AS the length of ITS2 region varies a lot across different species. It is not necessary to customize the reference, neither alignment.
* For UNITE reference, 1) it has both singleton included or singleton excluded version. For singleton excluded version, it has 22265 taxa, while 50494 taxa included when singleton were included. 2) It has developer version and regular version. They differ in the length of fasta data. Developer version has longer length of corresponding sequences.
* Always check the taxonomy output of classify.seqs. Summarize all of the non-fungi taxa by runing

```
awk '{print $2}' .taxonomy | awk -F ';' '{print $1}' | sort | uniq
```
* Based on above results, using remove.lineage to get rid of all non-fungi sequences.

**ITS sequence analysis MOTHUR vs QIIME2**

**QIIME2** are too aggressive when using data2 based screening and error correction. 

```
Bulk_1     Control_1     D14_1     Max1_1    Max2_1 
10539      8506          6901      7864      7153 
```
Fungi_assigned 41%


**MOTHUR**

```
Bulk_1     Control_1     D14_1     Max1_1    Max2_1 
41428      27620          22538    27592     26378 
```

Fungi_unclassified 30.48%



