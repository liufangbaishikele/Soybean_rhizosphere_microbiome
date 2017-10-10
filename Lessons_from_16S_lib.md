# Lessons from soybean 16S sequencing - protocol wise and strategy wise

## Protocol wise

The deeper I play with my data, the more flaws and problems I realized.

- Right now, I had sequenced my library at low clustering density due to several fators including, 

1) incorrect quantification of libray 

2) Potential primer competition because we load custom primer to there primer mix well on kit cassete 

3) Got more Phix reads (20%), which expected to be around 10% as we add 60uL 8pM Phix and 540uL "8pM" library. 

``
Suddenly, I got an idea to evaluate quantification of library based on the expected Phix read and output Phix read. It seems the concentration of our library is half of that of Phix (8pM).
``
4) Read2 quanlity did not pass QC analysis. 

    *Does this because poor library prep or sequencing setting problem?? How could I know??*
    
    ** Mock community is very important from evaluating library prep performance and subsequent sequencing process to sequence process**
    
    *Is there any way to evaluate sequencing performance?*
    Considerations of Miseq setup:
       * optimal loading concentration: too high, will decrease sequencing quality; too low, results low cluster density, few read yied 
       
## Mothur analysis wise

- Read preprocess and file rearrangement

1) keep a good file name habbit, using sample_01 instead of sample_1
2) Considering the subsequent sample_ID extraction, name sampleID in a wise way -- refer to mothur analysis documentation
3) After make.contigs, found that ambigous base pair is higher than MiSeq SOP. Due to lack of overal overlap between R1 and R2, AND lower sequencing quanlity at each end of sequencing.
4) Unique.seqs: this process could only reduce the contigs from 2# to #, whereas, populus dataset reduce from 7# to #.
5) This problem is severe when I process cluster.split process because there tons of unique sequence used to calculate distance matrix before real clustering process.It yields a distance file.3.dist being 700GB but still did not finish.

-  Realized the potential problem of my data, I went back checked the Read2 fastQC report of mine (before cutadapt2 trimming and after), Xiaolong's, populus projects. Found that my fastQC report is pretty bad and xiaolong's too. From my understanding, both of our sequences are V3_V4 amplicon sequences using MiSeq kitV3 300x2 cycle cassete. And soil microbiome, its self will be far complicated than that of gut microbiome.  But his sample number is small, only 15, but mine is 139+48. I asked Miriam about Rcorrector, although it is an RNA-seq corrector. Unsuperisingly, no improvement could be see from fastQC report.
``perl /lustre/projects/staton/software/rcorrector/run_rcorrector.pl -s trimmed_Ag_B_01_R2.fastq -k 31 
``
- From my understanding of Schloss mothur forum is: lack of complete overlap between read1 and read2, especially when we got low sequencing quanlity after 100bp. It will be a big problem when making contigs, because ambiguous base will be introduced because we are not sure which read should we trust if base calling is different between read1 and read2.

-  Current situation, I could not repeat the library prep process and together with mock community as strong standard for procedure evaluation. What I can do is: Go back to try precluster procedure and set diffs=5. And go through the process downward to OTU cluster to see what I will get.

- If modification of mothur pipeline does not work, I will definitely give up on this dataset. Repeat the whole experiment, with GR24, nor MS, nor germination evaluation just RT-qPCR expression level confirmation.

1) For the libray prep, IF I do by myself, NO melecule tags, NO frameshift as they are useless if we did no go with 1 cycle during molecular tagging PCR process. 13 treatment X 7 replicates=91 together with mock community and negative control. Dual indexing is not neccessary. But PNA blocking is still good to go.

2) IF Meg agree that we send our samples out for sequencing, that will be great too. It can save me buch of time. I could learn some cool stuff and prep article.
