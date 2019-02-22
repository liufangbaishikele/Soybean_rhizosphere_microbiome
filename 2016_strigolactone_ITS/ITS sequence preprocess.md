-----
ITS sequencing for fungi profiling
-----

### What is ITS? 

*ITS efers to the spacer DNA situated between the small-subunit ribosomal RNA (rRNA) and large-subunit rRNA genes in the chromosome*

### Primer sets

For ITS sequencing we targeted ITS2 region with customized primer sets composed of a mixture of 6 forward and 2 reverse primers designed to detect as much diverse taxa as possible- refer to this paper for detaield information (The Populus holobiont: dissecting the effects of plant niches and genotype on the microbiome)[https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0413-8].

### Mothur pipelines and trimming process

* Raw reads before any trimming - 2537994
* After the first step of trimming, all contigs with ambiguous callings as well as sequences that are too short or too long are removed - 2084687
* Chimera sequences that caused by the sequencing error or assembly error were detected and removed. After discarding all the chimeras, get 2069946 sequences left
* Based on the UNITE reference, sequences are classified to taxa. In our case, we have lots of unclassified sequences. So after remove all the non-fungi lineage (unknown-Metazoa-Plantae-Rhizaria), we got 1174908 sequences left.  
* At the end, needleman based alignment method were used to calculate the sequence distance. Then based on the distance, sequences were clustered into OTUs. So, all the sequences were clustered into 13204 OTUs.
* After analysis in R, I discovered that unclassified-fungi account for 30-50% of the sequences in my samples. Which is not common. After consulting with Dr. Mellisa in ORNL and discuss with my PI, I extracted all of the sequences that were classified to unclassified-fungi. Those extracted fungi-unclassified sequences were blasted against nematoda ITS reference. It turned out that most of the sequences were classified to two main nematode taxa (Pratylenchus and Acrobeloides) with e value equals to zero. Actually, for lots of the sequences could either be assigned to one of them with similar confidence (which indicate that ITS is not a very efficient tools for classification). 

```
remove.seqs(fasta=ST.trim.contigs.good.unique.precluster.pick.pick.fasta,taxonomy=ST.trim.contigs.good.unique.precluster.pick.UNITEv6_sh_97.wang.pick.taxonomy,count=ST.trim.contigs.good.unique.precluster.denovo.vsearch.pick.pick.count_table,accnos=nematoda_seq_ID.list)
```
* After remove-nematode sequences, got 971362 sequence across 36 samples, which belonged to 12087 OTUs.
* To eliminate the sequencing depth bias between samples, the shared file were rarefied to minimum sequencing depth across all samples (which is 1082) and all singleton were removed after this rarefaction. At the end, 385262 sequences belongs to 2970 OTUs were left across 36 samples.



 
