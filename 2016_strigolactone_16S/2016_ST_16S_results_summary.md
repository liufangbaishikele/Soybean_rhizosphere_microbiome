
##                                     2016 Strigolactone - 16S rRNA results
  
 
---
 This results are based on 2016 Strigolactone over-expression 16S rRNA data.
---

**Before OTU clustering**

1. Raw sequences before doing assembly- 2985235
2. After first step screening (remove contigs with ambiguous calls and those too short or too long) - 2267164
3. After alignment and screening - 2218723
4. After remove chimeara sequences - 1818421
5. After remove non-bacteria lineage - 1810039 

* In total, got 43163 OTUs left after the whole pipeline across 36 samples


**Rarefaction and singleton discarding**

* To reduce the sequencing depth bias between samples, all samples were rarefied to minimum sequencing depth - 26187 reads per sample
* After this rarefaction process, 30281 OTU left.
* As singleton are tend to be caused by sequencing error, we will remove all singletons in R. After remove all the singletons, got 13255 OTUs left across 36 samples.

**Phylotype based process in mothur in order to link to LefSe analysis**

* After subsample of the taxonomy, count_table and list files
* Subsample_taxonomy will be used to generate the list file using phylotype command
* This updated subsample.tx.list will be used to do remove.rare command
* Then the subsample.tx.pick.list will be used to generate the .shared and .constaxonomy file


### 16S mothur output summary - Based on 48 samples

1. 59196 OTUs across 48 sample

2. Sample Max2_10 was moved from the phyloseq object due to its low sequencing depth

3. Remove all singleton before downward analysis

4. Rarefaction at minimum sequencing depth across samples.  After rarefaction, got 21430 OTUs left across 47 samples.

