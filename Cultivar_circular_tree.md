## Cultivar project - circular phylogenetic tree

1. **Install anaconda**

```
wget https://repo.continuum.io/archive/Anaconda3-4.4.0-Linux-x86_64.sh
bash Anaconda3-4.4.0-Linux-x86_64.sh
``` 
2. **Download GraPhlAn python-based packge from Bitbucket cloud**

``wget https://bitbucket.org/nsegata/graphlan/downloads/``
``unzip`` downloaded zip file using ``unzip file.zip
``
3. **Test graphlan.py**
  * After download and unzip installation package, navigate to the directory including graphlan.py
  * Test command by running ``python graphlan.py --help``. It poped out error information
  ```
  Traceback (most recent call last):
  File "graphlan.py", line 23, in <module>
    from src.graphlan_lib import CircTree as CTree
  File "/lustre/haven/gamma/staton/projects/soybean_rhizosphere/05_final_run/16S_cultivar_proj/circular_tree/cultivar_circular_tree_2nd/src/graphlan_lib.py", line 1, in <module>
    from Bio import Phylo
ImportError: No module named Bio
  ```
  * It indicate that biopython module is not installed. To check if this is the problem, run python command ``import Bio``, it told me that I do not have biopython module. So install biopython through anaconda by ``conda install -c anaconda biopython``
  * Now if I type ``python graphlan.py`` it just output information of usage and too few arguments.
  * Now installation are done and we are ready for build tree and produce annotation file.
  
4. **R data pre-process**
For this circular phylogenetic tree, I will use not only taxanomy information but also genus level relative abunance information for each treatment.
  * First, prep-process phylotype-based shared file and cons taxonomy file ``cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.shared`` and ``cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.cons.taxonomy``
  * As sub.sample command in mothur could not combine subseted shared file with cons taxonomy file. So the below preprocessing will be done in R.
  * Rarefy OTU table to minimum sequencing depth
  * Prune OTU table to keep top 200 genus based on tax_sums results
  * Transform count phyloseq object to relative abundance
  * Merge samples based on treatment category. Here fun=mean does not work for OTU table, in fact, it use fun=sum (means it sums value of replications that belongs to each treatment)
  ```
  Agriculture_Bulk, Agriculture_CV1, Agriculture_CV2, Agriculture_CV3, Agriculture_CV4, Agriculture_CV5, Agriculture_CV6, Agriculture_Fresh, Forest_Bulk, Forest_CV1, Forest_CV2, Forest_CV3, Forest_CV4, Forest_CV5, Forest_CV6, Forest_Fresh  
  ```
  * Combine taxa table and otu table and write out the dataset to local computer for manual annotation and tree building
  * Generate node size 
  ```
  Genus_node_size<-data.frame(name=tax_table(merge_r_graphlan_phyloseq)[,6],size=taxa_sums(merge_r_graphlan_phyloseq))
  Family...
  Order...
  Class...
  Phylum...
  ```
 5. Generate tree file and annotation file - Now the downward process are done in ACF
  * 
 
 
 
 
 
 
 
 
 
 
 
 
  
  
  
  
  
  
  
  
  
