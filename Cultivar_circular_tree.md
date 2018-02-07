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
  
4. For this circular phylogenetic tree, I will use not only taxanomy information but also genus level relative abunance information for each treatment.
  * First, prep-process phylotype-based shared file and cons taxonomy file ``cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.shared`` and ``cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.cons.taxonomy``
  * As sub.sample command in mothur could not combine subseted shared file with cons taxonomy file. So the below preprocessing will be done in R.
  * Remove genus that has tax_sum smaller than 100
  * Rarefy shared file to minimum depth
  * Write out edited OTU table and taxonomy table for generating circular tree. 
