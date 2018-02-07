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
For this circular phylogenetic tree, I will use not only taxanomy information but also genus level relative abunance information for each treatment. About detailed R code, please find from this [link](https://github.com/liufangbaishikele/Soybean-rhizosphere-microbiome--16S-analysis/blob/master/cultivar_circular_tree_2nd.Rmd)
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
  * First, edit the taxanomy file to match with the tree file format required by GraPhlAn. Below is part of the tree file I built. **NOTE** Taxonomies are separated via ``.``.
  
  ```
Bacteria.Acidobacteria.Acidobacteria_Gp6.Acidobacteria_Gp6_order_incertae_sedis.Acidobacteria_Gp6_family_incertae_sedis.Gp6
Bacteria.Acidobacteria.Acidobacteria_Gp7.Acidobacteria_Gp7_order_incertae_sedis.Acidobacteria_Gp7_family_incertae_sedis.Gp7
Bacteria.Actinobacteria.Actinobacteria.Acidimicrobiales.Acidimicrobiaceae.Ilumatobacter
Bacteria.Actinobacteria.Actinobacteria.Acidimicrobiales.Acidimicrobiales_unclassified.Acidimicrobiales_unclassified
Bacteria.Actinobacteria.Actinobacteria.Acidimicrobiales.Acidimicrobineae_incertae_sedis.Aciditerrimonas
Bacteria.Actinobacteria.Actinobacteria.Acidimicrobiales.Iamiaceae.Iamia
Bacteria.Actinobacteria.Actinobacteria.Actinobacteria_unclassified.Actinobacteria_unclassified.Actinobacteria_unclassified
Bacteria.Actinobacteria.Actinobacteria.Actinomycetales.Actinomycetales_unclassified.Actinomycetales_unclassified
Bacteria.Actinobacteria.Actinobacteria.Actinomycetales.Actinospicaceae.Actinospica
Bacteria.Actinobacteria.Actinobacteria.Actinomycetales.Catenulisporaceae.Catenulispora
Bacteria.Actinobacteria.Actinobacteria.Actinomycetales.Geodermatophilaceae.Blastococcus
Bacteria.Actinobacteria.Actinobacteria.Actinomycetales.Glycomycetaceae.Glycomyces
Bacteria.Actinobacteria.Actinobacteria.Actinomycetales.Intrasporangiaceae.Intrasporangiaceae_unclassified
  ```
  * To avoid downward error during annotation, manually change:
   * ``Bacteria.Chloroflexi.Chloroflexi.Herpetosiphonales.Herpetosiphonaceae.Herpetosiphon`` to ``Bacteria.Chloroflexi.Chloroflexi_c.Herpetosiphonales.Herpetosiphonaceae.Herpetosiphon``
   * ``Bacteria.Nitrospira.Nitrospira.Nitrospirales.Nitrospiraceae.Nitrospira`` to ``Bacteria.Nitrospira_k.Nitrospira_c.Nitrospirales.Nitrospiraceae.Nitrospira_g``
   * ``Bacteria.Proteobacteria.Proteobacteria_unclassified.Proteobacteria_unclassified.Proteobacteria_unclassified.Proteobacteria_unclassified`` to ``Bacteria.Proteobacteria.Proteobacteria_unclassified_c.Proteobacteria_unclassified_o.Proteobacteria_unclassified_f.Proteobacteria_unclassified_g``
 * **NOTE** remember to change the corresponding taxa in node size file and otu_tax_combined_table

6. **Creat annotation file**
 * This is the most time_consuming part. Make sure during this process, always use **tab** between elements.
 
   i. title parameter
    
```
title   Rhizosphere microbiome between cultivars
```

   ii. ``total_plotted_degress`` - this is used to define how many degree of 360 will be used for plot tree. The left are used for ring label usage. In my case, I used 340 degree out of 360 degree
    
```
total_plotted_degrees   340
```
   iii. Global parameter set up
   * ``annotation_background_alpha`` were used to define the transparency of annotation background color. 
   * ``brantch_bracket_depth``  - is the ratio of bracket depth to the distance between two neigbour vertical nodes
   * ``annotation_font_size``
   * ``annotation_legend_font_size`` The will generate legend like the third column.
   * ``annotation_background_offset`` will determine the distance between the outside node and the first ring. Here I set it to be 0.1
   * ``title_font_size ``
   * ``brantch_bracket_width``  - ** still have not figure out.**
   * ``class_legend_font_size``  - ** not sure yet **
   * ``start_rotation `` - ** not sure yet **
   
   iV. Set up ``clade_marker_color``, ``clade_marker_size``  and ``clade_marker_shape``
   
   * Wildcard character could be used to define ``clade_marker_color`` and ``clade_marker_shape``. e.g., Acidobacteria* clade\_marker\_color #ff0000 will set up the color of all nodes that belongs to this phylum to red.
   * In terms of shapes, see the documentation from this [read.me](https://bitbucket.org/nsegata/graphlan/src/e91e79a421f96fdd28e8152b4de1c1b4e95ebb32/readme.txt?at=default&fileviewer=file-view-default)
   
   ** clade\_marker\_shape**
   ```
   Alphaproteobacteria*    clade_marker_shape	o
   Betaproteobacteria*     clade_marker_shape	*
   Deltaproteobacteria*    clade_marker_shape	p
   Gammaproteobacteria*    clade_marker_shape	s
   Proteobacteria_unclassified_c*  clade_marker_shape	D
   ```
  ** clade\_marker\_color**
  
  ```
  Acidobacteria*  clade_marker_color	#07aeba
  Actinobacteria* clade_marker_color	#b6cc0e
  Bacteroidetes*  clade_marker_color	#ed5567
  Chloroflexi*    clade_marker_color	#316022
  Firmicutes*     clade_marker_color	#f936f6
  Planctomycetes* clade_marker_color	#f2701a
  Proteobacteria* clade_marker_color	#ffee32
  Verrucomicrobia*        clade_marker_color	#3a44ff
  TM7*    clade_marker_color	#23b5ff
  ```
 ** clade\_marker\_size is generate from the node size files.
 
 ```
 Acidobacteria   clade_marker_size	981.6777199
 Proteobacteria  clade_marker_size	2302.54174
 Bacteroidetes   clade_marker_size	468.1885592
 Planctomycetes  clade_marker_size	283.5515276
 Verrucomicrobia clade_marker_size	366.2989901
 Actinobacteria  clade_marker_size	614.0824273
 ```
  V. Set up annotation background color correponding to clade marker colors.
  
  Vi. Using a trick to set up legend for different color of each Phylum.
  * Because their are Acidobacteria\_Phylum information inside of the tree file. So, it will treat is as isolated clade information. The color of each phylum is corresponding to ``clade_marker_color`` of each phylum.
  
  ```
  Acidobacteria_Phylum    clade_marker_color	#07aeba
  Actinobacteria_Phylum   clade_marker_color	#b6cc0e
  ```
  Combined with below
  ```
  Acidobacteria_Phylum    clade_marker_size	300
  Actinobacteria_Phylum   clade_marker_size	300
  ```
  Vii. Now ready for adding ring parameters. YAY!!!!
  
  1. Global options
  
    * ``ring_label_font_size``
    * ``ring_internal_separator_thickness``
    * ``ring_width``
    * ``ring_separator_color``
    * ``ring_label``
    * ``ring_label_color``
    
    ** NOTE ** I have no idea why when I copy the following contents from excel to my annotation file, it always got decode problems.
    Here are the error information
    
    ```
    Classes not implemented for external annotations
Traceback (most recent call last):
  File "/lustre/medusa/fliu21/anaconda2/bin/graphlan_annotate.py", line 56, in <module>
    ctree.annotate( args['annot'], args['outtree'] if args['outtree'] else args['intree'] ) # ,
  File "/lustre/medusa/fliu21/anaconda2/bin/src/graphlan_lib.py", line 320, in annotate
    Phylo.write( self.tree, out_file, "phyloxml")
  File "/lustre/medusa/fliu21/anaconda2/lib/python2.7/site-packages/Bio/Phylo/_io.py", line 83, in write
    n = getattr(supported_formats[format], 'write')(trees, fp, **kwargs)
  File "/lustre/medusa/fliu21/anaconda2/lib/python2.7/site-packages/Bio/Phylo/PhyloXMLIO.py", line 131, in write
    return Writer(obj).write(file, encoding=encoding, indent=indent)
  File "/lustre/medusa/fliu21/anaconda2/lib/python2.7/site-packages/Bio/Phylo/PhyloXMLIO.py", line 675, in __init__
    self._tree = ElementTree.ElementTree(self.phyloxml(phyloxml))
  File "/lustre/medusa/fliu21/anaconda2/lib/python2.7/site-packages/Bio/Phylo/PhyloXMLIO.py", line 688, in phyloxml
    elem.append(self.phylogeny(tree))
  File "/lustre/medusa/fliu21/anaconda2/lib/python2.7/site-packages/Bio/Phylo/PhyloXMLIO.py", line 652, in wrapped
    elem.append(getattr(self, method)(item))
  File "/lustre/medusa/fliu21/anaconda2/lib/python2.7/site-packages/Bio/Phylo/PhyloXMLIO.py", line 654, in wrapped
    elem.text = _serialize(obj.value)
  File "/lustre/medusa/fliu21/anaconda2/lib/python2.7/site-packages/Bio/Phylo/PhyloXMLIO.py", line 627, in _serialize
    return unicode(value)
UnicodeDecodeError: 'ascii' codec can't decode byte 0xc2 in position 1: ordinal not in range(128)

    ```
    The problems is that when I copy content from excel, their are some format that unix or python can not recognize
    
    ```
    ring_internal_separator_thickness	1	2 
    ring_width	1	0.8
    ring_separator_color	1	#888888
    ring_label	1	Ag_CV1
    ring_label_color	1	#ce0404
    ```
    
    When I head this part using ``head -38 cultivar.annot_04 | tail -l | od -c ``
    
    I found that at the end of the first line, it came out two numbers(302 240) that is not what I pasted into. See below:
    
    ```
    r   i   n   g   _   i   n   t   e   r   n   a
    l   _   s   e   p   a   r   a   t   o   r   _   t   h   i   c
    k   n   e   s   s  \t   2  \t   2 302 240
    ```
    So, I checked and editted the pasted content. It should be fine now.
  2. clade specific parameter in this format ``[clade_name] ring_option ring_level(integers)  parameter``
  
    * ``ring_alpha``
    * ``ring_height``
    * ``ring_color``
  ** Here, the ring alpha value were calculated as average relative abundance within treatment and multiply 10** 
  
  ** Note ** When you forgot to type in all the clade specific paramters and you set up more external ring levesl, during annotation, it will give me a warning that ``Classes not implemented for external annotations``
  Make sure add all of the clade ring information.
  
  3. ** NOTE ** Graphlan are very sensitive to extra characters or spaces or format 
    * ``head -n file | od -c `` to check is there are windows format that after pasting from excel
    * Using ``Vim`` to highlight trailing spaces and tabs in txt file in linux.
      1. ``Vim cultivar.annot_04``
      2. ``:set hlsearch ``
      3. ``/\s\+$``
      4. ``:%s/\s\+$//``
  
  4. Add rings for agriculture and forest fresh & bulk soil samples
    * Inside of R, subset ``r_filter_family_phyloseq`` object to fresh and bulk soil samples.
    * Merge phyloseq based on Treat factor and creat ``bulkfreshsub_r_filter_family_phyloseq``. Sample_sums 
    
    ```
    Agriculture_Bulk Agriculture_Fresh       Forest_Bulk      Forest_Fresh 
               12                 5                12                 5 
    ```
    * Combine OTU and taxonomy table to one table and write out as csv format
    * Inside of excel, make sure otu# order match with that of rhizosphere one.
    * Calculate alpha (average relative abundance within treatment * 10)
    
  5. Ready to run ``graphlan_annot.py`` and ``graphlan.py``
  
  ```
  python graphlan_annotate.py --annot cultivar.annot_07 cultivar.tree cultivar_annot_07.tree 
  python graphlan.py cultivar_annot_07.tree cultivar_tree_07.png --dpi 300 --pad 1 --size 16
  ```
 
 
 
 
 
 
 
 
 
 
 
  
  
  
  
  
  
  
  
  
