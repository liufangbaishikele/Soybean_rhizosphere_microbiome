#                    Cultivar 2016 project graphlan and igraph practice

-------
## GraPhlAn analysis practice

[Refered to: bitbucket/Nicola Segata repository/ GraPhlAn wiki page](https://bitbucket.org/nsegata/graphlan/wiki/browse/)

**Practice environment**: beacon server of the National Institute for Computational Sciences (NICS)

#### Raw data 

* OTU count table (generated follow mothur Miseq SOP pipeline) *cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.shared*
* Taxonomy file *cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.cons.taxonomy*
* Meta data 
```
Sample_names	SampleID	Soil_type	    Treatment	Compartment
Ag_B_1	      Ag_B_01	  Agriculture	  Bulk	      Bulk
Ag_B_10	      Ag_B_10	  Agriculture	  Bulk	      Bulk
Ag_B_11	      Ag_B_11	  Agriculture	  Bulk	      Bulk
Ag_B_12	      Ag_B_12	  Agriculture	  Bulk	      Bulk
Ag_B_2	      Ag_B_02	  Agriculture	  Bulk	      Bulk
Ag_B_3	      Ag_B_03	  Agriculture	  Bulk	      Bulk
Ag_B_4	      Ag_B_04	  Agriculture	  Bulk	      Bulk
Ag_B_5	      Ag_B_05	  Agriculture	  Bulk	      Bulk
Ag_B_6	      Ag_B_06	  Agriculture	  Bulk	      Bulk
```
### Prapare .biom format input for export2graphlan
* log into beacon using ssh NETID@duo.acf.tennessee.edu
* module load mothur
* Go to mothur environment and using make.biom command to generate biom file as input for export2graphlan
```
make.biom(shared=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.shared,constaxonomy=cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.cons.taxonomy,matrixtype=dense,label=0.03,metadata=cultivar_meta.txt)
```
Tips: Here if I want to use picrust parameter(used for bacterial community function prediction), several things is needed for this work. 
``1`` I need to set reftaxonomy parameter, which is the greengene version taxonomy for class.seqs process. Could be downloaded from ([greengenes reference taxonomy link](http://www.mothur.org/w/images/6/68/Gg_13_8_99.taxonomy.tgz)) 
``2`` At the same time, the constaxonomy file need to be classified using greengene reference taxonomy.
``3`` Picrust paramter needs green genes OTU IDs map table together with reference taxonomy , which could be downloaded from [this link](http://www.mothur.org/w/images/b/be/GG_13_5_otuMapTable.zip)

In my case, I will not do function prediction, so just skip using picrust paramter.

### Using export2graphlan to produce annocation and tree file for GraPhlAn

[export2graphlan](https://bitbucket.org/CibioCM/export2graphlan) This is a conversion software tool for producing both annotation and tree file for GraPhlAn.

1) Install export2graphlan
* First from **Anaconda**- Installation form anaconda cloud ``conda install -c bioconda export2graphlan``. It always conflict with other packages. 
* Then from **bitbucket**, I tried to follow [this documentation](https://bitbucket.org/CibioCM/export2graphlan) 
by using ``hg clone https://hg@bitbucket.org/CibioCM/export2graphlan``
BUT this ``hg`` command is a Mercurial command, which does not work on beacon.

---

**What is Mercurial?**

It is a DVCS that transfers code between your local system and Bitbucket Cloud.

---

* OK, no problem. Let's intall Mercurial following this documentation from [Bitbucket support](https://confluence.atlassian.com/get-started-with-bitbucket/mercurial-setup-860009660.html) , and then we can use ``hg`` command to install export2graphlan from bitbucket cloud. Again, got stucked at the first step ``cat /etc/apt/sources.list``. All right, I will just go to Mercurial website to download the source and install on beacon. It disappointed me because it does not have linux source. No problem, I can install it on my computer (Mac), however, Mercurial is written in Python (which means I also need to install pathon on my computer, and maybe other software that is needed for download pipeline)
* Finally I give up with direct intallation of export2graphlan

**HERE are my strategies**:

1) I explored in their [Bitbucket depository](https://bitbucket.org/CibioCM/export2graphlan/src/db0a809958d7ed860da44c9f9d51f2c9b068757f?at=default) AND looked at their source code. Found that this funciton is build using python. Great!

* Python is available on beacon. AND, export2graphlan requires the following additional library:

```
pandas ver. 0.13.1 (pandas)
BIOM ver. 2.0.1 (biom-format, only if you have input files in BIOM format) # installed by conda
SciPy (scipy, required by hclust2)
```
Then I added those library to my environment using conda

* When I run export2graphlan.py using its example data from export2graphlan ``/ examples / hmp_aerobiosis / `` directory, unfortunately, I got another error indicate that hclust2.hclust2 is no loaded.

* So, I looked into the source directory in hclust2 repository ``Nicola Segata/hclust2``. Then downloaded the [hclust2 directory](https://bitbucket.org/nsegata/hclust2/get/3d589ab2cb68.zip). By the way, this hclust2 is scripted with python language. 

## Here 








