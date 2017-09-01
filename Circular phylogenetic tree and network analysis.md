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

**Using export2graphlan to produce annocation.txt and tree.txt file for GraPhlAn**

[export2graphlan](https://bitbucket.org/CibioCM/export2graphlan) This is a conversion software tool for producing both annotation and tree file for GraPhlAn.

1) Install export2graphlan
* First from **Anaconda**- Installation form anaconda cloud ``conda install -c bioconda export2graphlan``. It always conflict with other packages. 

*Question: when I install anaconda, it asked me if I wish the installer to prepend the anacaonda2 install location to PATH in /nics/d/home/fliu21/.bash.rc

I chose no as I am running out of space in my home directory, so I want this PATH being /lustre/medusa/fliu21/anaconda2/bin
As a solution, everytime I want to intall any software from anaconda cloud, I first need to set up the environment using ``export PATH=/lustre/medusa/fliu21/anaconda2/bin:$PATH``


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

**Files inside of the directory of /lustre/medusa/fliu21/circular_phylogenetic_tree** 

```
-rw-r--r--. 1 fliu21 tug2004  89M Aug 30 18:23 cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.biom
-rwxr-xr-x. 1 fliu21 tug2004  29K Aug 30 16:23 export2graphlan.py
-rw-r--r--. 1 fliu21 tug2004 661K Aug 30 16:25 lefse_input.txt
-rw-r--r--. 1 fliu21 tug2004 101K Aug 30 16:26 lefse_output.txt
-rw-r--r--. 1 fliu21 tug2004  575 Aug 30 16:27 PIPELINE.sh
-rw-r--r--. 1 fliu21 tug2004 5.3K Aug 30 17:00 LICENSE.txt
drwxr-xr-x. 3 fliu21 tug2004 4.0K Aug 30 16:40 hclust2

```
The hclust2 directory has below files inside

```
drwxr-xr-x. 3 fliu21 tug2004 4.0K Aug 30 16:39 examples
-rwxr-xr-x. 1 fliu21 tug2004  35K Aug 30 16:39 hclust2.py
-rw-r--r--. 1 fliu21 tug2004  34K Aug 30 16:40 hclust2.pyc
-rwxr-xr-x. 1 fliu21 tug2004    0 Aug 30 16:39 __init__.py
-rw-r--r--. 1 fliu21 tug2004  138 Aug 30 16:40 __init__.pyc
-rw-r--r--. 1 fliu21 tug2004 9.5K Aug 30 16:39 README.md
```
Now preparation is done for runing export2graphlan, Let's do it!

```
module load python

python export2graphlan.py -i cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.biom  -t tree.txt -a annot.txt --title "HMP Aerobiosis" --annotations 2,3 --external_annotations 4,5,6 --fname_row 0 --skip_rows 1,2 --ftop 200 
```
After a while, it will generate the tree.txt and annotation.txt file

----

**Yay, ready to run GraPhlAn software and build circular tree

First, loading GraPhlAn to my environment using ``conda install -c bioconda graphlan  ``
Then navigate to /lustre/medusa/fliu21/circular_phylogenetic_tree ``cd /lustre/medusa/fliu21/circular_phylogenetic_tree``
Run graphlan
```
graphlan_annotate.py tree.txt tree_annot.txt --annot annot.txt
graphlan.py tree_annot.txt cultivar_tree.png --dpi 150 --size 14 

```
After several seconds I generated a circular tree. I can not wait to transfer to my local computer and look it. [Here](https://drive.google.com/drive/my-drive) it is. It is very wired looking no rings and useful information. I think the problem comes from the annotation from biom file. It can not be annotated properly, so only a tree is built but no abundance and metadata information.

## Now I will switch to LEfSe as this is what they used for export2graphlan demo.

1) First, use make.lefse command (inside of mothur) to create lefse formatted file. The LEfSe formatted file is can be used as an input to the LEfSe program. 

2) install lefse using conda 
* First need to set up channel following this [link](https://bioconda.github.io/recipes/lefse/README.html)
* Using  ``conda install lefse`` to finish installation

3) cp format\_input.py to my working directory (because when I install lefse, this function is not listed inside my anaconda2/bin directory)
```
python format_input.py cultivar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.lefse cultivar_lefse.in -c 1 -s 2 -u 3 -o 1000000
```
You know what? In fact, I found this funciton in my anaconda2/bin folder, but the function name is ``lefse-format_input.py``
Yeah, sometimes, the documentation may not updated with different version of their software

4) After formating the file in a way for lefse anaylysis. Here is a glance of the input data generated by ``format_input.py``

```
(dp0
S'class_sl'
p1
(dp2
S'Agriculture'
p3
(I0
I68
tp4
sS'Forest'
```
5) Now operate ``run_lefse.py``
* Again, I got stuck by the below error: 

        ```
        run_lefse.py
        Traceback (most recent call last):
          File "/lustre/medusa/fliu21/anaconda2/bin/run_lefse.py", line 4, in <module>
            from lefse import *
          File "/lustre/medusa/fliu21/anaconda2/share/lefse-1.0.8.post1-0/lefse.py", line 3, in <module>
            import rpy2.robjects as robjects
          File "/lustre/medusa/fliu21/anaconda2/lib/python2.7/site-packages/rpy2/robjects/__init__.py", line 15, in <module>
            import rpy2.rinterface as rinterface
          File "/lustre/medusa/fliu21/anaconda2/lib/python2.7/site-packages/rpy2/rinterface/__init__.py", line 99, in <module>
            from rpy2.rinterface._rinterface import *
        ImportError: libicuuc.so.58: cannot open shared object file: No such file or directory
        ```
* After google and try different solutions, finally solved by ``conda update --all``
* Here are the command I used for running ``run_lefse.py``
```
run_lefse.py cultivar_lefse.in cultivar_lefse.res
```
Got an error: 

```
/lustre/medusa/fliu21/anaconda2/lib/python2.7/site-packages/rpy2/rinterface/__init__.py:186: RRuntimeWarning: Error in kruskal.test.formula(y ~ x1 + x2 + x3, ) : 
  'formula' should be of the form response ~ group
```
* There maybe some format thing that is wrong when I generate the lefse.in file using my dataset. Let me try to run this command use their example data.


----




