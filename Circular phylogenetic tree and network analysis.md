#                    Cultivar 2016 project graphlan and igraph practice

-------
GraPhlAn analysis practice

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

[export2graphlan](https://bitbucket.org/CibioCM/export2graphlan) This is a conversion software tool for producing both annotation and tree file for GraPhlAn.

1) Install export2graphlan
* First try - Installation form anaconda cloud ``conda install -c bioconda export2graphlan``. It always conflict with other packages. 
* Then, I tried to follow [this documentation](https://bitbucket.org/CibioCM/export2graphlan) 
by using ``hg clone https://hg@bitbucket.org/CibioCM/export2graphlan``
BUT this ``hg`` command is a Mercurial command, which does not work on beacon.

---

**What is Mercurial?**


Well. From this [Bitbucket support](https://confluence.atlassian.com/get-started-with-bitbucket/mercurial-setup-860009660.html) webpage I found that it is a DVCS that transfers code between your local system and Bitbucket Cloud.

---

* OK, no problem. Let's intall Mercurial first, then using hg command to install export2graphlan from bitbucket cloud. Again, stucked at the first step ``cat /etc/apt/sources.list``






