## Tax4Fun used to do function and pathway annotation based on 16S sequencing data

---

Find the documentation [here](http://tax4fun.gobics.de/)

Download the [Readme_Tax4Fun.pdf](http://tax4fun.gobics.de/RPackage/Readme_Tax4Fun.pdf) file for detailed procesure for this analysis. 

---

* **Important notes before using Tax4Fun**

1. It accept SILVang output, QIIMEbiomoutput using SILVA as a reference either in biom or txt file using below commands:
  (if you have the below formate output, then life will be easier for you. If not, no worry, you can find the example input format by downloading the SILVang and QIIME output file and format your own in a good shape )
  
   * ``importQIIMEBiomData``
   * ``importQIIMEData``
   * ``importSilvaNgsData``
    
2. Here are my notes during format my mothur OTU and taxonomy output into a good shape like SILVang output format.
    
   * When I combined the otu table and taxonomy file into one and importSilvaNgsData it poped out mainly two types of error
    
     1. Error in data[[rlabp]] : subscript out of bounds -- this may caused by the wrong format, substitute all comma with tab in ACF.
      
     2. Error in read.table(file = file, header = header, sep = sep, quote = quote) : duplicate 'row.names' are not allowed - This one also make sense for me. Because this Tax4Fun software mainly use the taxonomy information and integrating to the precalculated pathway or KEGG orthologs. So OTU table does not make much sense. It will be good use genus level .shared and taxonomy file after rarefaction and remove rare genus with genus sum smaller than 50.
      
 3. Use the genus_based .shared file and taxonomy file, which could be generated using phylotype based clustering. Label=1 means at genus level.
 
 ```
  Phylotype(taxonomy=##.taxonomy) 
  make.shared=(list=##.list,count=##.count_table,label=1)
  classify.otu(list=##.list,count=##.count_table,taxonomy=##.taxonomy,label=1)
 ```
 4. Optional- if you want to do a rarefaction based on a specific sequencing depth and the same time want to keep otu table and taxonomy table consistent (corresponding otu labels). ``sub.sample`` command could be used using list and count table as well as taxonomy table.
 
 ```
sub.sample(list=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.pick.tx.list,taxonomy=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.pick.taxonomy,count=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table,label=1,size=21776)

make.shared(list=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.pick.tx.1.subsample.list,count=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.subsample.count_table,label=1)

classify.otu(list=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.pick.tx.1.subsample.list,count=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.subsample.count_table,taxonomy=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.pick.subsample.taxonomy,label=1)
 ```
5. Prepare input file by combing OTU table with taxonomy information in excel into below format for Tax4Fun, and save it in ``.txt`` format.

```
AR3	BZ1	CL1	DF1	EB020	
0	  0	  0	  1	  0	Bacteria;Acidobacteria;Acidobacteria;Acidobacteriales;Acidobacteriaceae (Subgroup 1);uncultured;
0 	0	  0	  2	  0	Bacteria;Acidobacteria;Acidobacteria;Subgroup 2;
0	  0	  0	  0	  2	Bacteria;Acidobacteria;Acidobacteria;Subgroup 4;Family Incertae Sedis;Blastocatella;
0	  0	  2	  1	  0	Bacteria;Acidobacteria;Acidobacteria;Subgroup 6;
0	  0	  0	  1	  0	Bacteria;Acidobacteria;Holophagae;Subgroup 7;
1	  0	  0	  0	  0	Bacteria;Actinobacteria;Actinobacteria;Corynebacteriales;Mycobacteriaceae;Mycobacterium;
```
6. Save the tab delimited ``.txt`` file into local computer and rename the file name to ``.csv``. I do not know why, but this is the way it works for the command of ``importSilvaNgsData``

**Download required [SILVA Reference data](http://tax4fun.gobics.de/Tax4Fun/ReferenceData/SILVA115.zip) ,no matter which version and unzip the file and put it into the same folder you put your input data**

**Start analysis inside of R **

```
OTU<-importSilvaNgsData(inputFiles = "Combined_otu_and_tax_for_Tax4Fun.csv")
str(OTU)
colnames(OTU$otuTable)<-OTU$sampleNames
function_profile<-Tax4Fun(OTU,"SILVA115",fctProfiling = TRUE,refProfile ="UProC",normCopyNo = TRUE)
str(function_profile)
```



