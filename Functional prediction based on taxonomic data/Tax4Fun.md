## Tax4Fun used to do function and pathway annotation based on 16S sequencing data

---

Find the documentation [here](http://tax4fun.gobics.de/)

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
      
     2. Error in read.table(file = file, header = header, sep = sep, quote = quote, : duplicate 'row.names' are not allowed - This one also make sense for me. Because this Tax4Fun software mainly use the taxonomy information and integrating to the precalculated pathway or KEGG orthologs. So OTU table does not make much sense. It will be good use genus level .shared and taxonomy file after rarefaction and remove rare genus with genus sum smaller than 50.
      
 3. Use the genus_based .shared file and taxonomy file
 
  * More details comes in
  

