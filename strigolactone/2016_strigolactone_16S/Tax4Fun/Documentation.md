## Tax4Fun used to predict KEGG pathway based on 16S rRNA amplicon sequencing data

Refer to this [Tax4Fun documentation](http://tax4fun.gobics.de/)

* In terms of the [input file](https://github.com/liufangbaishikele/Soybean_rhizosphere_microbiome/blob/master/strigolactone/2016_strigolactone_16S/Tax4Fun/Tax4Fun_input.csv) and [R code](https://github.com/liufangbaishikele/Soybean_rhizosphere_microbiome/blob/master/strigolactone/2016_strigolactone_16S/Tax4Fun/2016_strigolactone_16S_Tax4Fun_silva.Rmd) 

* I found Tax4Fun very strict with input format. The key is to format the input as the same as they used in their example.

* For me, I used the shared file and con.taxonomy mothur output file based on phylotype cluster. Which I think is more efficient as OTU input will not provide any information in terms of linking taxonomy information to precalculated KEGG profile.

* **This is how I formatted my input file**

1. Generate Otu count table in R
  
  setwd("/Users/fangliu/Documents/2016_strigolactone_project/16S_after_remove_Agrobacterium/silva_classification/Tax4Fun/Tax4Fun")

ST_genus_phyloseq<-import_mothur(mothur_shared_file = "/Users/fangliu/Documents/2016_strigolactone_project/16S_after_remove_Agrobacterium/silva_classification/strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.pick.subsample.tx.1.pick.shared",mothur_constaxonomy_file = "/Users/fangliu/Documents/2016_strigolactone_project/16S_after_remove_Agrobacterium/silva_classification/strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.pick.subsample.tx.1.pick.1.cons.taxonomy") 
  `# This shared file were exported from '/staton/projects/soybean_rhizosphere/2016_strigolactone/16S_2016_strigolactone/Mothur_analysis/02_mothur/rarefied_files'`

colnames(tax_table(ST_genus_phyloseq))<-c("Kingdom","Phylum","Class","Order","Family","Genus")

  `The original list file was subseted in mothur using sub.sample command to do rarefaction-based normalization. And then the list file were filtered to remove all singletons.`

ST_genus_otu_table<-otu_table(ST_genus_phyloseq)
dim(ST_genus_otu_table)
ST_genus_otu_table[1:5,1:5]
colnames(ST_genus_otu_table)
write.csv(ST_genus_otu_table,file="subsampled_and_rare_removed_ST_genus_shared.csv")
  
2. Combine otu count table with taxonomy table- this part is done using shell command


      * Generated subsampled_and_rare_removed_ST_genus_shared.csv file inside of R as below

      ```
      "","Bulk_1","Bulk_2","Bulk_3","Bulk_4","Bulk_5","Bulk_6","Bulk_7"
      "Otu0001",2162,2058,2022,2083,1926,2121,1993
      "Otu0002",440,440,409,374,447,426,521
      "Otu0003",807,722,698,733,956,752,849
      ```

      * Remove all "" and the first column as these Otu list is not allowed in the final Tax4Fun input.

      ```
      sed 's/"//g' subsampled_and_rare_removed_ST_genus_shared.csv | cut -f2- -d"," > ST_otu_for_Tax4Fun
      ```
      After this, the ST_otu_for_Tax4Fun file will looks like below.

      ```
      Bulk_1,Bulk_2,Bulk_3,Bulk_4,Bulk_5,Bulk_6
      2162,2058,2022,2083,1926,2121
      440,440,409,374,447,426
      807,722,698,733,956,752
      978,936,833,995,1009,921
      ```


      * Modify taxonomy file in below format using this command ``awk '{print $3}' strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.pick.subsample.tx.1.pick.1.cons.taxonomy | sed 's/(100)//g' > ST_tax_for_Tax4Fun``


      ```
      Taxonomy
      Bacteria;Acidobacteria;Subgroup_6;Subgroup_6_or;Subgroup_6_fa;Subgroup_6_ge;
      Bacteria;Planctomycetes;Planctomycetacia;Gemmatales;Gemmataceae;uncultured;
      Bacteria;Verrucomicrobia;Verrucomicrobiae;Pedosphaerales;Pedosphaeraceae;Pedosphaeraceae_ge;
      Bacteria;Planctomycetes;Phycisphaerae;Tepidisphaerales;WD2101_soil_group;WD2101_soil_group_ge;
      Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Rhizobiaceae;Rhizobiaceae_unclassified;
      Bacteria;Patescibacteria;Saccharimonadia;Saccharimonadales;Saccharimonadales_fa;Saccharimonadales_ge;
      Bacteria;Bacteria_unclassified;Bacteria_unclassified;Bacteria_unclassified;Bacteria_unclassified;Bacteria_unclassified;
      Bacteria;Proteobacteria;Gammaproteobacteria;Betaproteobacteriales;Burkholderiaceae;Burkholderiaceae_unclassified;
      Bacteria;Proteobacteria;Deltaproteobacteria;Myxococcales;Haliangiaceae;Haliangium;
      ```

      * Combine OTU count table with taxonomy table and save in txt format. Then save as csv file.

      ```
      paste -d "," ST_otu_for_Tax4Fun  ST_tax_for_Tax4Fun > Tax4Fun_input.csv
      sed 's/,/\t/g' Tax4Fun_input.csv >  Tax4Fun_input.txt
      cp Tax4Fun_input.txt Tax4Fun_input.csv
      ```

      * Remember to remove the Taxonomy column name in the first row

      * Now ready to transfer to local computer to do Tax4Fun
      
3. The remaining just follow the R code.
  
