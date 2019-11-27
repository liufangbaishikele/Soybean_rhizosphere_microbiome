Author|Date
---| ----
Fang Liu| 10/08/19

### Prediction of potential function of bacteria community based on 16S rRNA amplicon sequencing

**Personal user experience about Tax4Fun, Tax4Fun2 and PICRUST2**

* I used Tax4Fun and Tax4Fun2 with my other dataset before. Tax4Fun is good but very picky with the input file format and the documentation is not prepared for the potential errors. Another problem with Tax4Fun is that the SILVA Tax4Fun reference is not always consistently keep updated with the most updated SILVA database. As a results, the precalculated KEGG annotation file is not update either. Tax4Fun2 (which is still a beta test version) is too slow, and have so many other dependense. For user do not have much coding experience or R using experience, it will be hard to solve the error assoicated with install associated packages especially BLAST.

* I never used PICRUST(version1) because it required greengene based taxonomy information (which I think is not as informative as SILVA as the reference for 16S classification). As I have been using mothur, I have to deal with converting shared and taxonomy as well as my meta data into biom form (which is pretty tricky based on my experience because mothur only support biom1 format which is not jason format).

* However, PICRUST2 did great job by integrating more bacteria reference genome for gene function annotation. In addition, tehy keep updated with KEGG database.

* Allows output of MetaCyc ontology predictions that will be comparable with common shotgun metagenomics outputs


* PICRUST2 is more flexible with input and compatible with mothur OTU represent sequence together with count table, I changed Tax4Fun2 to PICRUST2 for bacteria function prediction

All the below codes followed [PICRUST2 workflow](https://github.com/picrust/picrust2/wiki)

**This is the workflow**

1. Generate the representative OTU sequence using mothur command ``get.oturep``

```
get.oturep(column=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist,list=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list,count=Seed.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table,fasta=Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta,label=0.03)
```

2. After generate the represent sequence for each OTU, I did some modification of the output to make it accomodate with PICRUST2 example sequence data. Below are the codes needed for file formating.

    1. Use mothur command to generated the OTU_rep.fasta file
    2. Edit the rep.fasta file to only keep SeqID and OTUID for the header.
    ```
    awk -F "|" '{print $1}' Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.rep.fasta_copy | sed '2~2 s/-//g' > Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.rep.edit.fasta  
    ```

    3.  Create the OTUID list

    ```
    head -1 Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.subsample.0.03.pick.shared | sed 's/\t/\n/g' | grep '^Otu' > Subsampled_OTUID

    ```

    4. Substract rep.edit.fasta to only include those contained in the subsample.pick shared file

    ```{r}
    for otu in $(cat Subsampled_OTUID); do echo $otu; grep -A 1 $otu Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.rep.edit.fasta >> Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.rep.edit.subsample.pick.fasta; done  
    ```

    5. update the header to only include OTUID and remove the seqID

    ```
    sed 's/^>M.*\(.Otu.*\)/\1/g'  Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.rep.edit.subsample.pick.fasta | sed 's/\tOtu/>Otu/g' > Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.rep.edit.subsample.pick_up.fasta
    ```

    6. Transpose shared file and replace space with tab

    ```
    awk '                                                                                                {                                                                                                                                 for (i=1; i<=NF; i++)  {
            a[NR,i] = $i
        }
    }
    NF>p { p = NF }
    END {    
        for(j=1; j<=p; j++) {
            str=a[1,j]
            for(i=2; i<=NR; i++){
                str=str" "a[i,j];
            }
            print str
        }
    }' Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.subsample.0.03.pick.shared > Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.subsample.0.03.pick_up.shared
    ```

    ```
    sed 's/ /\t/g' Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.subsample.0.03.pick_up.shared > Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.subsample.0.03.pick_up_tsv.shared
    ```

    7. Edit the above tsv file to the right format, only include OTUID and count in each sample using vim

    ```
    sed 's/ /\t/g' Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc_up.shared > Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc_up_tsv.shared
    ```
    
**NOTES** Please found the input [fasta sequence file and shared count files within this Archive](https://github.com/liufangbaishikele/Soybean_rhizosphere_microbiome/blob/master/2018_fungicide/Seed_treatment/16S/PICRUST/Archive.zip) attached as reference for guiding you to format the input in a correct way.


3. Once the above input rep.fasta and shared count files were ready, run PICRUST2 commands

```
picrust2_pipeline.py -s Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.rep.edit.subsample.pick_up.fasta  -i Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.subsample.0.03.pick.shared  -o picrust2_out_pipeline -p 30
```
4. For the following differential pathway abundance analysis, the output file from the pathways_out folder were copied to local computer for further analysis using R. Please find [detailed R command](https://github.com/liufangbaishikele/Soybean_rhizosphere_microbiome/blob/master/2018_fungicide/Seed_treatment/16S/PICRUST/PICRUST_KO_enrichment.Rmd) 


