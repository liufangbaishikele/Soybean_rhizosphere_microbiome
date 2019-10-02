Author|Date
---| ----
Fang | Liu

### Prediction of potential function of bacteria community based on 16S rRNA amplicon sequencing

**Personal user experience about Tax4Fun, Tax4Fun2 and PICRUST2**

* I used Tax4Fun and Tax4Fun2 with my other dataset before. Tax4Fun is good but very picky with the input file format and the documentation is not prepared for the potential errors. Another problem with Tax4Fun is that the SILVA Tax4Fun reference is not always consistently keep updated with the most updated SILVA database. As a results, the precalculated KEGG annotation file is not update either. Tax4Fun2 (which is still a beta test version) is too slow, and have so many other dependense. For user do not have much coding experience or R using experience, it will be hard to solve the error assoicated with install associated packages especially BLAST.

* I never used PICRUST(version1) because it required greengene based taxonomy information (which I think is not as informative as SILVA as the reference for 16S classification). As I have been using mothur, I have to deal with converting shared and taxonomy as well as my meta data into biom form (which is pretty tricky based on my experience because mothur only support biom1 format which is not jason format).

* However, PICRUST2 did great job by integrating more bacteria reference genome for gene function annotation. In addition, tehy keep updated with KEGG database.

* Allows output of MetaCyc ontology predictions that will be comparable with common shotgun metagenomics outputs


* PICRUST2 is more flexible with input and compatible with mothur OTU represent sequence together with count table, I changed Tax4Fun2 to PICRUST2 for bacteria function prediction

All the below codes followed [PICRUST2 workflow](https://github.com/picrust/picrust2/wiki)

**This is the workflow**

1. Generate 
