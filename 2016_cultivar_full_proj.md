#                                             2016 cultivar project

###                                               08_09_2017


-------

#### Raw data description: 
1. Six soybean cultivars were investigated in this study, including *Williams* (CV1), *Drought tolerant cultivar* (CV2), *soybean cyst nematode resistance cultivar* (CV3), *Williams non-nodulating mutant* (CV4), *Wild type - Glycine soja* (CV5) and *Williams 82* (CV6).
2. Two soil types were simutaneously investigated in this study, including agriculture and forest soil.

#### Upload raw_data into ACF server and get fastq files ready.

* Log into acf by ```ssh UserNetID@duo.acf.tennessee.edu``` with passcode and Duo mobile push code.
* Navigate to my project directory `` cd /nics/d/home/fliu21/16S_cultivar_project ``
* Creat a directory named ``raw_data``
* copy the folder ``Fang8-2-44439645``from local computer to acf project directory by:
```
scp -r -P 22 /g/UT/I_love_my_project/5_soybean_16S_data_and_analysis_fifth_run/cultivar_project/Fang8-2-44439645/ fliu21@duo.acf.tennessee.edu:/nics/d/home/fliu21/16S_cultivar_project/raw_data
```
Once secure copy is done, type in:

``ls``
```
A10-53762787  B6-53766768   D2-53761796   F10-53767768  G6-53753922
A11-53759854  B7-53753923   D3-53766765   F11-53761812  G7-53756890
A12-53759856  B8-53765776   D4-53768738   F12-53764779  G8-53767766
...           ...           ...           ...           ...
...           ...           ...           ...           ...
```
Each of the above folder consist of zipped fastq file of each sample.
* Move all of the fastq zip file from each folder to ``/nics/d/home/fliu21/16S_cultivar_project/raw_fastq/`` directory.

```
cd /nics/d/home/fliu21/16S_cultivar_project/raw_data
for file in *; do echo $file; cd $file; mv * ../../../raw_fastq/; cd ..;done
```
* Unzip all of the fastq file inside of raw_fastq direcotry by
```
gunzip /nics/d/home/fliu21/16S_cultivar_full_project/raw_fastq/*
```
* Change file names by:

```
for file in *; do echo $file; mv $file $(echo $file | sed 's/Ag/Ag_/g'); done
for file in *; do echo $file; mv $file $(echo $file | sed 's/For/For_/g'); done
for file in *; do echo $file; mv $file $(echo $file | sed 's/-/_/g'); done
```
* Copy all fastq file of bulk soils samples to raw_fastq directory
```
cp /nics/d/home/fliu21/16S_cultivar_proj/4th_run/raw_fastq/*.fastq  /nics/d/home/fliu21/16S_cultivar_full_project/raw_fastq/
```
* At the same time, cp all of the remaining bulk samples and several rhizosphere samples to this ``raw_fastq`` folder.

#### 
1. Make ``cultivar.file``, which is used for mothur to know what input fastq files are and which sample they belong to.

```
cd /nics/d/home/fliu21/16S_cultivar_proj/raw_fastq

```


.
.
.
.
.
.
.
.



----------------
# Strigolactone - tips for subset shared file based on sample list

----------------

* Make a list of the sampleID you want to keep in the subset shared file
```
awk '{print $1 }' strigolactone.file | paste -s -d-
```
The output looks likes this, which is used to run ``make.shared`` with ``groups= "below sample ID linked with -"``
```
CT_10-CT_1-CT_2-CT_3-CT_4-CT_5-CT_6-CT_7-CT_8-CT_9-Gene10_10-Gene10_1-Gene10_2-Gene10_3-Gene10_4-Gene10_5-Gene10_6-Gene10_7-Gene10_8-Gene10_9-Gene1_10-Gene1_1-Gene1_2-Gene1_3-Gene14_10-Gene14_1-Gene14_2-Gene14_3-Gene14_4-Gene14_5-Gene14_6-Gene14_7-Gene14_8-Gene14_9-Gene1_4-Gene1_5-Gene1_6-Gene1_7-Gene1_8-Gene1_9-HYPIII_B_1-HYPIII_B_2-HYPIII_B_3-HYPIII_B_4-HYPIII_B_5-HYPIII_B_6-HYPIII_B_7-HYPIII_B_8
```
* Run make.shared to subset OTU table to groups of sample (In my case, I want to filter PCR\_blanks)

```
make.shared(list=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list,count=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table,label=0.03,groups=CT_10-CT_1-CT_2-CT_3-CT_4-CT_5-CT_6-CT_7-CT_8-CT_9-Gene10_10-Gene10_1-Gene10_2-Gene10_3-Gene10_4-Gene10_5-Gene10_6-Gene10_7-Gene10_8-Gene10_9-Gene1_10-Gene1_1-Gene1_2-Gene1_3-Gene14_10-Gene14_1-Gene14_2-Gene14_3-Gene14_4-Gene14_5-Gene14_6-Gene14_7-Gene14_8-Gene14_9-Gene1_4-Gene1_5-Gene1_6-Gene1_7-Gene1_8-Gene1_9-HYPIII_B_1-HYPIII_B_2-HYPIII_B_3-HYPIII_B_4-HYPIII_B_5-HYPIII_B_6-HYPIII_B_7-HYPIII_B_8)
```
* After this process, all of the shared file including opti-mcc.shared and tx.shared are all filtered off PCR blank samples.

** make.lefse** & **make.biom** process could also filter off samples that we are not interested in for downward analysis. In my case, Gene14\_10 sample has tiny reads, So I want to filter off this outlier from the whole dataset.

1) make.lefse 
```
make.lefse(shared=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.shared,constaxonomy=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.cons.taxonomy,label=0.03,design=design.txt,groups=CT_10-CT_1-CT_2-CT_3-CT_4-CT_5-CT_6-CT_7-CT_8-CT_9-Gene10_10-Gene10_1-Gene10_2-Gene10_3-Gene10_4-Gene10_5-Gene10_6-Gene10_7-Gene10_8-Gene10_9-Gene1_10-Gene1_1-Gene1_2-Gene1_3-Gene14_1-Gene14_2-Gene14_3-Gene14_4-Gene14_5-Gene14_6-Gene14_7-Gene14_8-Gene14_9-Gene1_4-Gene1_5-Gene1_6-Gene1_7-Gene1_8-Gene1_9-HYPIII_B_1-HYPIII_B_2-HYPIII_B_3-HYPIII_B_4-HYPIII_B_5-HYPIII_B_6-HYPIII_B_7-HYPIII_B_8)
```
2) make.biom
```
make.biom(shared=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.shared,constaxonomy=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.cons.taxonomy,metadata=design.txt,groups=CT_10-CT_1-CT_2-CT_3-CT_4-CT_5-CT_6-CT_7-CT_8-CT_9-Gene10_10-Gene10_1-Gene10_2-Gene10_3-Gene10_4-Gene10_5-Gene10_6-Gene10_7-Gene10_8-Gene10_9-Gene1_10-Gene1_1-Gene1_2-Gene1_3-Gene14_1-Gene14_2-Gene14_3-Gene14_4-Gene14_5-Gene14_6-Gene14_7-Gene14_8-Gene14_9-Gene1_4-Gene1_5-Gene1_6-Gene1_7-Gene1_8-Gene1_9-HYPIII_B_1-HYPIII_B_2-HYPIII_B_3-HYPIII_B_4-HYPIII_B_5-HYPIII_B_6-HYPIII_B_7-HYPIII_B_8,label=0.03)
```
Both ``make.biom`` and ``make.lefse`` are used to combine ``shared`` and ``constaxonomy`` information together with meta data (which is used for differential abundance analysis, biomarker detection and circular taxonomy tree as well as network analysis- Cytoskape).  







