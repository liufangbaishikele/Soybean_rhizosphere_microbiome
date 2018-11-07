##                            2016 cultivar project

> Date: 08_09_2017  
> Analyzer: Fang Liu  

-------

### Raw data description:
 
1. Six soybean cultivars were investigated in this study, including *Williams* (WIL), *Drought tolerant cultivar* (DRT), *soybean cyst nematode resistance cultivar* (CNR), *Williams non-nodulating mutant* (NNW), *Wild type - Glycine soja* (SOJ) and *Williams 82* (W82).  

2. Two soil types were simutaneously investigated in this study, including agriculture and forest soil.  

### Upload raw_data into ACF server and get fastq files ready.

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

###  Mothur process
 
1. Make ``cultivar.file``, which is used for mothur to know what input fastq files are and which sample they belong to.

```
cd /nics/d/home/fliu21/16S_cultivar_proj/raw_fastq
ls *R1_001.fastq > forward_ID
ls *R1_001.fastq > reverse_ID
sed 's/\([^.*]\)_S.*/\1/' forward_ID > treatment_ID
paste --delimiters="\t" treatment_ID forward_ID reverse)_ID > cultivar.files
```
2. make.contigs together with oligo file will only make contigs using extactly match of both forward and reverse primer

* Here is a quick example illustrating how oligo works
  * There are 167705 raw reads in this sample
  * make.contigs(ffastq=Ag_B_10_S13_L001_R1_001.fastq, rfastq=Ag_B_10_S13_L001_R2_001.fastq,processors=16)
  * the output trim.contigs file has 167,705 contigs.
 
  *BUT when I make.contigs together using oligo file*
  * make.contigs(ffastq=Ag_B_10_S13_L001_R1_001.fastq, rfastq=Ag_B_10_S13_L001_R2_001.fastq,oligos=cultivar.oligo,processors=16)
  * It turned out that only 149,099 contigs are left after make.contigs as those with mismatch with primer sequences were automatically discarded
  
3. The following are a quick summary of lost of reads along each step of mothur pipeline 

-----------------------------
**Sequence number summary**
-----------------------------
  Total read from 136 samples - 19,358,039
  
  After make.contigs - 15,946,467
  
  After trimming based on ambiguous base and read length - 12,247,497 
  
  After further screen using aligned reads – 11,773,645
  
  Romove Chimera – 9,967,343 
  
  Remove non-bacteria – 9,946,720
  
  OTUs- 175,957 (9,945,986 reads = 9,945,986/ 15,946,467 = 62.37% ; 9,945,986/19,358,039=51.37%)
  
.
.
.
.
.
.
.
.

### In terms of the detailed documentation, please find in the mothur.logfile**



