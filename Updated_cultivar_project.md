#                                                Cultivar Project from 2016 Data
###                                               08_09_2017


-------

#### Raw data description: 
1. Six soybean cultivars were investigated in this study, including *Williams* (CV1), *Drought tolerant cultivar* (CV2), *soybean cyst nematode resistance cultivar* (CV3), *Williams non-nodulating mutant* (CV4), *Wild type - Glycine soja* (CV5) and *Williams 82* (CV6).
2. 


#### Upload raw_data into ACF server

* Log into acf by ```ssh UserNetID@duo.acf.tennessee.edu``` with passcode and Duo mobile push code.
* Navigate to my project directory `` cd /nics/d/home/fliu21/16S_cultivar_full_project ``
* Creat a directory named ``raw_data``
* copy the folder ``Fang8-2-44439645``from local computer to acf project directory by:
```
scp -r -P 22 /g/UT/I_love_my_project/5_soybean_16S_data_and_analysis_fifth_run/cultivar_project/Fang8-2-44439645/ fliu21@duo.acf.tennessee.edu:/nics/d/home/fliu21/16S_cultivar_full_project/raw_data
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
* Move all fastq file of bulk soils samples to raw_fastq directory

```
cp /nics/d/home/fliu21/16S_cultivar_full_project/raw_fastq/*.fastq  /nics/d/home/fliu21/16S_cultivar_full_project/raw_fastq/
```

