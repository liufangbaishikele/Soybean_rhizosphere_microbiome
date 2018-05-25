#                               2016 strigolactone LefSe analysis documentation
-----------------

## Data manipulation within mothur -
The reason that I manipulate shared file inside of mothur instead of within R is that make.lefse command need .shared file as input.

-------------------------

### Genus level shared file -- taxlevel=1

Forder path

```
/nics/d/home/fliu21/16S_strigolactone_proj/LefSe/LefSe_st_genus_update
```

1. Rarefy\_based normalization - Subset ``.shared`` using rarefaction method to minimum read depth (26050) along all samples using ``sub_sample``

```
sub.sample(shared=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.shared_original,size=26050,groups=CT_1-CT_2-CT_3-CT_4-CT_5-CT_6-CT_7-Gene10_1-Gene10_2-Gene10_3-Gene10_4-Gene10_5-Gene10_6-Gene10_7-Gene1_1-Gene1_2-Gene1_3-Gene14_1-Gene14_2-Gene14_3-Gene14_4-Gene14_5-Gene14_6-Gene14_7-Gene1_4-Gene1_5-Gene1_6-Gene1_7-HYPIII_B_1-HYPIII_B_2-HYPIII_B_3-HYPIII_B_4-HYPIII_B_5-HYPIII_B_6-HYPIII_B_7-HYPIII_B_8)
```
2. Remove genus with total abundance less than 50. Using ``remove.rare`` function in mothur


```
remove.rare(shared=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.1.subsample.shared_original,nseqs=50)
```
3. Creat lefse object using ``make.lefse``

```
make.lefse(shared=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.1.subsample.1.pick.shared_original,constaxonomy=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.cons.taxonomy_original,design=design.txt,scale=totalgroup)
```
4. Now, change to beacon from mothur. Vim generated *.lefse* file, there are lots of double quotes inside. Using vim replacement to replace all doucle quote with empty
```
$%s/"//g
```
5. Add path to the environment ``export PATH=/lustre/medusa/fliu21/anaconda2/bin:$PATH``

6. Run ``lefse-format_lefse.py`` to change format for running real lefse analysis

**At this step, -o pamameter is very import** 

``-o : set the normalization value. The default set is -1 means no normalization``

  * 1st (-o 1000000)

   ```
   lefse-format_input.py strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.1.subsample.1.pick.1.lefse strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.1.subsample.1.pick.1.lefse.in_1 -c 1 -s -1 -u 2 -o 1000000
   ```
  * 2nd (-o 100)

   ```
   lefse-format_input.py strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.1.subsample.1.pick.1.lefse strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.1.subsample.1.pick.1.in_2 -c 1 -s -1 -u 2 -o 100
   ```
   * 3rd (-o=10000)
   
   ```
   lefse-format_input.py strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.1.subsample.1.pick.1.lefse strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.1.subsample.1.pick.1.in_3 -c 1 -s -1 -u 2 -o 10000
   ```
   * 4th (-o=100000)
   ```
   lefse-format_input.py strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.1.subsample.1.pick.1.lefse strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.1.subsample.1.pick.1.in_4 -c 1 -s -1 -u 2 -o 100000
   ```
   
   

7. Doing lefse analysis using ``run_lefse.py``

  * y=1 (parameter is used to set whether the test is performed in one against one or in one again all setting. default is one against all y = 0)
  
   ```
   (for multiclass tasks) set whether the test is performed in a one-against-one ( 1 - more strict!) or in a one-against-all setting ( 0 - less strict) (default 0)
   ```
   
   Inside of the paper: here are some description of this parameter:
   
   ```
   LEfSe allows ordinal classes with more than two levels to be analyzed in two different stringencies. The first requires significant taxa to differ between every pair of class values ; the discovered biomarkers must accurately distinguish all individual classes - strict version. Alternatively, LEfSe can determine significant taxa differing in at least one class value(s) (non-strict version).
   ```
  
    1. 
    ```
    run_lefse.py   strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.1.subsample.1.pick.1.lefse.in_1 strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.1.subsample.1.pick.1.lefse.res_1_l2 -l 2 -b 100 -t strigolactone_genus_lefse -y 1

    Number of significantly discriminative features: 719 ( 756 ) before internal wilcoxon
    Number of discriminative features with abs LDA score > 2.0 : 719
    
    run_lefse.py   strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.1.subsample.1.pick.1.lefse.in_1 strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.1.subsample.1.pick.1.lefse.res_1_l4 -l 4 -b 100 -t strigolactone_genus_lefse -y 1
    Number of significantly discriminative features: 719 ( 756 ) before internal wilcoxon
    Number of discriminative features with abs LDA score > 4.0 : 50
    ```
    2. 
    ```
    run_lefse.py strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.1.subsample.1.pick.1.in_2 strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.1.subsample.1.pick.1.res_2 -l 2 -b 100 -t strigolactone_genus_LefSe -y 1
    Number of significantly discriminative features: 719 ( 756 ) before internal wilcoxon
    Number of discriminative features with abs LDA score > 2.0 : 0
    ```
    3.
    ```
    run_lefse.py strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.1.subsample.1.pick.1.in_3 strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.1.subsample.1.pick.1.res_3 -l 2 -b 100 -t strigolactone_genus_LefSe -y 1
    
    Number of significantly discriminative features: 719 ( 756 ) before internal wilcoxon
    Number of discriminative features with abs LDA score > 2.0 : 50
    ```
    4. plot LDA stem\_leaf plot using ``lefse-plot_res.py``
    
    ```
    lefse-plot_res.py  strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.1.subsample.1.pick.1.res_3 strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.1.subsample.1.pick.1.plot_3.png --dpi 300 --title strigolactone_genus_lefse
    Number of significantly discriminative features: 719 ( 756 ) before internal wilcoxon
    Number of discriminative features with abs LDA score > 2.0 : 394
    ```












