#                                     Community analysis summary


----
Fang Liu
08/02/2019


**Overal description**

This is the results of Combined run of ITS sequencing targeted ITS2 region. For the first run, we got very limited reads left for endosphere and some rhizosphere samples. So, we design ITS blocker to block off soybean ITS, but it did not work as good as we expected. However, for most of the samples we got enough reads for analysis.

## PERMANOVA results

### All samples without rarefaction, but singletons were removed

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8360 taxa and 164 samples ]
sample_data() Sample Data:       [ 164 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 8360 taxa by 7 taxonomic ranks ]
```

* **Read_depth**

```
adonis2(formula = t(otu_table(r_rms_seed_phyloseq)) ~ Read_depth, data = data.frame(sample_data(r_rms_seed_phyloseq)), permutations = 999, by = "margin")
            Df SumOfSqs      R2      F Pr(>F)    
Read_depth   1    3.792 0.08548 15.142  0.001 ***
Residual   162   40.566 0.91452                  
Total      163   44.358 1.00000                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

* **Compartment** --- significant p=0.001 and R2=38.00%

```
aadonis2(formula = t(otu_table(r_rms_seed_phyloseq)) ~ Compartment, data = data.frame(sample_data(r_rms_seed_phyloseq)), permutations = 999, by = "margin")
             Df SumOfSqs      R2     F Pr(>F)    
Compartment   4   17.472 0.37988 24.35  0.001 ***
Residual    159   28.521 0.62012                 
Total       163   45.993 1.00000  
```
* **Read_depth vs Compartment** -- Read_depth is still significant, but Compartment is the driving factor after parsing out the impact of sequencing depth 

```
adonis2(formula = t(otu_table(r_rms_seed_phyloseq)) ~ Read_depth + Compartment, data = data.frame(sample_data(r_rms_seed_phyloseq)), permutations = 999, by = "margin")
             Df SumOfSqs      R2       F Pr(>F)    
Read_depth    1    0.407 0.00884  2.2856  0.009 ** 
Compartment   4   11.037 0.23998 15.5069  0.001 ***
Residual    158   28.115 0.61128                   
Total       163   45.993 1.00000 
```


### After remove all the seed samples: non_seed samples

* **Read_Depth**

```
adonis2(formula = t(otu_table(r_rms_non_Seed_up)) ~ Read_depth, data = data.frame(sample_data(r_rms_non_Seed_up)), permutations = 999, by = "margin")
            Df SumOfSqs     R2      F Pr(>F)    
Read_depth   1    2.644 0.0749 11.902  0.001 ***
Residual   147   32.661 0.9251                  
Total      148   35.306 1.0000
```

### Subset to rms_Soil_up 
In order to check the clustering of plots before sowing soybean seeds

**Treatment impact** -- insignificant

```
adonis2(formula = t(otu_table(r_rms_Soil_up)) ~ Treatment, data = data.frame(sample_data(r_rms_Soil_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)
Treatment  2  0.19993 0.14733 1.0367  0.359
Residual  12  1.15709 0.85267              
Total     14  1.35702 1.00000
```
### Subset to BRE - Bulk, rhizosphere and endosphere 

* **Read_depth**

```
adonis2(formula = t(otu_table(r_rms_BRE_up)) ~ Read_depth, data = data.frame(sample_data(r_rms_BRE_up)), permutations = 999, by = "margin")
            Df SumOfSqs     R2      F Pr(>F)    
Read_depth   1    2.554 0.0789 11.307  0.001 ***
Residual   132   29.818 0.9211                  
Total      133   32.372 1.0000
```

* **Plot impacts**

```
adonis2(formula = t(otu_table(r_rms_BRE_up)) ~ Plot, data = data.frame(sample_data(r_rms_BRE_up)), permutations = perm, by = "margin")
Model: adonis0(formula = lhs ~ Plot, data = data, method = method)
          Df SumOfSqs      R2      F Pr(>F)    
Plot      14    5.192 0.16037 1.6236  0.001 ***
Residual 119   27.180 0.83963                  
Total    133   32.372 1.00000 
```
* **Compartment**

```
adonis2(formula = t(otu_table(r_rms_BRE_up)) ~ Compartment, data = data.frame(sample_data(r_rms_BRE_up)), permutations = perm, by = "margin")
Model: adonis0(formula = lhs ~ Compartment, data = data, method = method)
             Df SumOfSqs      R2      F Pr(>F)    
Compartment   2   10.182 0.31452 30.054  0.001 ***
Residual    131   22.190 0.68548                  
Total       133   32.372 1.00000
```

* **Time impact**

```
Marginal effects of terms
Plots: paste(data.frame(sample_data(r_rms_BRE_up))$Compartment, data.frame(sample_data(r_rms_BRE_up))$Plot, , plot permutation: none
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(r_rms_BRE_up)) ~ Time, data = data.frame(sample_data(r_rms_BRE_up)), permutations = perm, by = "margin")
Model: adonis0(formula = lhs ~ Time, data = data, method = method)
          Df SumOfSqs      R2      F Pr(>F)    
Time       2    1.885 0.05824 4.0504  0.001 ***
Residual 131   30.487 0.94176                  
Total    133   32.372 1.00000 
```

* **Read_depth vs Compartment vs Time vs Plot**

```
adonis2(formula = t(otu_table(r_rms_BRE_up)) ~ Read_depth + Compartment + Time + Plot, data = data.frame(sample_data(r_rms_BRE_up)), permutations = 999, by = "margin")
             Df SumOfSqs      R2       F Pr(>F)    
Read_depth    1    0.274 0.00846  2.1010  0.033 *  
Compartment   2    7.789 0.24060 29.8663  0.001 ***
Time          2    1.887 0.05830  7.2367  0.001 ***
Plot         14    5.179 0.15997  2.8368  0.001 ***
Residual    114   14.865 0.45919                   
Total       133   32.372 1.00000 
```
* **Read_depth vs Compartment vs Time vs Treatment**

```
adonis2(formula = t(otu_table(r_rms_BRE_up)) ~ Read_depth + Compartment + Time + Treatment, data = data.frame(sample_data(r_rms_BRE_up)), permutations = 999, by = "margin")
             Df SumOfSqs      R2       F Pr(>F)    
Read_depth    1    0.273 0.00842  1.8023  0.059 .  
Compartment   2    7.853 0.24259 25.9574  0.001 ***
Time          2    1.893 0.05849  6.2580  0.001 ***
Treatment     2    0.983 0.03038  3.2506  0.001 ***
Residual    126   19.060 0.58878                   
Total       133   32.372 1.00000
```

### subset to Bulk -- r_rms_Bulk vs r_rf_rms_Bulk vs r_norm_Bulk. Here, rms means singletons were removed; rf means samples were rarefied to minimum read depth after remove singletons; 

**r_rms_Bulk**

* --- Read_depth ----
```
adonis2(formula = t(otu_table(r_rms_Bulk_up)) ~ Read_depth, data = data.frame(sample_data(r_rms_Bulk_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)  
Read_depth  1   0.2703 0.04328 1.9454  0.014 *
Residual   43   5.9736 0.95672                
Total      44   6.2439 1.00000
```

* ---- Time ---

```
Permutation test for adonis under NA model
Marginal effects of terms
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(r_rms_Bulk_up)) ~ Time, data = data.frame(sample_data(r_rms_Bulk_up)), permutations = 999, by = "margin")
         Df SumOfSqs      R2      F Pr(>F)    
Time      2   1.3833 0.22155 5.9768  0.001 ***
Residual 42   4.8605 0.77845                  
Total    44   6.2439 1.00000
```
* ---- Plot -----

```
Permutation test for adonis under NA model
Marginal effects of terms
Plots: data.frame(sample_data(r_rms_Bulk_up))$Time, plot permutation: none
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(r_rms_Bulk_up)) ~ Plot, data = data.frame(sample_data(r_rms_Bulk_up)), permutations = perm, by = "margin")
         Df SumOfSqs      R2      F Pr(>F)    
Plot     14   2.4409 0.39092 1.3753  0.001 ***
Residual 30   3.8030 0.60908                  
Total    44   6.2439 1.00000
```
* ---- Treatment ----

```
Permutation test for adonis under NA model
Marginal effects of terms
Plots: data.frame(sample_data(r_rms_Bulk_up))$Time, plot permutation: none
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(r_rms_Bulk_up)) ~ Treatment, data = data.frame(sample_data(r_rms_Bulk_up)), permutations = perm, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)   
Treatment  2   0.3617 0.05793 1.2913  0.003 **
Residual  42   5.8822 0.94207                 
Total     44   6.2439 1.00000                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
* ----Read_depth vs Time vs Plot

```
adonis2(formula = t(otu_table(r_rms_Bulk_up)) ~ Read_depth + Plot + Time, data = data.frame(sample_data(r_rms_Bulk_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)    
Read_depth  1   0.0995 0.01593 1.1573  0.222    
Plot       14   2.4008 0.38450 1.9955  0.001 ***
Time        2   1.2131 0.19428 7.0581  0.001 ***
Residual   27   2.3202 0.37160                  
Total      44   6.2439 1.00000
```

* ---- Read_depth vs Time vs Treatment

```
adonis2(formula = t(otu_table(r_rms_Bulk_up)) ~ Read_depth + Treatment + Time, data = data.frame(sample_data(r_rms_Bulk_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)    
Read_depth  1   0.1505 0.02410 1.3496  0.112    
Treatment   2   0.3726 0.05968 1.6711  0.012 *  
Time        2   1.2374 0.19817 5.5489  0.001 ***
Residual   39   4.3484 0.69642                  
Total      44   6.2439 1.00000
```
                        
### Subset to Rhizosphere -- before rarefaction

* --- Plot----

```
adonis2(formula = t(otu_table(Rhi_up)) ~ Plot, data = data.frame(sample_data(Rhi_up)), permutations = 999, by = "margin")
         Df SumOfSqs      R2      F Pr(>F)   
Plot     14   3.2632 0.38842 1.3156  0.002 **
Residual 29   5.1380 0.61158                 
Total    43   8.4012 1.00000
```
* -----Read_depth---

```
adonis2(formula = t(otu_table(Rhi_up)) ~ Read_depth, data = data.frame(sample_data(Rhi_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)    
Read_depth  1   1.1971 0.14249 6.9791  0.001 ***
Residual   42   7.2041 0.85751                  
Total      43   8.4012 1.00000
```

* ---- Time ----

```
adonis2(formula = t(otu_table(Rhi_up)) ~ Time, data = data.frame(sample_data(Rhi_up)), permutations = 999, by = "margin")
         Df SumOfSqs      R2      F Pr(>F)    
Time      2   1.6706 0.19885 5.0883  0.001 ***
Residual 41   6.7306 0.80115                  
Total    43   8.4012 1.00000
```

* --- Read_depth vs Time vs Plot 

```
adonis2(formula = t(otu_table(Rhi_up)) ~ Read_depth + Time + Plot, data = data.frame(sample_data(Rhi_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)    
Read_depth  1   0.3648 0.04343 3.0389  0.001 ***
Time        2   1.0049 0.11961 4.1852  0.001 ***
Plot       14   2.9480 0.35090 1.7540  0.001 ***
Residual   26   3.1213 0.37153                  
Total      43   8.4012 1.00000
```
* ---- Read_depth vs Time vs Treatment

```
adonis2(formula = t(otu_table(Rhi_up)) ~ Read_depth + Time + Treatment, data = data.frame(sample_data(Rhi_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)    
Read_depth  1   0.6631 0.07893 4.6724  0.001 ***
Time        2   1.1289 0.13438 3.9774  0.001 ***
Treatment   2   0.6764 0.08051 2.3831  0.001 ***
Residual   38   5.3929 0.64192                  
Total      43   8.4012 1.00000
```

* ---- Treatment impact along Rhi_wk1, Rhi_wk3 and Rhi_wk4 ---

```
>>>> Rhi_1wk --- p=0.09 and R2=18.58%

adonis2(formula = t(otu_table(r_Rhi_1wk_up)) ~ Treatment, data = data.frame(sample_data(r_Rhi_1wk_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)  
Treatment  2  0.34943 0.18584 1.2554   0.09 .
Residual  11  1.53089 0.81416                
Total     13  1.88032 1.00000

>>>> Rhi_3wk --- significant and R=27.99%

adonis2(formula = t(otu_table(Rhi_3wk_up)) ~ Treatment, data = data.frame(sample_data(Rhi_3wk_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)   
Treatment  2  0.61565 0.27987 2.3319  0.002 **
Residual  12  1.58409 0.72013                 
Total     14  2.19974 1.00000

>>>> Rhi_4wk ---

adonis2(formula = t(otu_table(Rhi_4wk_up)) ~ Treatment, data = data.frame(sample_data(Rhi_4wk_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)
Treatment  2  0.38076 0.16041 1.1463   0.23
Residual  12  1.99293 0.83959              
Total     14  2.37369 1.00000 
```

### Subset to Rhizosphere -- After rarefaction

* ---- Read_depth ---

```
adonis2(formula = t(otu_table(rf_Rhi_up)) ~ Read_depth, data = data.frame(sample_data(rf_Rhi_up)), permutations = 999, by = "margin")
           Df SumOfSqs     R2      F Pr(>F)   
Read_depth  1   0.4240 0.0563 2.5055  0.003 **
Residual   42   7.1084 0.9437                 
Total      43   7.5325 1.0000
```
* ----- Plot ----

```
adonis2(formula = t(otu_table(rf_Rhi_up)) ~ Plot, data = data.frame(sample_data(rf_Rhi_up)), permutations = 999, by = "margin")
         Df SumOfSqs      R2      F Pr(>F)    
Plot     14   3.0188 0.40077 1.3854  0.001 ***
Residual 29   4.5137 0.59923                  
Total    43   7.5325 1.00000
```

* ---- Time ---

```
adonis2(formula = t(otu_table(rf_Rhi_up)) ~ Time, data = data.frame(sample_data(rf_Rhi_up)), permutations = 999, by = "margin")
         Df SumOfSqs      R2      F Pr(>F)    
Time      2   1.4676 0.19484 4.9607  0.001 ***
Residual 41   6.0649 0.80516                  
Total    43   7.5325 1.00000
```

* ---- Read_depth vs Time vs Plot --

```
adonis2(formula = t(otu_table(rf_Rhi_up)) ~ Read_depth + Time + Plot, data = data.frame(sample_data(rf_Rhi_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)    
Read_depth  1   0.1652 0.02193 1.4850  0.086 .  
Time        2   1.0548 0.14004 4.7407  0.001 ***
Plot       14   3.0343 0.40283 1.9481  0.001 ***
Residual   26   2.8926 0.38402                  
Total      43   7.5325 1.00000
```

* --- Read_depth vs Time vs Treatment ---

```
adonis2(formula = t(otu_table(rf_Rhi_up)) ~ Read_depth + Time + Treatment, data = data.frame(sample_data(rf_Rhi_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)    
Read_depth  1   0.1527 0.02027 1.1083  0.286    
Time        2   1.1945 0.15858 4.3350  0.001 ***
Treatment   2   0.6914 0.09179 2.5091  0.001 ***
Residual   38   5.2355 0.69506                  
Total      43   7.5325 1.00000
```

* --- Treatment impact along Rhi_1wk, Rhi_3wk and Rhi_4wk ----

```
>>> rf_Rhi_1wk

adonis2(formula = t(otu_table(r_rf_Rhi_1wk_up)) ~ Treatment, data = data.frame(sample_data(r_rf_Rhi_1wk_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)
Treatment  2  0.35167 0.18056 1.2119  0.127
Residual  11  1.59606 0.81944              
Total     13  1.94773 1.00000

>>> rf_Rhi_3wk

adonis2(formula = t(otu_table(rf_Rhi_3wk_up)) ~ Treatment, data = data.frame(sample_data(rf_Rhi_3wk_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)   
Treatment  2  0.56049 0.27167 2.2381  0.003 **
Residual  12  1.50261 0.72833                 
Total     14  2.06310 1.00000

>>> rf_Rhi_4wk

adonis2(formula = t(otu_table(rf_Rhi_4wk_up)) ~ Treatment, data = data.frame(sample_data(rf_Rhi_4wk_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)
Treatment  2   0.3185 0.15982 1.1413  0.225
Residual  12   1.6744 0.84018              
Total     14   1.9929 1.00000

```

### Subset to Endosphere - without rarefaction due to the very limited sequencing depth

* ---- Plot ----

```
adonis2(formula = t(otu_table(Endo)) ~ Plot, data = data.frame(sample_data(Endo)), permutations = 999, by = "margin")
         Df SumOfSqs    R2      F Pr(>F)    
Plot     14   4.7811 0.382 1.3246  0.001 ***
Residual 30   7.7348 0.618                  
Total    44  12.5159 1.000
```
* --- Read_depth ---

```
adonis2(formula = t(otu_table(Endo)) ~ Read_depth, data = data.frame(sample_data(Endo)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)    
Read_depth  1   1.4133 0.11292 5.4738  0.001 ***
Residual   43  11.1026 0.88708                  
Total      44  12.5159 1.00000
```
* --- Time ----

```
adonis2(formula = t(otu_table(Endo)) ~ Time, data = data.frame(sample_data(Endo)), permutations = 999, by = "margin")
         Df SumOfSqs      R2      F Pr(>F)    
Time      2   1.6394 0.13098 3.1652  0.001 ***
Residual 42  10.8766 0.86902                  
Total    44  12.5159 1.00000
```


* ---- Read_depth vs Time vs Plot

```
adonis2(formula = t(otu_table(Endo)) ~ Plot + Time + Read_depth, data = data.frame(sample_data(Endo)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)    
Plot       14   4.2099 0.33637 1.4506  0.001 ***
Time        2   1.1834 0.09455 2.8542  0.001 ***
Read_depth  1   0.4983 0.03981 2.4037  0.001 ***
Residual   27   5.5972 0.44720                  
Total      44  12.5159 1.00000
```

* --- Read_depth vs Time vs Treatment

```
Treatment   2   0.6602 0.05275 1.4074  0.021 *  
Time        2   1.2930 0.10331 2.7565  0.001 ***
Read_depth  1   0.9004 0.07194 3.8392  0.001 ***
Residual   39   9.1469 0.73082                  
Total      44  12.5159 1.00000
```

* ---- Treatment impact along soybean development stages

```

>>> Endo_1wk

adonis2(formula = t(otu_table(r_Endo_1wk_up)) ~ Treatment, data = data.frame(sample_data(r_Endo_1wk_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)
Treatment  2   0.4300 0.13369 0.9259  0.681
Residual  12   2.7865 0.86631              
Total     14   3.2165 1.00000

>>> Endo_3wk

adonis2(formula = t(otu_table(Endo_3wk_up)) ~ Treatment, data = data.frame(sample_data(Endo_3wk_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)
Treatment  2   0.5600 0.17161 1.2429  0.129
Residual  12   2.7035 0.82839              
Total     14   3.2635 1.00000

>>> Endo_4wk

adonis2(formula = t(otu_table(Endo_4wk_up)) ~ Treatment, data = data.frame(sample_data(Endo_4wk_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)
Treatment  2   0.4844 0.13864 0.9658  0.546
Residual  12   3.0092 0.86136              
Total     14   3.4936 1.00000
```


### Subset to SEED samples

* ---- Treatment impact ----

```
adonis2(formula = t(otu_table(SEED)) ~ Treatment, data = data.frame(sample_data(SEED)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)
Treatment  2   0.5141 0.16021 1.1446   0.25
Residual  12   2.6949 0.83979              
Total     14   3.2090 1.00000
```


### Marjor conclusions

* For fungi community soybean development impacts are comparatively low compared with bacteria community, with 16%, 12% and 9% variations could be explained by collection time for Bulk, Rhizosphere and Endosphere samples.

* Fungi community of presowing soil samples are similar to 1wk/4wk bulk soil samples. 

* Collection time impact could be indicated by fungi community difference between different collection time, including presowing, Bulk_1wk and Bulk_3wk as well as Bulk_4wk. We found distinct microbial community of Bulk_3wk compared with Bulk_1wk, Bulk_4wk and Pre-sowing samples.

* The community difference of bulk soil was strongly impacted by read_depth. Fungi community difference along soybean development stages changes and this change pattern differs between pre-rarefied data and rarefied data.

* Soil, Bulk and Rhizosphere were subseted to quantify the Read_depth vs Time vs Plot/Treatment impacts. Below are the results: Time impacts are decreasing from bulk soil to rhizosphere to endosphere samples. 

**Fungi** --- Decreasing time impacts from bulk to rhizosphere and to endosphere compartment but similar between plots difference.
```
Variations explained by factors			
	Bulk	Rhizosphere	Endosphere
Read_depth	3.50%	4.34%	3.98%
Time	16.50%	12%	9.45%
Plot	34.80%	35.10%	33.60%
```
**Bacteria** --- Increasing time impacts but decreading between plots difference.

```
Variations could be explained by each factors			
	Bulk	Rhizosphere	Endosphere
Read_depth	1.26%	0.75%	2.01%
Plots	48.80%	27.70%	17.80%
Time	11.30%	38.40%	58.40%
```

* The above comparision indicated that soybean development impacts fungi community to a similar degree across bulk, rhizosphere and endosphere compartments. 

* We saw significantly difference fungi community between bulk, rhizosphere and endosphere compartments. But soybean development did not change rhizosphere and endosphere microbiome composition that much.

* In contrast, soybean development impacts on rhizosphere and endosphere bacteria community are way stronger compared with bulk samples. This indicate that soybean development strongly modulate soybean rhizosphere and endosphere bacteria community.

* Plots impacts on bulk fungi community is larger than bacteria community. Which indicate more divergent bacteria community composition between plots. However, this plots difference decreasingly very significantly for bacteria community for rhizosphere and endosphere samples. This indicate a very consistent and robust selection of bacterial community occurs in rhizosphere and endosphere despite the very divergent indigenous bacteria pools. 

* In comparison, Plots impacts on soybean fungi community are similar between different compartments. This indirectly indicate that the modulation force of soybean rhizosphere and endosphere on fungi community are very limited and probably more depend on the indigenous fungi pool instead of a strong and consistent selection on specific taxa. More random selection of fungi taxa occurs in soybean rhizosphere and endospehre compartment maybe.  


* Plot difference of bulk bacteria community are highere than that on fungi community. However, this Plot impacts decreased for bacteria rhizospehre and endosphere impacts (which indicate more robust rhizosphere and endosphere bacteria recruitment despites the indigenous bacteria community difference).For fungi community, Plots difference on bulk soil are less compared with that on bacteria community. But this plots differences did not change for rhizosphere and endospehre compartments. Which indicate that their are less robust recruitment of fungi taxa occured in rhizosphere and endosphere. 


