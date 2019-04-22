#                                     Community analysis summary


----
Fang Liu
04/13/2019


**Overal description**

This is the first run of ITS sequencing targeted ITS2 region. Due to the large propotion of soybean derived ITS, endosphere samples and seed samples were left with less informative sequences. Anyway, here are the summary of this run results.

## PERMANOVA results

### All samples without rarefaction

* **Read_depth**

```
adonis2(formula = t(otu_table(r_rms_seed_phyloseq)) ~ Read_depth, data = data.frame(sample_data(r_rms_seed_phyloseq)), permutations = 999, by = "margin")
            Df SumOfSqs      R2      F Pr(>F)    
Read_depth   1    6.841 0.14874 28.307  0.001 ***
Residual   162   39.152 0.85126                  
Total      163   45.993 1.00000                  
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

```
capscale(formula = t(otu_table(r_rms_seed_phyloseq)) ~ Read_depth + Condition(Compartment), data = data.frame(sample_data(r_rms_seed_phyloseq)), distance = "bray", dfun = vegdist) 

----- summary ------
Partitioning of squared Bray distance:
              Inertia Proportion
Total         46.8384   1.000000
Conditioned   17.4739   0.373067
Constrained    0.4088   0.008727
Unconstrained 28.9557   0.618205

----- anova results ----

Df SumOfSqs      F Pr(>F)    
Model      1   0.4088 2.2305  0.001 ***
Residual 158  28.9557 



capscale(formula = t(otu_table(r_rms_seed_phyloseq)) ~ Compartment + Condition(Read_depth), data = data.frame(sample_data(r_rms_seed_phyloseq)),distance = "bray", dfun = vegdist) 

 ----- summary -----
Partitioning of squared Bray distance:
              Inertia Proportion
Total          46.838     1.0000
Conditioned     6.842     0.1461
Constrained    11.041     0.2357
Unconstrained  28.956     0.6182

----  anova results -----
Df SumOfSqs      F Pr(>F)    
Model      4   11.041 15.062  0.001 ***
Residual 158   28.956
```

### After remove all the seed samples: non_seed samples

* **Plot impact**

```
          Df SumOfSqs     R2      F Pr(>F)    
Plot      14    5.501 0.1488 1.6732  0.001 ***
Residual 134   31.470 0.8512                  
Total    148   36.972 1.0000
```
* **Read_Depth**

```
            Df SumOfSqs      R2      F Pr(>F)    
Read_depth   1    6.630 0.17933 32.121  0.001 ***
Residual   147   30.342 0.82067                  
Total      148   36.972 1.00000
```
* **Compartment impact**

```
              Df SumOfSqs      R2      F Pr(>F)    
Compartment   3   10.938 0.29585 20.307  0.001 ***
Residual    145   26.034 0.70415                  
Total       148   36.972 1.00000
```

* **Compartment vs Read_Depth**

```
             Df SumOfSqs      R2      F Pr(>F)    
Read_depth    1    0.407 0.01100 2.2861  0.012 *  
Compartment   3    4.715 0.12753 8.8310  0.001 ***
Residual    144   25.627 0.69315                  
Total       148   36.972 1.00000
```

* **Read_depth vs Plot vs Compartment vs Time** -- Plot is a significant factor

```
adonis2(formula = t(otu_table(r_non_Seed_up)) ~ Read_depth + Compartment + Plot, data = data.frame(sample_data(r_non_Seed_up)), permutations = 999, by = "margin")
             Df SumOfSqs      R2      F Pr(>F)    
Read_depth    1    0.425 0.01149 2.7448  0.003 ** 
Compartment   3    4.511 0.12200 9.7149  0.001 ***
Plot         14    5.507 0.14895 2.5416  0.001 ***
Residual    130   20.120 0.54419                  
Total       148   36.972 1.00000
```
* **Read_depth vs  Compartment vs Time vs Treatment**

```
adonis2(formula = t(otu_table(r_non_Seed_up)) ~ Read_depth + Compartment + Treatment + Time, data = data.frame(sample_data(r_non_Seed_up)), permutations = 999, by = "margin")
             Df SumOfSqs      R2       F Pr(>F)    
Read_depth    1    0.234 0.00634  1.4257  0.103    
Compartment   2    4.060 0.10983 12.3449  0.001 ***
Treatment     2    0.931 0.02518  2.8306  0.001 ***
Time          2    1.675 0.04530  5.0916  0.001 ***
Residual    140   23.024 0.62275                   
Total       148   36.972 1.00000
```

### Subset to bulk+rhizosphere+endosphere samples, in short of BRE

* **Plot impacts**

```
adonis2(formula = t(otu_table(r_BRE_up)) ~ Plot, data = data.frame(sample_data(r_BRE_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)    
Plot      14    5.498 0.16102 1.6313  0.001 ***
Residual 119   28.648 0.83898                  
Total    133   34.146 1.00000
```

* **Read_depth**

```
adonis2(formula = t(otu_table(r_BRE_up)) ~ Read_depth, data = data.frame(sample_data(r_BRE_up)), permutations = 999, by = "margin")
            Df SumOfSqs      R2      F Pr(>F)    
Read_depth   1    6.140 0.17983 28.942  0.001 ***
Residual   132   28.006 0.82017                  
Total      133   34.146 1.00000
```

* **Compartment**

```
adonis2(formula = t(otu_table(r_BRE_up)) ~ Compartment, data = data.frame(sample_data(r_BRE_up)), permutations = 999, by = "margin")
             Df SumOfSqs      R2      F Pr(>F)    
Compartment   2    9.413 0.27567 24.929  0.001 ***
Residual    131   24.733 0.72433                  
Total       133   34.146 1.00000 
```

* **Time impact**

```
adonis2(formula = t(otu_table(r_BRE_up)) ~ Time, data = data.frame(sample_data(r_BRE_up)), permutations = 999, by = "margin")
          Df SumOfSqs    R2      F Pr(>F)    
Time       2    1.844 0.054 3.7389  0.001 ***
Residual 131   32.302 0.946                  
Total    133   34.146 1.000
```

* **Read_depth vs Compartment vs Time vs Plot**

```
adonis2(formula = t(otu_table(r_BRE_up)) ~ Read_depth + Compartment + Time + Plot, data = data.frame(sample_data(r_BRE_up)), permutations = 999, by = "margin")
             Df SumOfSqs      R2       F Pr(>F)    
Read_depth    1    0.326 0.00955  2.1751  0.014 *  
Compartment   2    3.570 0.10454 11.9073  0.001 ***
Time          2    1.597 0.04676  5.3266  0.001 ***
Plot         14    5.496 0.16095  2.6190  0.001 ***
Residual    114   17.088 0.50043                   
Total       133   34.146 1.00000
```
* **Read_depth vs Compartment vs Time vs Treatment**

```
adonis2(formula = t(otu_table(r_BRE_up)) ~ Read_depth + Compartment + Time + Treatment, data = data.frame(sample_data(r_BRE_up)), permutations = 999, by = "margin")
             Df SumOfSqs      R2       F Pr(>F)    
Read_depth    1    0.308 0.00902  1.7948  0.029 *  
Compartment   2    3.738 0.10947 10.8873  0.001 ***
Time          2    1.614 0.04728  4.7023  0.001 ***
Treatment     2    0.954 0.02793  2.7781  0.001 ***
Residual    126   21.630 0.63345                   
Total       133   34.146 1.00000
```

### subset to Bulk -- before rarefaction and after rarefaction, there are not pattern-changing difference


* --- Plot ----
```
adonis2(formula = t(otu_table(rf_bulk_up)) ~ Plot, data = data.frame(sample_data(rf_bulk_up)), permutations = 999, by = "margin")
         Df SumOfSqs      R2      F Pr(>F)   
Plot     14   2.4372 0.38695 1.3525  0.004 **
Residual 30   3.8613 0.61305                 
Total    44   6.2985 1.00000 
```

* ---- Read_depth ---

```
adonis2(formula = t(otu_table(rf_bulk_up)) ~ Read_depth, data = data.frame(sample_data(rf_bulk_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)   
Read_depth  1   0.3712 0.05893 2.6928  0.002 **
Residual   43   5.9273 0.94107                 
Total      44   6.2985 1.00000
```
* ---- Time -----

```
adonis2(formula = t(otu_table(rf_bulk_up)) ~ Time, data = data.frame(sample_data(rf_bulk_up)), permutations = 999, by = "margin")
         Df SumOfSqs      R2      F Pr(>F)    
Time      2   1.3968 0.22178 5.9845  0.001 ***
Residual 42   4.9016 0.77822                  
Total    44   6.2985 1.00000
```
* ---- Treatment ----

```
adonis2(formula = t(otu_table(rf_bulk_up)) ~ Treatment, data = data.frame(sample_data(rf_bulk_up)), permutations = 999, by = "margin")
          Df SumOfSqs     R2      F Pr(>F)
Treatment  2   0.3577 0.0568 1.2646  0.122
Residual  42   5.9407 0.9432              
Total     44   6.2985 1.0000
```
* ----Read_depth vs Time vs Plot

```
adonis2(formula = t(otu_table(rf_bulk_up)) ~ Read_depth + Plot + Time, data = data.frame(sample_data(rf_bulk_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)    
Read_depth  1   0.1217 0.01932 1.4025  0.106    
Plot       14   2.3727 0.37671 1.9532  0.001 ***
Time        2   1.1127 0.17666 6.4117  0.001 ***
Residual   27   2.3427 0.37195                  
Total      44   6.2985 1.00000
```

* ---- Read_depth vs Time vs Treatment

```
adonis2(formula = t(otu_table(rf_bulk_up)) ~ Read_depth + Treatment + Time, data = data.frame(sample_data(rf_bulk_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)    
Read_depth  1   0.1665 0.02644 1.4835  0.070 .  
Treatment   2   0.3380 0.05367 1.5059  0.038 *  
Time        2   1.1956 0.18982 5.3259  0.001 ***
Residual   39   4.3774 0.69499                  
Total      44   6.2985 1.00000
```

* ---- Treatment impact on Bulk_1wk, Bulk_3wk and Bulk_4wk

```
>>> Bulk_1wk --- insignificant

adonis2(formula = t(otu_table(r_rf_bulk_1wk_up)) ~ Treatment, data = data.frame(sample_data(r_rf_bulk_1wk_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)
Treatment  2  0.23437 0.13121 0.9061  0.728
Residual  12  1.55191 0.86879              
Total     14  1.78628 1.00000

>>> Bulk_3wk --- insignificant

adonis2(formula = t(otu_table(rf_bulk_3wk_up)) ~ Treatment, data = data.frame(sample_data(rf_bulk_3wk_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)
Treatment  2  0.25781 0.14583 1.0243  0.379
Residual  12  1.51014 0.85417              
Total     14  1.76795 1.00000

>>> Bulk_4wk ---insignificant

adonis2(formula = t(otu_table(rf_bulk_4wk_up)) ~ Treatment, data = data.frame(sample_data(rf_bulk_4wk_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)
Treatment  2  0.21858 0.14487 1.0165  0.419
Residual  12  1.29021 0.85513              
Total     14  1.50880 1.00000
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


