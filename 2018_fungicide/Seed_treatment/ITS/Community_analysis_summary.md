##                                      Community analysis summary


----
Fang Liu
04/13/2019


**Overal description**

This is the first run of ITS sequencing targeted ITS2 region. Due to the large propotion of soybean derived ITS, endosphere samples and seed samples were left with less informative sequences. Anyway, here are the summary of this run results.

### PERMANOVA results

#### All samples without rarefaction

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

#### After remove all the seed samples: non_seed samples

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

#### Subset to bulk+rhizosphere+endosphere samples, in short of BRE

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

#### subset to Bulk -- before rarefaction and after rarefaction, there are not pattern-changing difference


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





