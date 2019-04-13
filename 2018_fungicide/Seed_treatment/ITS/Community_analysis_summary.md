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

* **Compartment**

```
adonis2(formula = t(otu_table(r_rms_seed_phyloseq)) ~ Read_depth, data = data.frame(sample_data(r_rms_seed_phyloseq)), permutations = 999, by = "margin")
            Df SumOfSqs      R2      F Pr(>F)    
Read_depth   1    6.841 0.14874 28.307  0.001 ***
Residual   162   39.152 0.85126                  
Total      163   45.993 1.00000  
```
* **Read_depth vs Compartment**

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

#### non_seed samples

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

* **Compartment vs Read_depth vs Plot**

```
adonis2(formula = t(otu_table(r_BRE_up)) ~ Compartment + Read_depth + Plot, data = data.frame(sample_data(r_BRE_up)), permutations = 999, by = "margin")
             Df SumOfSqs      R2       F Pr(>F)    
Compartment   2    3.635 0.10647 11.2852  0.001 ***
Read_depth    1    0.565 0.01654  3.5070  0.001 ***
Plot         14    5.510 0.16137  2.4435  0.001 ***
Residual    116   18.684 0.54719                   
Total       133   34.146 1.00000
```
* **Read_depth vs Plot vs Compartment vs Time vs Treatment**

```
adonis2(formula = t(otu_table(r_BRE_up)) ~ Read_depth + Plot + Compartment + Time + Treatment, data = data.frame(sample_data(r_BRE_up)), permutations = 999, by = "margin")
             Df SumOfSqs      R2       F Pr(>F)    
Read_depth    1    0.326 0.00955  2.1751  0.014 *  
Plot         12    4.542 0.13302  2.5252  0.001 ***
Compartment   2    3.570 0.10454 11.9073  0.001 ***
Time          2    1.597 0.04676  5.3266  0.001 ***
Treatment     0    0.000 0.00000    -Inf           
Residual    114   17.088 0.50043                   
Total       133   34.146 1.00000
```

#### subset to SBR (soil and bulk and rhizosphere samples)

* **Read_depth before rarefaction**

```
adonis2(formula = t(otu_table(SBR_up)) ~ Read_depth, data = data.frame(sample_data(SBR_up)), permutations = 999, by = "margin")
            Df SumOfSqs      R2      F Pr(>F)    
Read_depth   1   3.0231 0.14749 17.647  0.001 ***
Residual   102  17.4738 0.85251                  
Total      103  20.4969 1.00000
```
* **Read_depth after rarefaction**

```
adonis2(formula = t(otu_table(r_rf_SBR_up)) ~ Read_depth, data = data.frame(sample_data(r_rf_SBR_up)), permutations = 999, by = "margin")
            Df SumOfSqs    R2      F Pr(>F)    
Read_depth   1   1.4586 0.078 8.6288  0.001 ***
Residual   102  17.2419 0.922                  
Total      103  18.7005 1.000
```



