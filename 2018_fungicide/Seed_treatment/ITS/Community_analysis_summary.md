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

