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



