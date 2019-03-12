## Bacterial community analysis were conducted in R

### Overal impacts

1. **Plot variations** - Subset phyloseq object to only include non-seed samples: PERMANOVA results indicate an non-significant plots impacts on overal bacterial community composition.

```
Permutation test for adonis under NA model
Marginal effects of terms
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(r_non_Seed_22017_up)) ~ Plot, data = data.frame(sample_data(r_non_Seed_22017_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)
Plot      14    3.079 0.08886 0.9195  0.643
Residual 132   31.569 0.91114              
Total    146   34.648 1.00000
```

2. **Read_depth** impacts

```
Permutation test for adonis under NA model
Marginal effects of terms
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(r_Seed_22017_up)) ~ Read_depth, data = data.frame(sample_data(r_Seed_22017_up)), permutations = 999, by = "margin")
            Df SumOfSqs      R2      F Pr(>F)      
Read_depth   1    0.899 0.02098 3.3649  0.004 **
Residual   157   41.950 0.97902                 
Total      158   42.849 1.00000  
```
3. Overal impacts of **Read_depth, Compartment, Time and Treatments**

```
Permutation test for adonis under reduced model
Marginal effects of terms
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(r_Seed_22017_up)) ~ Read_depth + Compartment + Treatment + Time, data = data.frame(sample_data(r_Seed_22017_up)), permutations = 999, by = "margin")
             Df SumOfSqs      R2       F Pr(>F)    
Read_depth    1    0.113 0.00264  1.1577  0.300    
Compartment   2   14.684 0.34269 75.1768  0.001 ***
Treatment     2    0.404 0.00944  2.0707  0.026 *  
Time          2    3.151 0.07355 16.1345  0.001 ***
Residual    149   14.552 0.33961                   
Total       158   42.849 1.00000    
```


**Considering the significant impacts of compartment, soybean development (Time) and seed fungicide treatment impacts will be analyzed after subset phyloseq to corresponding compartment**

### Compare the degree of compartment impacts along soybean development 

1. **1 week old soybeans**

```
> adonis2(t(otu_table(wk1_up ))~Compartment,data=data.frame(sample_data(wk1_up )),permutations=999,by="margin")
Permutation test for adonis under NA model
Marginal effects of terms
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(wk1_up)) ~ Compartment, data = data.frame(sample_data(wk1_up)), permutations = 999, by = "margin")
            Df SumOfSqs      R2      F Pr(>F)    
Compartment  2   5.9422 0.61282 32.447  0.001 ***
Residual    41   3.7543 0.38718                  
Total       43   9.6965 1.00000                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

2. **3 weeks old soybeans**

```
Permutation test for adonis under NA model
Marginal effects of terms
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(wk3_up)) ~ Compartment, data = data.frame(sample_data(wk3_up)), permutations = 999, by = "margin")
            Df SumOfSqs      R2      F Pr(>F)    
Compartment  2   6.8324 0.71468 52.602  0.001 ***
Residual    42   2.7276 0.28532                  
Total       44   9.5600 1.00000                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

3. **4 weeks old soybeans**

```
adonis2(t(otu_table(wk4_up ))~Compartment,data=data.frame(sample_data(wk4_up )),permutations=999,by="margin")
Permutation test for adonis under NA model
Marginal effects of terms
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(wk4_up)) ~ Compartment, data = data.frame(sample_data(wk4_up)), permutations = 999, by = "margin")
            Df SumOfSqs      R2     F Pr(>F)    
Compartment  2   6.3334 0.67224 41.02  0.001 ***
Residual    40   3.0879 0.32776                 
Total       42   9.4213 1.00000                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

## Treatment and soybean development impacts difference among compartments

1. **Bulk samples**

```
adonis2(t(otu_table(r_Bulk_22017_up))~Treatment+Time,data=data.frame(sample_data(r_Bulk_22017_up)),permutations=999,by="margin") # Time is significant, R2=0.12 and p=0.001
Permutation test for adonis under reduced model
Marginal effects of terms
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(r_Bulk_22017_up)) ~ Treatment + Time, data = data.frame(sample_data(r_Bulk_22017_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)    
Treatment  2   0.2691 0.07207 1.7497  0.005 ** 
Time       2   0.4649 0.12451 3.0228  0.001 ***
Residual  39   2.9992 0.80323                  
Total     43   3.7340 1.00000                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Interaction effects
adonis2(t(otu_table(r_Bulk_22017_up))~Treatment*Time,data=data.frame(sample_data(r_Bulk_22017_up)),permutations=999,by="margin") # no interaction impacts
Permutation test for adonis under reduced model
Marginal effects of terms
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(r_Bulk_22017_up)) ~ Treatment * Time, data = data.frame(sample_data(r_Bulk_22017_up)), permutations = 999, by = "margin")
               Df SumOfSqs      R2      F Pr(>F)
Treatment:Time  4   0.2040 0.05464 0.6387      1
Residual       35   2.7952 0.74859              
Total          43   3.7340 1.00000    
```
2. **Rhizosphere**

```
adonis2(t(otu_table(Rhi_22017_up))~Treatment+Time,data=data.frame(sample_data(Rhi_22017_up)),permutations=999,by="margin") # Time is significant, R2=0.42 and p=0.001
Permutation test for adonis under reduced model
Marginal effects of terms
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(Rhi_22017_up)) ~ Treatment + Time, data = data.frame(sample_data(Rhi_22017_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2       F Pr(>F)    
Treatment  2   0.2985 0.04465  1.6478  0.102    
Time       2   2.8398 0.42473 15.6741  0.001 ***
Residual  39   3.5330 0.52841                   
Total     43   6.6861 1.00000                   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```



3. **Endosphere**







