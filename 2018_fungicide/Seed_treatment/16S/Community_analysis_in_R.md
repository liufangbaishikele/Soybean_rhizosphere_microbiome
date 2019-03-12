## Bacterial community analysis were conducted in R

**Highligts **

1. Read_depth impacts

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
2. Overal evaluation of compartment, Time and Treatment impacts

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
3. Plot variations - Subset phyloseq object to only include non-seed samples: PERMANOVA results indicate an non-significant plots impacts on overal bacterial community composition.

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

4. 
