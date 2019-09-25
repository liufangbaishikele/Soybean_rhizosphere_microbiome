## Due to the very limited ITS2 sequencing depth for endosphere and seed samples, we designed PNA and resequenced those with low sequencing depth.

The below results are based on combined sequence data.

**Summary between normalization strategy**

1. No normalization vs DEseq2 normalization --- compartment impacts

  *No rarefaction*

  ```
  Permutation test for adonis under reduced model
  Marginal effects of terms
  Permutation: free
  Number of permutations: 999

  adonis2(formula = t(otu_table(r_rms_seed_phyloseq)) ~ Read_depth + Compartment, data = data.frame(sample_data(r_rms_seed_phyloseq)), permutations = 999, by = "margin")
               Df SumOfSqs      R2       F Pr(>F)    
  Read_depth    1    0.260 0.00587  1.6088  0.073 .  
  Compartment   4   15.010 0.33837 23.1983  0.001 ***
  Residual    158   25.557 0.57615                   
  Total       163   44.358 1.00000                   
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  ```
  *DESeq2 based rarefaction*

  ```
  Permutation test for adonis under reduced model
  Marginal effects of terms
  Plots: paste(data.frame(sample_data(r_norm_seed_phyloseq))$Plot, data.frame(sample_data(r_norm_seed_phyloseq))$Time, , plot permutation: none
  Permutation: free
  Number of permutations: 999

  adonis2(formula = t(otu_table(r_norm_seed_phyloseq)) ~ Compartment + Read_depth, data = data.frame(sample_data(r_norm_seed_phyloseq)), permutations = perm, by = "margin")
  Model: adonis0(formula = lhs ~ Compartment + Read_depth, data = data, method = method)
               Df SumOfSqs      R2      F Pr(>F)    
  Compartment   4   15.066 0.33885 23.261  0.001 ***
  Read_depth    1    0.259 0.00582  1.598  0.058 .  
  Residual    158   25.583 0.57540                  
  Total       163   44.461 1.00000
  ```
2. No rarefaction versus minimum depth based rarefaction

*No rarefaction*

```
Permutation test for adonis under reduced model
Marginal effects of terms
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(r_rms_Rhi_up)) ~ Read_depth + Time + Plot, data = data.frame(sample_data(r_rms_Rhi_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)    
Read_depth  1   0.1177 0.01570 1.0813  0.355    
Time        2   1.2127 0.16174 5.5714  0.001 ***
Plot       14   3.0498 0.40675 2.0015  0.001 ***
Residual   26   2.8298 0.37741                  
Total      43   7.4979 1.00000                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

*Minimum read depth based rarefaction*

```
Permutation test for adonis under reduced model
Marginal effects of terms
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(r_rf_rms_Rhi_up)) ~ Read_depth + Time + Plot, data = data.frame(sample_data(r_rf_rms_Rhi_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)    
Read_depth  1   0.1165 0.01544 1.0596  0.374    
Time        2   1.2187 0.16158 5.5436  0.001 ***
Plot       14   3.0620 0.40598 1.9898  0.001 ***
Residual   26   2.8579 0.37891                  
Total      43   7.5424 1.00000                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

**Conclusion**: Rarefaction / normalization trategy did not change PERMANOVA results. 


