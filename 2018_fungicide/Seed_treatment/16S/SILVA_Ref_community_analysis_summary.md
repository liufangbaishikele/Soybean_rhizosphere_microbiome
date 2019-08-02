## Soybean root associated bacterial community dynamics in response to plant development

author | Date
-----| -----
Fang Liu | 08/02/2019

### Introduction

For this experiment, soybean bulk, rhizosphere and endosphere microbes were characterized based on 16S rRNA V3-V4 sequencing data. MOTHUR was used to process the sequencing data to generate community abundance dataframe. R packages and Python based softwares were used to analyze bacteria abundance data. [R Code couls be found here](https://github.com/liufangbaishikele/Soybean_rhizosphere_microbiome/blob/master/2018_fungicide/Seed_treatment/16S/2018_seed_fungicide_16S_13021.Rmd)

## Results summary


### Before rarefaction at 13021**

           * Summary of the otu_table, tax_table and the sample_data
           ```
           phyloseq-class experiment-level object
           otu_table()   OTU Table:         [ 144821 taxa and 164 samples ]
           sample_data() Sample Data:       [ 164 samples by 7 sample variables ]
           tax_table()   Taxonomy Table:    [ 144821 taxa by 6 taxonomic ranks ]
           ```
           * In total, it has 13079622 reads after mothur pipeline before doing any rarefaction

### r_Seed_13021_up

   1. Summary of the otu_table, tax_table and sample_data
           
           ```
           phyloseq-class experiment-level object
           otu_table()   OTU Table:         [ 21355 taxa and 149 samples ]
           sample_data() Sample Data:       [ 149 samples by 7 sample variables ]
           tax_table()   Taxonomy Table:    [ 21355 taxa by 6 taxonomic ranks ]
           ```
   2. Total read - 1913852

   3. PERMANOVA results
           
        * Read-depth

            ```
            adonis2(formula = t(otu_table(r_Seed_13021_up)) ~ Read_depth, data = data.frame(sample_data(r_Seed_13021_up)), permutations = 999, by = "margin")
                        Df SumOfSqs      R2      F Pr(>F)  
            Read_depth   1    0.478 0.01323 1.9714  0.067 .
            Residual   147   35.644 0.98677                
            Total      148   36.122 1.00000
            ```

### r_Soil_up --check for field plot design

   *Plot design impact check up*

We have 15 plots that randomly applied CM, EE or CT fungicide treatment. So, we want to test if there were any unintended clustering of plots that accidently overlap with the fungicide treatment. For this purpose, we analized the bacteria community difference between samples before we do any fungicide treatment.

   * Treatment impact
   
           ```
           adonis2(formula = t(otu_table(r_Soil_up)) ~ Treatment, data = data.frame(sample_data(r_Soil_up)), permutations = 999, by = "margin")
                     Df SumOfSqs      R2      F Pr(>F)
           Treatment  2  0.17841 0.14635 1.0286  0.377
           Residual  12  1.04069 0.85365              
           Total     14  1.21909 1.00000
           ```

### r_SB_CT_up - more exploration of the data

   * Time impact
           ```
                      adonis2(formula = t(otu_table(r_SB_CT_up)) ~ Time, data = data.frame(sample_data(r_SB_CT_up)), permutations = perm, by = "margin")
                    Df SumOfSqs      R2      F Pr(>F)    
           Time      3  0.42175 0.23849 1.6703  0.001 
           Residual 16  1.34668 0.76151                  
           Total    19  1.76842 1.00000
           ```
   * Read_depth impact
           
           ```
                      adonis2(formula = t(otu_table(r_SB_CT_up)) ~ Read_depth, data = data.frame(sample_data(r_SB_CT_up)), permutations = 999, by = "margin")
                      Df SumOfSqs      R2      F Pr(>F)  
           Read_depth  1  0.13545 0.07659 1.4931  0.042 *
           Residual   18  1.63297 0.92341                
           Total      19  1.76842 1.00000
           ```
### r_BRE_up - including bulk, rhizosphere and endosphere samples

   * Summary of the otu_table, tax_table and metadata

           ```
           phyloseq-class experiment-level object
           otu_table()   OTU Table:         [ 19681 taxa and 134 samples ]
           sample_data() Sample Data:       [ 134 samples by 7 sample variables ]
           tax_table()   Taxonomy Table:    [ 19681 taxa by 6 taxonomic ranks ]
           ```
* **PERMANOVA** 

   * ---- Read_depth - insignificant

           ```
           adonis2(formula = t(otu_table(r_BRE_up)) ~ Read_depth, data = data.frame(sample_data(r_BRE_up)), permutations = 999, by = "margin")
                       Df SumOfSqs      R2      F Pr(>F)
           Read_depth   1    0.373 0.01121 1.4958  0.144
           Residual   132   32.886 0.98879              
           Total      133   33.258 1.00000
           ```


   * ---- Compartment impact

           ```
           Permutation test for adonis under NA model
           Marginal effects of terms
           Plots: paste(data.frame(sample_data(r_BRE_up))$Time, data.frame(sample_data(r_BRE_up))$Plot, , plot permutation: none
           Permutation: free
           Number of permutations: 999

           adonis2(formula = t(otu_table(r_BRE_up)) ~ Compartment, data = data.frame(sample_data(r_BRE_up)), permutations = perm, by = "margin")
           Model: adonis0(formula = lhs ~ Compartment, data = data, method = method)
                        Df SumOfSqs      R2      F Pr(>F)    
           Compartment   2   15.215 0.45748 55.233  0.001 ***
           Residual    131   18.043 0.54252                  
           Total       133   33.258 1.00000                  
           ---
           Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
           ```

   * ---- Time impact

           ```
           Permutation test for adonis under NA model
           Marginal effects of terms
           Plots: paste(data.frame(sample_data(r_BRE_up))$Compartment, data.frame(sample_data(r_BRE_up))$Plot, , plot permutation: none
           Permutation: free
           Number of permutations: 999

           adonis2(formula = t(otu_table(r_BRE_up)) ~ Time, data = data.frame(sample_data(r_BRE_up)), permutations = perm, by = "margin")
           Model: adonis0(formula = lhs ~ Time, data = data, method = method)
                     Df SumOfSqs      R2      F Pr(>F)    
           Time       2    3.183 0.09572 6.9332  0.001 ***
           Residual 131   30.075 0.90428                  
           Total    133   33.258 1.00000                  
           ---
           Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
           ```

   * ---- Plot impact vs Treatment impact

           ```
           1.  Treatment
           Permutation test for adonis under NA model
           Marginal effects of terms
           Plots: paste(data.frame(sample_data(r_BRE_up))$Compartment, data.frame(sample_data(r_BRE_up))$Time, , plot permutation: none
           Permutation: free
           Number of permutations: 999

           adonis2(formula = t(otu_table(r_BRE_up)) ~ Treatment, data = data.frame(sample_data(r_BRE_up)), permutations = perm, by = "margin")
           Model: adonis0(formula = lhs ~ Treatment, data = data, method = method)
                      Df SumOfSqs      R2      F Pr(>F)    
           Treatment   2    0.411 0.01235 0.8188  0.001 ***
           Residual  131   32.848 0.98765                  
           Total     133   33.258 1.00000

           2. Plot 
           Permutation test for adonis under NA model
           Marginal effects of terms
           Plots: paste(data.frame(sample_data(r_BRE_up))$Compartment, data.frame(sample_data(r_BRE_up))$Time, , plot permutation: none
           Permutation: free
           Number of permutations: 999

           adonis2(formula = t(otu_table(r_BRE_up)) ~ Plot, data = data.frame(sample_data(r_BRE_up)), permutations = perm, by = "margin")
           Model: adonis0(formula = lhs ~ Plot, data = data, method = method)
                     Df SumOfSqs      R2      F Pr(>F)    
           Plot      14    2.880 0.08659 0.8058  0.001 ***
           Residual 119   30.378 0.91341                  
           Total    133   33.258 1.00000                  
           ---
           Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
           > 

           ```

   * ---- Read_Depth vs Compartment vs Time vs Plot

           ```
           Terms added sequentially (first to last)

                        Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
           Compartment   2    15.215  7.6075  73.138 0.45748  0.001 ***
           Plot         14     2.867  0.2048   1.969 0.08619  0.001 ***
           Time          2     3.199  1.5996  15.379 0.09619  0.001 ***
           Read_depth    1     0.120  0.1195   1.149 0.00359  0.291    
           Residuals   114    11.858  0.1040         0.35654           
           Total       133    33.258                 1.00000
           ```

           * ---- Read_depth vs Compartment vs Time vs Treatment

           ```
           Terms added sequentially (first to last)

                        Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
           Compartment   2    15.215  7.6075  66.937 0.45748  0.001 ***
           Treatment     2     0.407  0.2034   1.789 0.01223  0.057 .  
           Time          2     3.206  1.6030  14.104 0.09639  0.001 ***
           Read_depth    1     0.110  0.1104   0.972 0.00332  0.402    
           Residuals   126    14.320  0.1137         0.43057           
           Total       133    33.258                 1.00000
           ```

### Compare Time impacts between Bulk, Rhizosphere and Endosphere samples**

**Bulk Compartment**

           * --- Read_depth -- insignificant

           ```
           adonis2(formula = t(otu_table(r_Bulk_up)) ~ Read_depth, data = data.frame(sample_data(r_Bulk_up)), permutations = 999, by = "margin")
                      Df SumOfSqs      R2      F Pr(>F)
           Read_depth  1   0.1149 0.02677 1.1828   0.19
           Residual   43   4.1754 0.97323              
           Total      44   4.2903 1.00000
           ```
           * --- Time impact

           ```
           adonis2(formula = t(otu_table(r_Bulk_up)) ~ Time, data = data.frame(sample_data(r_Bulk_up)), permutations = perm, by = "margin")
                    Df SumOfSqs     R2     F Pr(>F)    
           Time      2   0.5144 0.1199 2.861  0.001 ***
           Residual 42   3.7759 0.8801                 
           Total    44   4.2903 1.0000
           ```
           * --- Plot impact

           ```
           adonis2(formula = t(otu_table(r_Bulk_up)) ~ Plot, data = data.frame(sample_data(r_Bulk_up)), permutations = perm, by = "margin")
                    Df SumOfSqs      R2      F Pr(>F)    
           Plot     14   2.0667 0.48172 1.9917  0.001 ***
           Residual 30   2.2235 0.51828                  
           Total    44   4.2903 1.00000
           ```

           * --- Treatment impact

           ```
           adonis2(formula = t(otu_table(r_Bulk_up)) ~ Treatment, data = data.frame(sample_data(r_Bulk_up)), permutations = 999, by = "margin")
                     Df SumOfSqs      R2      F Pr(>F)  
           Treatment  2   0.2661 0.06201 1.3884  0.038 *
           Residual  42   4.0242 0.93799                
           Total     44   4.2903 1.00000
           ```

           * --- Read_depth vs Time vs Plot

           ```
           Permutation test for adonis under reduced model
           Marginal effects of terms
           Permutation: free
           Number of permutations: 999

           adonis2(formula = t(otu_table(r_Bulk_up)) ~ Read_depth + Plot + Time, data = data.frame(sample_data(r_Bulk_up)), permutations = 999, by = "margin")
                      Df SumOfSqs      R2      F Pr(>F)    
           Read_depth  1   0.0552 0.01286 0.9007  0.577    
           Plot       14   2.0061 0.46760 2.3392  0.001 ***
           Time        2   0.4696 0.10945 3.8326  0.001 ***
           Residual   27   1.6540 0.38552                  
           Total      44   4.2903 1.00000                  
           ---
           Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
           ```
           * --- Read_depth vs Time vs Treatment

           ```
           adonis2(formula = t(otu_table(r_Bulk_up)) ~ Read_depth + Treatment + Time, data = data.frame(sample_data(r_Bulk_up)), permutations = 999, by = "margin")
                      Df SumOfSqs      R2      F Pr(>F)    
           Read_depth  1   0.0780 0.01819 0.8867  0.595    
           Treatment   2   0.2283 0.05322 1.2973  0.078 .  
           Time        2   0.5033 0.11732 2.8601  0.001 ***
           Residual   39   3.4318 0.79990                  
           Total      44   4.2903 1.00000                  
           ---
           Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
           ```
**Rhizosphere Compartment**

           * --- Read_depth impact
           ```
           adonis2(formula = t(otu_table(r_Rhi_up)) ~ Read_depth, data = data.frame(sample_data(r_Rhi_up)), permutations = 999, by = "margin")
                      Df SumOfSqs      R2      F Pr(>F)
           Read_depth  1   0.2444 0.03495 1.5209  0.126
           Residual   42   6.7492 0.96505              
           Total      43   6.9936 1.00000
           ```

           * --- Time
           ```
           adonis2(formula = t(otu_table(r_Rhi_up)) ~ Time, data = data.frame(sample_data(r_Rhi_up)), permutations = perm, by = "margin")
                    Df SumOfSqs      R2      F Pr(>F)    
           Time      2   2.8993 0.41457 14.517  0.001 ***
           Residual 41   4.0943 0.58543                  
           Total    43   6.9936 1.00000                  
           ---
           Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
           ```

           * --- Plot

           ```
           adonis2(formula = t(otu_table(r_Rhi_up)) ~ Plot, data = data.frame(sample_data(r_Rhi_up)), permutations = perm, by = "margin")
                    Df SumOfSqs      R2      F Pr(>F)    
           Plot     14   2.1208 0.30325 0.9016  0.001 ***
           Residual 29   4.8728 0.69675                  
           Total    43   6.9936 1.00000 
           ```

           * --- Treatment
           ```
           adonis2(formula = t(otu_table(r_Rhi_up)) ~ Treatment, data = data.frame(sample_data(r_Rhi_up)), permutations = perm, by = "margin")
                     Df SumOfSqs      R2     F Pr(>F)   
           Treatment  2   0.3209 0.04589 0.986  0.008 **
           Residual  41   6.6727 0.95411                
           Total     43   6.9936 1.00000                
           ---
           Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
           ```

           * --- Read_depth vs Time vs Plot

           ```
           adonis2(formula = t(otu_table(r_Rhi_up)) ~ Read_depth + Plot + Time, data = data.frame(sample_data(r_Rhi_up)), permutations = 999, by = "margin")
                      Df SumOfSqs      R2       F Pr(>F)    
           Read_depth  1   0.0565 0.00807  0.7386  0.614    
           Plot       14   1.9242 0.27514  1.7979  0.002 ** 
           Time        2   2.6118 0.37345 17.0819  0.001 ***
           Residual   26   1.9877 0.28421                   
           Total      43   6.9936 1.00000                   
           ---
           Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
           ```

           * --- Read_depth vs Time vs Treatment

           ```
           adonis2(formula = t(otu_table(r_Rhi_up)) ~ Read_depth + Treatment + Time, data = data.frame(sample_data(r_Rhi_up)), permutations = 999, by = "margin")
                      Df SumOfSqs      R2       F Pr(>F)    
           Read_depth  1   0.1315 0.01881  1.3668  0.178    
           Treatment   2   0.2549 0.03645  1.3245  0.177    
           Time        2   2.7381 0.39151 14.2262  0.001 ***
           Residual   38   3.6569 0.52289                   
           Total      43   6.9936 1.00000                   
           ---
           Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
           ```
**Endosphere Compartment**

          * --- Read_depth

          ```
          adonis2(formula = t(otu_table(r_Endo_up)) ~ Read_depth, data = data.frame(sample_data(r_Endo_up)), permutations = 999, by = "margin")
                     Df SumOfSqs      R2      F Pr(>F)
          Read_depth  1   0.1217 0.01845 0.8081  0.453
          Residual   43   6.4770 0.98155              
          Total      44   6.5987 1.00000
          ```

          * --- Time

          ```
          Permutation test for adonis under NA model
          Marginal effects of terms
          Plots: data.frame(sample_data(r_Endo_up))$Plot, plot permutation: none
          Permutation: free
          Number of permutations: 999

          adonis2(formula = t(otu_table(r_Endo_up)) ~ Time, data = data.frame(sample_data(r_Endo_up)), permutations = perm, by = "margin")
                   Df SumOfSqs      R2      F Pr(>F)    
          Time      2   3.9672 0.60121 31.659  0.001 ***
          Residual 42   2.6315 0.39879                  
          Total    44   6.5987 1.00000                  
          ---
          Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
          ```

          * --- Plot

          ```
          Permutation test for adonis under NA model
          Marginal effects of terms
          Plots: data.frame(sample_data(r_Endo_up))$Time, plot permutation: none
          Permutation: free
          Number of permutations: 999

          adonis2(formula = t(otu_table(r_Endo_up)) ~ Plot, data = data.frame(sample_data(r_Endo_up)), permutations = perm, by = "margin")
                   Df SumOfSqs      R2      F Pr(>F)    
          Plot     14   1.0839 0.16426 0.4212  0.001 ***
          Residual 30   5.5148 0.83574                  
          Total    44   6.5987 1.00000                  
          ---
          Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
          ```

          * --- Treatment

          ```
          Permutation test for adonis under NA model
          Marginal effects of terms
          Plots: data.frame(sample_data(r_Endo_up))$Time, plot permutation: none
          Permutation: free
          Number of permutations: 999

          adonis2(formula = t(otu_table(r_Endo_up)) ~ Treatment, data = data.frame(sample_data(r_Endo_up)), permutations = perm, by = "margin")
                    Df SumOfSqs      R2      F Pr(>F)  
          Treatment  2   0.1590 0.02409 0.5185  0.095 .
          Residual  42   6.4397 0.97591                
          Total     44   6.5987 1.00000                
          ---
          Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
          ```

          * --- Read_depth vs Time vs Plot

          ```
          Permutation test for adonis under reduced model
          Marginal effects of terms
          Permutation: free
          Number of permutations: 999

          adonis2(formula = t(otu_table(r_Endo_up)) ~ Read_depth + Time + Plot, data = data.frame(sample_data(r_Endo_up)), permutations = 999, by = "margin")
                     Df SumOfSqs      R2       F Pr(>F)    
          Read_depth  1   0.1292 0.01957  2.4585  0.062 .  
          Time        2   3.8746 0.58718 36.8758  0.001 ***
          Plot       14   1.1595 0.17572  1.5765  0.043 *  
          Residual   27   1.4185 0.21496                   
          Total      44   6.5987 1.00000                   
          ---
          Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
          ```

          * --- Read_depth vs Time vs Treatment

          ```
          Permutation test for adonis under reduced model
          Marginal effects of terms
          Permutation: free
          Number of permutations: 999

          adonis2(formula = t(otu_table(r_Endo_up)) ~ Read_depth + Time + Treatment, data = data.frame(sample_data(r_Endo_up)), permutations = 999, by = "margin")
                     Df SumOfSqs      R2       F Pr(>F)    
          Read_depth  1   0.0646 0.00979  1.0462  0.318    
          Time        2   3.8962 0.59045 31.5524  0.001 ***
          Treatment   2   0.1701 0.02577  1.3772  0.195    
          Residual   39   2.4079 0.36491                   
          Total      44   6.5987 1.00000                   
          ---
          Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
          ```

### Comparison the compartment impact between different time points

**>>>>Week1>>>>**

           ```
           phyloseq-class experiment-level object
           otu_table()   OTU Table:         [ 11316 taxa and 44 samples ]
           sample_data() Sample Data:       [ 44 samples by 7 sample variables ]
           tax_table()   Taxonomy Table:    [ 11316 taxa by 6 taxonomic ranks ]
           ```

           * Read_depth

           ```
           adonis2(formula = t(otu_table(r_wk1_up)) ~ Read_depth, data = data.frame(sample_data(r_wk1_up)), permutations = 999, by = "margin")
                      Df SumOfSqs      R2      F Pr(>F)
           Read_depth  1   0.1056 0.01061 0.4503  0.885
           Residual   42   9.8447 0.98939              
           Total      43   9.9503 1.00000
           ```

           * Compartment

           ```
           Permutation test for adonis under NA model
           Marginal effects of terms
           Plots: data.frame(sample_data(r_wk1_up))$Plot, plot permutation: none
           Permutation: free
           Number of permutations: 999

           adonis2(formula = t(otu_table(r_wk1_up)) ~ Compartment, data = data.frame(sample_data(r_wk1_up)), permutations = perm, by = "margin")
                       Df SumOfSqs      R2      F Pr(>F)    
           Compartment  2   5.9584 0.59882 30.599  0.001 ***
           Residual    41   3.9919 0.40118                  
           Total       43   9.9503 1.00000                  
           ---
           Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
           ```
           * Plot

           ```
           Permutation test for adonis under NA model
           Marginal effects of terms
           Plots: data.frame(sample_data(r_wk1_up))$Compartment, plot permutation: none
           Permutation: free
           Number of permutations: 999

           adonis2(formula = t(otu_table(r_wk1_up)) ~ Plot, data = data.frame(sample_data(r_wk1_up)), permutations = perm, by = "margin")
                    Df SumOfSqs      R2      F Pr(>F)    
           Plot     14   1.9146 0.19241 0.4935  0.001 ***
           Residual 29   8.0357 0.80759                  
           Total    43   9.9503 1.00000                  
           ---
           Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
           ```

           * Treatment

           ```
           Permutation test for adonis under NA model
           Marginal effects of terms
           Plots: data.frame(sample_data(r_wk1_up))$Compartment, plot permutation: none
           Permutation: free
           Number of permutations: 999

           adonis2(formula = t(otu_table(r_wk1_up)) ~ Treatment, data = data.frame(sample_data(r_wk1_up)), permutations = perm, by = "margin")
                     Df SumOfSqs      R2      F Pr(>F)
           Treatment  2   0.2325 0.02337 0.4906  0.161
           Residual  41   9.7177 0.97663              
           Total     43   9.9503 1.00000
           ```

           * Read_depth vs Compartment vs Plot

           ```
           adonis2(formula = t(otu_table(r_wk1_up)) ~ Compartment + Read_depth + Plot, data = data.frame(sample_data(r_wk1_up)), permutations = 999, by = "margin")
                       Df SumOfSqs      R2       F Pr(>F)    
           Compartment  2   5.8181 0.58471 36.1110  0.001 ***
           Read_depth   1   0.0626 0.00629  0.7775  0.529    
           Plot        14   1.7920 0.18010  1.5889  0.019 *  
           Residual    26   2.0945 0.21050                   
           Total       43   9.9503 1.00000                   
           ---
           Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
           ```

           * Read_depth vs Compartment vs Treatment

           ```
           adonis2(formula = t(otu_table(r_wk1_up)) ~ Compartment + Read_depth + Treatment, data = data.frame(sample_data(r_wk1_up)), permutations = 999, by = "margin")
                       Df SumOfSqs      R2       F Pr(>F)    
           Compartment  2   5.9410 0.59707 30.7394  0.001 ***
           Read_depth   1   0.1011 0.01016  1.0462  0.348    
           Treatment    2   0.2144 0.02155  1.1095  0.342    
           Residual    38   3.6721 0.36905                   
           Total       43   9.9503 1.00000
           ```

**>>>>Week3>>>>**

           ```
           phyloseq-class experiment-level object
           otu_table()   OTU Table:         [ 11403 taxa and 45 samples ]
           sample_data() Sample Data:       [ 45 samples by 7 sample variables ]
           tax_table()   Taxonomy Table:    [ 11403 taxa by 6 taxonomic ranks ]
           ```
           * ---- Read_depth -- this is cause by one specific sample with super high read depth - 105E_3wk sample

           ```
           adonis2(formula = t(otu_table(r_wk3_up)) ~ Read_depth, data = data.frame(sample_data(r_wk3_up)), permutations = 999, by = "margin")
                      Df SumOfSqs      R2      F Pr(>F)    
           Read_depth  1   1.2177 0.12417 6.0962  0.001 ***
           Residual   43   8.5889 0.87583                  
           Total      44   9.8065 1.00000
           ```

           * ---- Compartment

           ```
           Permutation test for adonis under NA model
           Marginal effects of terms
           Plots: data.frame(sample_data(r_wk3_up))$Plot, plot permutation: none
           Permutation: free
           Number of permutations: 999

           adonis2(formula = t(otu_table(r_wk3_up)) ~ Compartment, data = data.frame(sample_data(r_wk3_up)), permutations = perm, by = "margin")
                       Df SumOfSqs      R2      F Pr(>F)    
           Compartment  2   6.8091 0.69435 47.705  0.001 ***
           Residual    42   2.9974 0.30565                  
           Total       44   9.8065 1.00000                  
           ---
           Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
           ```

           * ---- Plot

           ```
           Permutation test for adonis under NA model
           Marginal effects of terms
           Plots: data.frame(sample_data(r_wk3_up))$Compartment, plot permutation: none
           Permutation: free
           Number of permutations: 999

           adonis2(formula = t(otu_table(r_wk3_up)) ~ Plot, data = data.frame(sample_data(r_wk3_up)), permutations = perm, by = "margin")
                    Df SumOfSqs      R2      F Pr(>F)    
           Plot     14   1.4567 0.14855 0.3738  0.001 ***
           Residual 30   8.3498 0.85145                  
           Total    44   9.8065 1.00000
           ```

           * ---- Treatment 

           ```
           Permutation test for adonis under NA model
           Marginal effects of terms
           Plots: data.frame(sample_data(r_wk3_up))$Compartment, plot permutation: none
           Permutation: free
           Number of permutations: 999

           adonis2(formula = t(otu_table(r_wk3_up)) ~ Treatment, data = data.frame(sample_data(r_wk3_up)), permutations = perm, by = "margin")
                     Df SumOfSqs      R2      F Pr(>F)  
           Treatment  2   0.2066 0.02106 0.4519  0.018 *
           Residual  42   9.6000 0.97894                
           Total     44   9.8065 1.00000                
           ---
           Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
           ```

           * ---- Read_depth vs compartment vs Plot

           ```
           Permutation test for adonis under reduced model
           Marginal effects of terms
           Permutation: free
           Number of permutations: 999

           adonis2(formula = t(otu_table(r_wk3_up)) ~ Compartment + Read_depth + Plot, data = data.frame(sample_data(r_wk3_up)), permutations = 999, by = "margin")
                       Df SumOfSqs      R2       F Pr(>F)    
           Compartment  2   5.2200 0.53230 47.4055  0.001 ***
           Read_depth   1   0.0541 0.00552  0.9833  0.346    
           Plot        14   1.4832 0.15124  1.9242  0.006 ** 
           Residual    27   1.4865 0.15159                   
           Total       44   9.8065 1.00000                   
           ---
           Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
           ```

           * ---- Read_depth vs compartment vs Treatment

           ```
           Permutation test for adonis under reduced model
           Marginal effects of terms
           Permutation: free
           Number of permutations: 999

           adonis2(formula = t(otu_table(r_wk3_up)) ~ Compartment + Read_depth + Treatment, data = data.frame(sample_data(r_wk3_up)), permutations = 999, by = "margin")
                       Df SumOfSqs      R2       F Pr(>F)    
           Compartment  2   5.5587 0.56684 39.2318  0.001 ***
           Read_depth   1   0.0279 0.00284  0.3937  0.942    
           Treatment    2   0.2068 0.02109  1.4594  0.152    
           Residual    39   2.7629 0.28175                   
           Total       44   9.8065 1.00000                   
           ---
           Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
           ```
**>>>>Week4>>>>**

          ```
          phyloseq-class experiment-level object
          otu_table()   OTU Table:         [ 11513 taxa and 45 samples ]
          sample_data() Sample Data:       [ 45 samples by 7 sample variables ]
          tax_table()   Taxonomy Table:    [ 11513 taxa by 6 taxonomic ranks ]
          ```

          * ---- Read_depth
          ```
          Permutation test for adonis under NA model
          Marginal effects of terms
          Permutation: free
          Number of permutations: 999

          adonis2(formula = t(otu_table(r_wk4_up)) ~ Read_depth, data = data.frame(sample_data(r_wk4_up)), permutations = 999, by = "margin")
                     Df SumOfSqs      R2      F Pr(>F)
          Read_depth  1   0.1900 0.01873 0.8206  0.472
          Residual   43   9.9576 0.98127              
          Total      44  10.1476 1.00000
          ```

          * ---- Compartment

          ```
          Permutation test for adonis under NA model
          Marginal effects of terms
          Plots: data.frame(sample_data(r_wk4_up))$Plot, plot permutation: none
          Permutation: free
          Number of permutations: 999

          adonis2(formula = t(otu_table(r_wk4_up)) ~ Compartment, data = data.frame(sample_data(r_wk4_up)), permutations = perm, by = "margin")
                      Df SumOfSqs      R2      F Pr(>F)    
          Compartment  2   6.6819 0.65847 40.489  0.001 ***
          Residual    42   3.4657 0.34153                  
          Total       44  10.1476 1.00000                  
          ---
          Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
          > 
          ```

          * ---- Plot

          ```
          Permutation test for adonis under NA model
          Marginal effects of terms
          Plots: data.frame(sample_data(r_wk4_up))$Compartment, plot permutation: none
          Permutation: free
          Number of permutations: 999

          adonis2(formula = t(otu_table(r_wk4_up)) ~ Plot, data = data.frame(sample_data(r_wk4_up)), permutations = perm, by = "margin")
                   Df SumOfSqs      R2      F Pr(>F)    
          Plot     14   1.7794 0.17535 0.4556  0.001 ***
          Residual 30   8.3682 0.82465                  
          Total    44  10.1476 1.00000
          ```

          * ---- Treatment

          ```
          Permutation test for adonis under NA model
          Marginal effects of terms
          Plots: data.frame(sample_data(r_wk4_up))$Compartment, plot permutation: none
          Permutation: free
          Number of permutations: 999

          adonis2(formula = t(otu_table(r_wk4_up)) ~ Treatment, data = data.frame(sample_data(r_wk4_up)), permutations = perm, by = "margin")
                    Df SumOfSqs      R2      F Pr(>F)   
          Treatment  2   0.2676 0.02637 0.5688  0.003 **
          Residual  42   9.8800 0.97363                 
          Total     44  10.1476 1.00000                 
          ---
          Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
          ```

          * ---- Read_depth vs Compartment vs Plot

          ```
          Permutation test for adonis under reduced model
          Marginal effects of terms
          Permutation: free
          Number of permutations: 999

          adonis2(formula = t(otu_table(r_wk4_up)) ~ Compartment + Read_depth + Plot, data = data.frame(sample_data(r_wk4_up)), permutations = 999, by = "margin")
                      Df SumOfSqs      R2       F Pr(>F)    
          Compartment  2   6.5726 0.64770 55.1227  0.001 ***
          Read_depth   1   0.0766 0.00755  1.2851  0.229    
          Plot        14   1.7326 0.17074  2.0758  0.001 ***
          Residual    27   1.6097 0.15863                   
          Total       44  10.1476 1.00000                   
          ---
          Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
          ```

          * ---- Read_depth vs Compartment vs Treatment

          ```
          Permutation test for adonis under reduced model
          Marginal effects of terms
          Permutation: free
          Number of permutations: 999

          adonis2(formula = t(otu_table(r_wk4_up)) ~ Compartment + Read_depth + Plot, data = data.frame(sample_data(r_wk4_up)), permutations = 999, by = "margin")
                      Df SumOfSqs      R2       F Pr(>F)    
          Compartment  2   6.5726 0.64770 55.1227  0.001 ***
          Read_depth   1   0.0766 0.00755  1.2851  0.229    
          Plot        14   1.7326 0.17074  2.0758  0.001 ***
          Residual    27   1.6097 0.15863                   
          Total       44  10.1476 1.00000                   
          ---
          Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
          ```
