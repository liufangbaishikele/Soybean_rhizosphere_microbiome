## Soybean root associated bacterial community dynamics in response to plant development

author | Date
-----| -----
Fang Liu | 08/02/2019

### Introduction

For this experiment, soybean bulk, rhizosphere and endosphere microbes were characterized based on 16S rRNA V3-V4 sequencing data. MOTHUR was used to process the sequencing data to generate community abundance dataframe. R packages and Python based softwares were used to analyze bacteria abundance data. [R Code couls be found here](https://github.com/liufangbaishikele/Soybean_rhizosphere_microbiome/blob/master/2018_fungicide/Seed_treatment/16S/2018_seed_fungicide_16S_13021.Rmd)

### Results summary


**BEfore rarefaction at 13021**

* Summary of the otu_table, tax_table and the sample_data
```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 144821 taxa and 164 samples ]
sample_data() Sample Data:       [ 164 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 144821 taxa by 6 taxonomic ranks ]
```
* In total, it has 13079622 reads after mothur pipeline before doing any rarefaction

**r_Seed_13021_up**

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

**r_Soil_up**

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

**More exploration with r_SB_CT_up**

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
**Subset phyloseq to BRE - including bulk, rhizosphere and endosphere samples**

* Summary of the otu_table, tax_table and metadata

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 19681 taxa and 134 samples ]
sample_data() Sample Data:       [ 134 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 19681 taxa by 6 taxonomic ranks ]
```
* PERMANOVA 

* Read_depth - insignificant

```
adonis2(formula = t(otu_table(r_BRE_up)) ~ Read_depth, data = data.frame(sample_data(r_BRE_up)), permutations = 999, by = "margin")
            Df SumOfSqs      R2      F Pr(>F)
Read_depth   1    0.373 0.01121 1.4958  0.144
Residual   132   32.886 0.98879              
Total      133   33.258 1.00000
```


* Compartment impact

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

* Time impact

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

* Plot impact vs Treatment impact

```

```





