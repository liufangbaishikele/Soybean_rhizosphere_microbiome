#                                     Community analysis summary


----
Fang Liu
08/02/2019


**Overal description**

This is the results of Combined run of ITS sequencing targeted ITS2 region. For the first run, we got very limited reads left for endosphere and some rhizosphere samples. So, we design ITS blocker to block off soybean ITS, but it did not work as good as we expected. However, for most of the samples we got enough reads for analysis.


### All samples without rarefaction, but singletons were removed

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8360 taxa and 164 samples ]
sample_data() Sample Data:       [ 164 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 8360 taxa by 7 taxonomic ranks ]
```

* **Read_depth**

```
adonis2(formula = t(otu_table(r_rms_seed_phyloseq)) ~ Read_depth, data = data.frame(sample_data(r_rms_seed_phyloseq)), permutations = 999, by = "margin")
            Df SumOfSqs      R2      F Pr(>F)    
Read_depth   1    3.792 0.08548 15.142  0.001 ***
Residual   162   40.566 0.91452                  
Total      163   44.358 1.00000                 
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


### After remove all the seed samples: non_seed samples

* **Read_Depth**

```
adonis2(formula = t(otu_table(r_rms_non_Seed_up)) ~ Read_depth, data = data.frame(sample_data(r_rms_non_Seed_up)), permutations = 999, by = "margin")
            Df SumOfSqs     R2      F Pr(>F)    
Read_depth   1    2.644 0.0749 11.902  0.001 ***
Residual   147   32.661 0.9251                  
Total      148   35.306 1.0000
```

### Subset to rms_Soil_up 
In order to check the clustering of plots before sowing soybean seeds

**Treatment impact** -- insignificant

```
adonis2(formula = t(otu_table(r_rms_Soil_up)) ~ Treatment, data = data.frame(sample_data(r_rms_Soil_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)
Treatment  2  0.19993 0.14733 1.0367  0.359
Residual  12  1.15709 0.85267              
Total     14  1.35702 1.00000
```
### Subset to BRE - Bulk, rhizosphere and endosphere 

* **Read_depth**

```
adonis2(formula = t(otu_table(r_rms_BRE_up)) ~ Read_depth, data = data.frame(sample_data(r_rms_BRE_up)), permutations = 999, by = "margin")
            Df SumOfSqs     R2      F Pr(>F)    
Read_depth   1    2.554 0.0789 11.307  0.001 ***
Residual   132   29.818 0.9211                  
Total      133   32.372 1.0000
```

* **Plot impacts**

```
adonis2(formula = t(otu_table(r_rms_BRE_up)) ~ Plot, data = data.frame(sample_data(r_rms_BRE_up)), permutations = perm, by = "margin")
Model: adonis0(formula = lhs ~ Plot, data = data, method = method)
          Df SumOfSqs      R2      F Pr(>F)    
Plot      14    5.192 0.16037 1.6236  0.001 ***
Residual 119   27.180 0.83963                  
Total    133   32.372 1.00000 
```
* **Compartment**

```
adonis2(formula = t(otu_table(r_rms_BRE_up)) ~ Compartment, data = data.frame(sample_data(r_rms_BRE_up)), permutations = perm, by = "margin")
Model: adonis0(formula = lhs ~ Compartment, data = data, method = method)
             Df SumOfSqs      R2      F Pr(>F)    
Compartment   2   10.182 0.31452 30.054  0.001 ***
Residual    131   22.190 0.68548                  
Total       133   32.372 1.00000
```

* **Time impact**

```
Marginal effects of terms
Plots: paste(data.frame(sample_data(r_rms_BRE_up))$Compartment, data.frame(sample_data(r_rms_BRE_up))$Plot, , plot permutation: none
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(r_rms_BRE_up)) ~ Time, data = data.frame(sample_data(r_rms_BRE_up)), permutations = perm, by = "margin")
Model: adonis0(formula = lhs ~ Time, data = data, method = method)
          Df SumOfSqs      R2      F Pr(>F)    
Time       2    1.885 0.05824 4.0504  0.001 ***
Residual 131   30.487 0.94176                  
Total    133   32.372 1.00000 
```

* **Read_depth vs Compartment vs Time vs Plot**

```
adonis2(formula = t(otu_table(r_rms_BRE_up)) ~ Read_depth + Compartment + Time + Plot, data = data.frame(sample_data(r_rms_BRE_up)), permutations = 999, by = "margin")
             Df SumOfSqs      R2       F Pr(>F)    
Read_depth    1    0.274 0.00846  2.1010  0.033 *  
Compartment   2    7.789 0.24060 29.8663  0.001 ***
Time          2    1.887 0.05830  7.2367  0.001 ***
Plot         14    5.179 0.15997  2.8368  0.001 ***
Residual    114   14.865 0.45919                   
Total       133   32.372 1.00000 
```
* **Read_depth vs Compartment vs Time vs Treatment**

```
adonis2(formula = t(otu_table(r_rms_BRE_up)) ~ Read_depth + Compartment + Time + Treatment, data = data.frame(sample_data(r_rms_BRE_up)), permutations = 999, by = "margin")
             Df SumOfSqs      R2       F Pr(>F)    
Read_depth    1    0.273 0.00842  1.8023  0.059 .  
Compartment   2    7.853 0.24259 25.9574  0.001 ***
Time          2    1.893 0.05849  6.2580  0.001 ***
Treatment     2    0.983 0.03038  3.2506  0.001 ***
Residual    126   19.060 0.58878                   
Total       133   32.372 1.00000
```

### subset to Bulk -- r_rms_Bulk vs r_rf_rms_Bulk vs r_norm_Bulk. Here, rms means singletons were removed; rf means samples were rarefied to minimum read depth after remove singletons; 

**r_rms_Bulk**

* --- Read_depth ----
```
adonis2(formula = t(otu_table(r_rms_Bulk_up)) ~ Read_depth, data = data.frame(sample_data(r_rms_Bulk_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)  
Read_depth  1   0.2703 0.04328 1.9454  0.014 *
Residual   43   5.9736 0.95672                
Total      44   6.2439 1.00000
```

* ---- Time ---

```
Permutation test for adonis under NA model
Marginal effects of terms
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(r_rms_Bulk_up)) ~ Time, data = data.frame(sample_data(r_rms_Bulk_up)), permutations = 999, by = "margin")
         Df SumOfSqs      R2      F Pr(>F)    
Time      2   1.3833 0.22155 5.9768  0.001 ***
Residual 42   4.8605 0.77845                  
Total    44   6.2439 1.00000
```
* ---- Plot -----

```
Permutation test for adonis under NA model
Marginal effects of terms
Plots: data.frame(sample_data(r_rms_Bulk_up))$Time, plot permutation: none
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(r_rms_Bulk_up)) ~ Plot, data = data.frame(sample_data(r_rms_Bulk_up)), permutations = perm, by = "margin")
         Df SumOfSqs      R2      F Pr(>F)    
Plot     14   2.4409 0.39092 1.3753  0.001 ***
Residual 30   3.8030 0.60908                  
Total    44   6.2439 1.00000
```
* ---- Treatment ----

```
Permutation test for adonis under NA model
Marginal effects of terms
Plots: data.frame(sample_data(r_rms_Bulk_up))$Time, plot permutation: none
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(r_rms_Bulk_up)) ~ Treatment, data = data.frame(sample_data(r_rms_Bulk_up)), permutations = perm, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)   
Treatment  2   0.3617 0.05793 1.2913  0.003 **
Residual  42   5.8822 0.94207                 
Total     44   6.2439 1.00000                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
* ----Read_depth vs Time vs Plot

```
adonis2(formula = t(otu_table(r_rms_Bulk_up)) ~ Read_depth + Plot + Time, data = data.frame(sample_data(r_rms_Bulk_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)    
Read_depth  1   0.0995 0.01593 1.1573  0.222    
Plot       14   2.4008 0.38450 1.9955  0.001 ***
Time        2   1.2131 0.19428 7.0581  0.001 ***
Residual   27   2.3202 0.37160                  
Total      44   6.2439 1.00000
```

* ---- Read_depth vs Time vs Treatment

```
adonis2(formula = t(otu_table(r_rms_Bulk_up)) ~ Read_depth + Treatment + Time, data = data.frame(sample_data(r_rms_Bulk_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)    
Read_depth  1   0.1505 0.02410 1.3496  0.112    
Treatment   2   0.3726 0.05968 1.6711  0.012 *  
Time        2   1.2374 0.19817 5.5489  0.001 ***
Residual   39   4.3484 0.69642                  
Total      44   6.2439 1.00000
```
--------------------------------------------------------------------------------                        


### Subset to Rhizosphere -- before rarefaction

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 5020 taxa and 44 samples ]
sample_data() Sample Data:       [ 44 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 5020 taxa by 7 taxonomic ranks ]
````
* -----Read_depth---

```
adonis2(formula = t(otu_table(r_rms_Rhi_up)) ~ Read_depth, data = data.frame(sample_data(r_rms_Rhi_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2    F Pr(>F)  
Read_depth  1   0.3343 0.04459 1.96  0.016 *
Residual   42   7.1636 0.95541              
Total      43   7.4979 1.00000
```

* --- Plot----

```
Permutation test for adonis under NA model
Marginal effects of terms
Plots: data.frame(sample_data(r_rms_Rhi_up))$Time, plot permutation: none
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(r_rms_Rhi_up)) ~ Plot, data = data.frame(sample_data(r_rms_Rhi_up)), permutations = perm, by = "margin")
         Df SumOfSqs      R2      F Pr(>F)    
Plot     14   3.1112 0.41495 1.4691  0.001 ***
Residual 29   4.3867 0.58505                  
Total    43   7.4979 1.00000
```

* ---- Time ----

```
aPermutation test for adonis under NA model
Marginal effects of terms
Plots: data.frame(sample_data(r_rms_Rhi_up))$Plot, plot permutation: none
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(r_rms_Rhi_up)) ~ Time, data = data.frame(sample_data(r_rms_Rhi_up)), permutations = perm, by = "margin")
         Df SumOfSqs      R2      F Pr(>F)    
Time      2   1.4557 0.19414 4.9387  0.001 ***
Residual 41   6.0423 0.80586                  
Total    43   7.4979 1.00000
```

* --- Treatment

```
Permutation test for adonis under NA model
Marginal effects of terms
Plots: data.frame(sample_data(r_rms_Rhi_up))$Time, plot permutation: none
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(r_rms_Rhi_up)) ~ Treatment, data = data.frame(sample_data(r_rms_Rhi_up)), permutations = perm, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)    
Treatment  2   0.7743 0.10326 2.3606  0.001 ***
Residual  41   6.7237 0.89674                  
Total     43   7.4979 1.00000
```

* --- Read_depth vs Time vs Plot 

```
           Df SumOfSqs      R2      F Pr(>F)    
Read_depth  1   0.1177 0.01570 1.0813  0.355    
Time        2   1.2127 0.16174 5.5714  0.001 ***
Plot       14   3.0498 0.40675 2.0015  0.001 ***
Residual   26   2.8298 0.37741                  
Total      43   7.4979 1.00000
```
* ---- Read_depth vs Time vs Treatment

```
term          df SumOfSqs     R2 statistic p.value
  <chr>      <dbl>    <dbl>  <dbl>     <dbl>   <dbl>
1 Read_depth     1    0.149 0.0199      1.11   0.306
2 Time           2    1.28  0.171       4.76   0.001
3 Treatment      2    0.768 0.102       2.86   0.001
4 Residual      38    5.11  0.682      NA     NA    
5 Total         43    7.50  1          NA     NA 
```


### Subset to Endosphere - without rarefaction due to the very limited sequencing depth

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 2391 taxa and 45 samples ]
sample_data() Sample Data:       [ 45 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 2391 taxa by 7 taxonomic ranks ]
```

* --- Read_depth ---

```
adonis2(formula = t(otu_table(r_rms_Endo_up)) ~ Read_depth, data = data.frame(sample_data(r_rms_Endo_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)  
Read_depth  1   0.3568 0.04236 1.9021  0.012 *
Residual   43   8.0660 0.95764                
Total      44   8.4228 1.00000
```

* ---- Plot ----

```
ermutation test for adonis under NA model
Marginal effects of terms
Plots: data.frame(sample_data(r_rms_Endo_up))$Time, plot permutation: none
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(r_rms_Endo_up)) ~ Plot, data = data.frame(sample_data(r_rms_Endo_up)), permutations = perm, by = "margin")
         Df SumOfSqs      R2      F Pr(>F)    
Plot     14   3.3155 0.39364 1.3911  0.001 ***
Residual 30   5.1073 0.60636                  
Total    44   8.4228 1.00000 
```
* ---- Treatment ---

```
Permutation test for adonis under NA model
Marginal effects of terms
Plots: data.frame(sample_data(r_rms_Endo_up))$Time, plot permutation: none
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(r_rms_Endo_up)) ~ Treatment, data = data.frame(sample_data(r_rms_Endo_up)), permutations = perm, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)   
Treatment  2   0.5881 0.06983 1.5764  0.002 **
Residual  42   7.8347 0.93017                 
Total     44   8.4228 1.00000
```
* --- Time ----

```
Permutation test for adonis under NA model
Marginal effects of terms
Plots: data.frame(sample_data(r_rms_Endo_up))$Plot, plot permutation: none
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(r_rms_Endo_up)) ~ Time, data = data.frame(sample_data(r_rms_Endo_up)), permutations = perm, by = "margin")
         Df SumOfSqs      R2      F Pr(>F)    
Time      2   1.2603 0.14963 3.6951  0.001 ***
Residual 42   7.1625 0.85037                  
Total    44   8.4228 1.00000                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```


* ---- Read_depth vs Time vs Plot

```
term          df SumOfSqs     R2 statistic p.value
  <chr>      <dbl>    <dbl>  <dbl>     <dbl>   <dbl>
1 Plot          14    3.23  0.383       1.71   0.001
2 Time           2    1.17  0.138       4.32   0.001
3 Read_depth     1    0.204 0.0242      1.51   0.065
4 Residual      27    3.64  0.433      NA     NA    
5 Total         44    8.42  1          NA     NA 
```

* --- Read_depth vs Time vs Treatment

```
term          df SumOfSqs     R2 statistic p.value
  <chr>      <dbl>    <dbl>  <dbl>     <dbl>   <dbl>
1 Treatment      2    0.539 0.0639      1.66   0.011
2 Time           2    1.18  0.140       3.64   0.001
3 Read_depth     1    0.242 0.0288      1.49   0.068
4 Residual      39    6.33  0.752      NA     NA    
5 Total         44    8.42  1          NA     NA
```


### Subset to SEED samples

* ---- Treatment impact ----

```
Permutation test for adonis under NA model
Marginal effects of terms
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(r_rms_SEED)) ~ Treatment, data = data.frame(sample_data(r_rms_SEED)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)  
Treatment  2  0.45624 0.20242 1.5228  0.079 .
Residual  12  1.79765 0.79758                
Total     14  2.25389 1.00000
```

### Subset based on development stage - week1, week3 and week4




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


