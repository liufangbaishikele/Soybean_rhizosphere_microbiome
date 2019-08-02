## Bacterial community analysis were conducted in R

### All samples - r_Seed_22017_up

1. **Read_depth** --- significant but explain only **2%** of variations

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

2. **Compartment** --- Significant and can explain **57.46%** variations

```
adonis2(formula = t(otu_table(r_Seed_22017_up)) ~ Compartment, data = data.frame(sample_data(r_Seed_22017_up)), permutations = 999, by = "margin")
             Df SumOfSqs     R2      F Pr(>F)    
Compartment   4   24.621 0.5746 52.004  0.001 ***
Residual    154   18.228 0.4254                  
Total       158   42.849 1.0000  
```

3. **Read_depth vs Compartment** -- Read_depth no more significant, but compartment still significant, explaining **55.64%** variations

```
Permutation test for adonis under reduced model
Marginal effects of terms
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(r_Seed_22017_up)) ~ Read_depth + Compartment, data = data.frame(sample_data(r_Seed_22017_up)), permutations = 999, by = "margin")
             Df SumOfSqs      R2       F Pr(>F)    
Read_depth    1    0.117 0.00273  0.9883  0.402    
Compartment   4   23.839 0.55635 50.3479  0.001 ***
Residual    153   18.111 0.42267                   
Total       158   42.849 1.00000  
```

### Remove seed samples -- non_Seed_22017_up

1. **Plot variations** --- insignificant p=0.643

- Subset phyloseq object to only include non-seed samples: PERMANOVA results indicate an non-significant plots impacts on overal bacterial community composition.

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
2. **Read_depth** -- Significant p=0.033, but only explain **1.69%** variations

```
adonis2(formula = t(otu_table(r_non_Seed_22017_up)) ~ Read_depth, data = data.frame(sample_data(r_non_Seed_22017_up)), permutations = 999, by = "margin")
            Df SumOfSqs      R2      F Pr(>F)  
Read_depth   1    0.586 0.01691 2.4934  0.033 *
Residual   145   34.062 0.98309                
Total      146   34.648 1.00000    
```

3. **Compartment** -- significant and can explain **47.83%** of variations

```
adonis2(formula = t(otu_table(r_non_Seed_22017_up)) ~ Compartment, data = data.frame(sample_data(r_non_Seed_22017_up)), permutations = 999, by = "margin")
             Df SumOfSqs      R2      F Pr(>F)    
Compartment   3   16.573 0.47832 43.705  0.001 ***
Residual    143   18.075 0.52168                  
Total       146   34.648 1.00000 
```
4. **Time** -- significant and could explain 13.43% of variations

```
adonis2(formula = t(otu_table(r_non_Seed_22017_up)) ~ Time, data = data.frame(sample_data(r_non_Seed_22017_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)    
Time       3    4.652 0.13426 7.3922  0.001 ***
Residual 143   29.996 0.86574                  
Total    146   34.648 1.00000                  
```

5. **Treatment* -- 

```
adonis2(formula = t(otu_table(r_non_Seed_22017_up)) ~ Treatment, data = data.frame(sample_data(r_non_Seed_22017_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)
Treatment   2    0.464 0.01339 0.9769  0.425
Residual  144   34.184 0.98661              
Total     146   34.648 1.00000    
```
6. **Overal Read_depth vs Compartment vs Time vs Treatment** Read_depth no more significant. Compartment, Time and Treatment all significant, but the impact of treatment is very small.

```
adonis2(formula = t(otu_table(r_non_Seed_22017_up)) ~ Read_depth + Treatment + Time + Compartment, data = data.frame(sample_data(r_non_Seed_22017_up)), permutations = 999, by = "margin")
             Df SumOfSqs      R2       F Pr(>F)    
Read_depth    1    0.112 0.00324  1.0787  0.309    
Treatment     2    0.429 0.01238  2.0592  0.029 *  
Time          2    3.151 0.09095 15.1276  0.001 ***
Compartment   2   14.680 0.42369 70.4694  0.001 ***
Residual    138   14.374 0.41485                   
Total       146   34.648 1.00000
```


### Compartment impacts changes along soybean development

1. **1 week old soybeans**  --- significant and could explain **61.28%** of the total variations

```
#-----Read_depth-----insignificant

adonis2(formula = t(otu_table(wk1_up)) ~ Read_depth, data = data.frame(sample_data(wk1_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)
Read_depth  1   0.0973 0.01003 0.4256  0.897
Residual   42   9.5993 0.98997              
Total      43   9.6965 1.00000

#-----Plot----- insignificant
set.seed(1013) #- insignificant
>adonis2(t(otu_table(wk1_up ))~Plot,data=data.frame(sample_data(wk1_up )),permutations=999,by="margin")
Permutation test for adonis under NA model
Marginal effects of terms
Permutation: free
Number of permutations: 999

adonis2(formula = t(otu_table(wk1_up)) ~ Plot, data = data.frame(sample_data(wk1_up)), permutations = 999, by = "margin")
         Df SumOfSqs      R2      F Pr(>F)
Plot     14   1.8483 0.19062 0.4878      1
Residual 29   7.8482 0.80938              
Total    43   9.6965 1.00000

# ----Compartment ----- significant and could explain 61.28% of the total varations

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

Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# ---- Fungicide treatment -----
adonis2(formula = t(otu_table(wk1_up)) ~ Treatment, data = data.frame(sample_data(wk1_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)
Treatment  2   0.2239 0.02309 0.4845  0.922
Residual  41   9.4726 0.97691              
Total     43   9.6965 1.00000

# ---- Compartment vs Read_depth--- significant and explain 61.34% variance
adonis2(formula = t(otu_table(wk1_up)) ~ Compartment + Read_depth, data = data.frame(sample_data(wk1_up)), permutations = 999, by = "margin")
            Df SumOfSqs      R2       F Pr(>F)    
Compartment  2   5.9479 0.61340 32.5788  0.001 ***
Read_depth   1   0.1029 0.01061  1.1274  0.337    
Residual    40   3.6514 0.37657                   
Total       43   9.6965 1.00000

```

2. **3 weeks old soybeans**

```
#-------Plot -------

adonis2(formula = t(otu_table(wk3_up)) ~ Plot, data = data.frame(sample_data(wk3_up)), permutations = 999, by = "margin")
         Df SumOfSqs      R2      F Pr(>F)
Plot     14   1.3592 0.14218 0.3552      1
Residual 30   8.2008 0.85782              
Total    44   9.5600 1.00000 

#-----Read_depth ----- significant, but this is due to the relatively high sequencing depth for endo but small for rhizosphere

adonis2(formula = t(otu_table(wk3_up)) ~ Read_depth, data = data.frame(sample_data(wk3_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)    
Read_depth  1   1.2708 0.13293 6.5921  0.001 ***
Residual   43   8.2892 0.86707                  
Total      44   9.5600 1.00000   

#-----Compartment-----
adonis2(formula = t(otu_table(wk3_up)) ~ Compartment, data = data.frame(sample_data(wk3_up)), permutations = 999, by = "margin")
            Df SumOfSqs      R2      F Pr(>F)    
Compartment  2   6.8324 0.71468 52.602  0.001 ***
Residual    42   2.7276 0.28532                  
Total       44   9.5600 1.00000  

# --- Compartment vs Read_depth --- significant compartment, 58.45% variance

adonis2(formula = t(otu_table(wk3_up)) ~ Read_depth + Compartment, data = data.frame(sample_data(wk3_up)), permutations = 999, by = "margin")
            Df SumOfSqs      R2       F Pr(>F)    
Read_depth   1   0.0260 0.00272  0.3941  0.928    
Compartment  2   5.5876 0.58447 42.3981  0.001 ***
Residual    41   2.7016 0.28260                   
Total       44   9.5600 1.00000

```

3. **4 weeks old soybeans**

```
# --- Plot ---- insignificant
adonis2(formula = t(otu_table(wk4_up)) ~ Plot, data = data.frame(sample_data(wk4_up)), permutations = 999, by = "margin")
         Df SumOfSqs      R2      F Pr(>F)
Plot     14   1.8590 0.19732 0.4916  0.998
Residual 28   7.5623 0.80268              
Total    42   9.4213 1.00000

#---- Read_depth ---- insignificant

adonis2(formula = t(otu_table(wk4_up)) ~ Read_depth, data = data.frame(sample_data(wk4_up)), permutations = 999, by = "margin")
           Df SumOfSqs     R2      F Pr(>F)
Read_depth  1   0.1733 0.0184 0.7684  0.505
Residual   41   9.2480 0.9816              
Total      42   9.4213 1.0000

# ---- Treatment ----- insignificant
adonis2(formula = t(otu_table(wk4_up)) ~ Treatment, data = data.frame(sample_data(wk4_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)
Treatment  2   0.2972 0.03155 0.6515  0.689
Residual  40   9.1241 0.96845              
Total     42   9.4213 1.00000 

# ---- Compartment ---- significant and 67.22% could be explained by compartment difference

adonis2(formula = t(otu_table(wk4_up)) ~ Compartment, data = data.frame(sample_data(wk4_up)), permutations = 999, by = "margin")
            Df SumOfSqs      R2     F Pr(>F)    
Compartment  2   6.3334 0.67224 41.02  0.001 ***
Residual    40   3.0879 0.32776                 
Total       42   9.4213 1.00000                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# ---- Compartment vs Read_depth --- significant and could explain 66.7% variance

adonis2(formula = t(otu_table(wk4_up)) ~ Compartment + Read_depth, data = data.frame(sample_data(wk4_up)), permutations = 999, by = "margin")
            Df SumOfSqs      R2       F Pr(>F)    
Compartment  2   6.2860 0.66722 41.3843  0.001 ***
Read_depth   1   0.1260 0.01337  1.6587  0.132    
Residual    39   2.9619 0.31439                   
Total       42   9.4213 1.00000

```
**Considering the significant impacts of compartment, soybean development (Time) and seed fungicide treatment impacts will be analyzed after subset phyloseq to corresponding compartment**

## Treatment and soybean development impacts difference among compartments

1. **Bulk samples**

```
# --- Plot impacts --- significant, explaining 51.08% variance

adonis2(formula = t(otu_table(r_Bulk_up)) ~ Plot, data = data.frame(sample_data(r_Bulk_up)), permutations = 999, by = "margin")
         Df SumOfSqs      R2      F Pr(>F)    
Plot     14   1.9072 0.51077 2.1626  0.001 ***
Residual 29   1.8268 0.48923                  
Total    43   3.7340 1.00000                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#---- Read_depth --- non significant

adonis2(formula = t(otu_table(r_Bulk_up)) ~ Read_depth, data = data.frame(sample_data(r_Bulk_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)  
Read_depth  1   0.1249 0.03346 1.4541  0.071 .
Residual   42   3.6090 0.96654                
Total      43   3.7340 1.00000

# ---- Time impact ----- significant and could explain 12.47% variance

adonis2(formula = t(otu_table(r_Bulk_up)) ~ Time, data = data.frame(sample_data(r_Bulk_up)), permutations = 999, by = "margin")
         Df SumOfSqs     R2      F Pr(>F)    
Time      2   0.4656 0.1247 2.9205  0.001 ***
Residual 41   3.2683 0.8753                  
Total    43   3.7340 1.0000 

# --- Treatment ----- significant but does this because the plot impact or itself?
adonis2(formula = t(otu_table(r_Bulk_up)) ~ Treatment, data = data.frame(sample_data(r_Bulk_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)  
Treatment  2   0.2698 0.07226 1.5967  0.014 *
Residual  41   3.4642 0.92774                
Total     43   3.7340 1.00000

# --- Treament vs Time ---

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

## --- Treament vs Time --- Interaction effects

adonis2(formula = t(otu_table(r_Bulk_22017_up)) ~ Treatment * Time, data = data.frame(sample_data(r_Bulk_22017_up)), permutations = 999, by = "margin")
               Df SumOfSqs      R2      F Pr(>F)
Treatment:Time  4   0.2040 0.05464 0.6387      1
Residual       35   2.7952 0.74859              
Total          43   3.7340 1.00000    

# ----- Read_depth vs Time vs Treatment ---

adonis2(formula = t(otu_table(r_Bulk_up)) ~ Read_depth + Treatment + Time, data = data.frame(sample_data(r_Bulk_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)    
Read_depth  1   0.0728 0.01950 0.9457  0.499    
Treatment   2   0.2151 0.05762 1.3968  0.067 .  
Time        2   0.4540 0.12159 2.9478  0.001 ***
Residual   38   2.9264 0.78373                  
Total      43   3.7340 1.00000 

# ----- Read_depth vs Time vs Plot ---

adonis2(formula = t(otu_table(r_Bulk_up)) ~ Read_depth + Plot + Time, data = data.frame(sample_data(r_Bulk_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)    
Read_depth  1   0.0470 0.01259 0.9249  0.527    
Plot       14   1.8204 0.48753 2.5590  0.001 ***
Time        2   0.4219 0.11298 4.1511  0.001 ***
Residual   26   1.3211 0.35381                  
Total      43   3.7340 1.00000

# -------- Treatment impact of Bulk_wk1, Bulk_wk3 and Bulk_wk4 -----

>>>>  Bulk_1wk

adonis2(formula = t(otu_table(r_Bulk_1wk_up)) ~ Treatment, data = data.frame(sample_data(r_Bulk_1wk_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2     F Pr(>F)
Treatment  2  0.13076 0.12397 0.849  0.732
Residual  12  0.92405 0.87603             
Total     14  1.05482 1.00000   


>>> Bulk_3wk 
adonis2(formula = t(otu_table(r_Bulk_3wk_up)) ~ Treatment, data = data.frame(sample_data(r_Bulk_3wk_up)), permutations = 999, by = "margin")
          Df SumOfSqs     R2      F Pr(>F)
Treatment  2  0.15486 0.1394 0.9718  0.507
Residual  12  0.95610 0.8606              
Total     14  1.11096 1.0000 

>>> Bulk_4wk 

adonis2(formula = t(otu_table(r_Bulk_4wk_up)) ~ Treatment, data = data.frame(sample_data(r_Bulk_4wk_up)), permutations = 999, by = "margin")
          Df SumOfSqs     R2      F Pr(>F)
Treatment  2  0.17109 0.1727 1.1481  0.205
Residual  11  0.81960 0.8273              
Total     13  0.99068 1.0000 
```

2. **Rhizosphere**

```
# ------ Plot ---- insignificant
adonis2(formula = t(otu_table(Rhi_up)) ~ Plot, data = data.frame(sample_data(Rhi_up)), permutations = 999, by = "margin")
         Df SumOfSqs      R2      F Pr(>F)
Plot     14   2.0502 0.30664 0.9161  0.664
Residual 29   4.6358 0.69336              
Total    43   6.6861 1.00000

# ---- Read_depth ----- insignificant
adonis2(formula = t(otu_table(Rhi_up)) ~ Read_depth, data = data.frame(sample_data(Rhi_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)
Read_depth  1   0.2380 0.03559 1.5501  0.126
Residual   42   6.4481 0.96441              
Total      43   6.6861 1.00000

# ----- Time ------ significant, R2=42.69%

adonis2(formula = t(otu_table(Rhi_up)) ~ Time, data = data.frame(sample_data(Rhi_up)), permutations = 999, by = "margin")
         Df SumOfSqs      R2      F Pr(>F)    
Time      2   2.8546 0.42694 15.273  0.001 ***
Residual 41   3.8315 0.57306                  
Total    43   6.6861 1.00000

# ------ Treatment ------ insignificant
adonis2(formula = t(otu_table(Rhi_up)) ~ Treatment, data = data.frame(sample_data(Rhi_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)
Treatment  2   0.3133 0.04686 1.0078  0.371
Residual  41   6.3728 0.95314              
Total     43   6.6861 1.00000

# ---- Treatment vs Time ---

adonis2(formula = t(otu_table(Rhi_22017_up)) ~ Treatment + Time, data = data.frame(sample_data(Rhi_22017_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2       F Pr(>F)    
Treatment  2   0.2985 0.04465  1.6478  0.102    
Time       2   2.8398 0.42473 15.6741  0.001 ***
Residual  39   3.5330 0.52841                   
Total     43   6.6861 1.00000                   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# --- Interactions between Treatment*Time ---

adonis2(formula = t(otu_table(Rhi_up)) ~ Treatment * Time, data = data.frame(sample_data(Rhi_up)), permutations = 999, by = "margin")
               Df SumOfSqs      R2      F Pr(>F)
Treatment:Time  4   0.2617 0.03914 0.6999  0.836
Residual       35   3.2713 0.48927              
Total          43   6.6861 1.00000

# --- Read_depth vs Time vs Treatment ---
adonis2(formula = t(otu_table(Rhi_up)) ~ Read_depth + Treatment + Time, data = data.frame(sample_data(Rhi_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2       F Pr(>F)    
Read_depth  1   0.1235 0.01848  1.3769  0.177    
Treatment   2   0.2451 0.03666  1.3659  0.184    
Time        2   2.6960 0.40323 15.0242  0.001 ***
Residual   38   3.4094 0.50993                   
Total      43   6.6861 1.00000

# Treatment impacts of Rhi_wk1, Rhi_wk3 and Rhi_wk4

>>> Rhi_1wk
adonis2(formula = t(otu_table(r_Rhi_1wk_up)) ~ Treatment, data = data.frame(sample_data(r_Rhi_1wk_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)
Treatment  2  0.14647 0.12323 0.7731  0.796
Residual  11  1.04211 0.87677              
Total     13  1.18858 1.00000

>>> Rhi_3wk --- insignificant
adonis2(formula = t(otu_table(r_Rhi_3wk_up)) ~ Treatment, data = data.frame(sample_data(r_Rhi_3wk_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)
Treatment  2  0.16711 0.15673 1.1152  0.267
Residual  12  0.89910 0.84327              
Total     14  1.06621 1.00000

>>> Rhi_4wk -- insignificant
adonis2(formula = t(otu_table(r_Rhi_4wk_up)) ~ Treatment, data = data.frame(sample_data(r_Rhi_4wk_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2     F Pr(>F)
Treatment  2  0.23425 0.15671 1.115  0.234
Residual  12  1.26053 0.84329             
Total     14  1.49478 1.00000

```

3. **Endosphere**

```
# --- Plot ---
adonis2(formula = t(otu_table(Endo_up)) ~ Plot, data = data.frame(sample_data(Endo_up)), permutations = 999, by = "margin")
         Df SumOfSqs      R2     F Pr(>F)
Plot     14   1.0953 0.17258 0.432      1
Residual 29   5.2512 0.82742             
Total    43   6.3465 1.00000

# --- Read_depth ---
adonis2(formula = t(otu_table(Endo_up)) ~ Read_depth, data = data.frame(sample_data(Endo_up)), permutations = 999, by = "margin")
           Df SumOfSqs      R2      F Pr(>F)
Read_depth  1   0.1416 0.02232 0.9587  0.376
Residual   42   6.2049 0.97768              
Total      43   6.3465 1.00000

# --- Time ----
adonis2(formula = t(otu_table(Endo_up)) ~ Time, data = data.frame(sample_data(Endo_up)), permutations = 999, by = "margin")
         Df SumOfSqs      R2      F Pr(>F)    
Time      2   3.8475 0.60624 31.562  0.001 ***
Residual 41   2.4990 0.39376                  
Total    43   6.3465 1.0000

# --- Treatment ----
adonis2(formula = t(otu_table(Endo_up)) ~ Treatment, data = data.frame(sample_data(Endo_up)), permutations = 999, by = "margin")
          Df SumOfSqs     R2      F Pr(>F)
Treatment  2   0.1656 0.0261 0.5494  0.762
Residual  41   6.1809 0.9739              
Total     43   6.3465 1.0000

# ---- Treatment + Time ---
adonis2(formula = t(otu_table(Endo_up)) ~ Treatment + Time, data = data.frame(sample_data(Endo_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2       F Pr(>F)    
Treatment  2   0.1546 0.02435  1.2856  0.233    
Time       2   3.8364 0.60449 31.9095  0.001 ***
Residual  39   2.3445 0.36941                   
Total     43   6.3465 1.00000

# --- interactions between Treatment*Time----
adonis2(formula = t(otu_table(Endo_up)) ~ Treatment * Time, data = data.frame(sample_data(Endo_up)), permutations = 999, by = "margin")
               Df SumOfSqs      R2      F Pr(>F)
Treatment:Time  4   0.1736 0.02736 0.6999   0.72
Residual       35   2.1708 0.34205              
Total          43   6.3465 1.00000

# Treatment impact along Endo_1wk, Endo_3wk and Endo_4wk 

>>> Endo_1wk 
adonis2(formula = t(otu_table(r_Endo_1wk_up)) ~ Treatment, data = data.frame(sample_data(r_Endo_1wk_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)
Treatment  2  0.19072 0.13268 0.9179  0.607
Residual  12  1.24667 0.86732              
Total     14  1.43739 1.00000

>>> Endo_3wk 
adonis2(formula = t(otu_table(r_Endo_3wk_up)) ~ Treatment, data = data.frame(sample_data(r_Endo_3wk_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)
Treatment  2  0.06663 0.13765 0.9578  0.509
Residual  12  0.41741 0.86235              
Total     14  0.48404 1.00000

>>> Endo_4wk 
adonis2(formula = t(otu_table(r_Endo_4wk_up)) ~ Treatment, data = data.frame(sample_data(r_Endo_4wk_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)
Treatment  2  0.06446 0.11923 0.7445   0.82
Residual  11  0.47620 0.88077              
Total     13  0.54066 1.00000
```

4. **Soil** -- Soil samples were collected before sowing soybean


```
# ---- Treatment ---- insignificant

adonis2(formula = t(otu_table(r_Soil_up)) ~ Treatment, data = data.frame(sample_data(r_Soil_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)
Treatment  2  0.16189 0.14502 1.0177  0.405
Residual  12  0.95439 0.85498              
Total     14  1.11628 1.00000 
```

The insignificant fungicide treatment impact could because fungicide does not have impact on fungal community. Alternatively, this could be due to the driving impact of Plot difference as shown across all bulk soil samples. From this point, we realized that we need replicates within each plot for future fungicide treatment experiment.



5. **SEED**

For the SEED samples, due to the limits of PNA in blocking the amplification of SEED mitochondria, we have lots of bacteria_unclassified in the SEED samples. So, before analysis of the fungicide treatment, we removed those bacteria_unclassified from the OTU table.

* PERMAOVA analysis 
```
# ------ Treatment ----

adonis2(formula = t(otu_table(r_SEED4_up)) ~ Treatment, data = data.frame(sample_data(r_SEED4_up)), permutations = 999, by = "margin")
          Df SumOfSqs      R2      F Pr(>F)    
Treatment  2   1.7623 0.46659 5.2483  0.001 ***
Residual  12   2.0148 0.53341                  
Total     14   3.7771 1.00000

```



