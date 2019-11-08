LefSe based differential abundance analysis

-- Usually, after a global glance of community composition between treatment, we will have some idea how different the communities are between the treatments. Usually, in this step we want to project our wide data (which means tens of thousands OTUs across around 100 or 200 samples) onto a reduced dimention (usually 2 or 3 dimentions). By doing this, the variations between samples could be maximally projected from tons of dimentions to 2/3 dimentions. This process could be accomplished by PCoA, NMDS or CCA. For microbial community analysis, the PCoA visualization method is the mostly used based on my feeling.

-- PCoA with weighted Uni-Frac distance was suggested to be the best visualization strategy for microbial community data - 

* Now the next step is to focus on each taxa to see if they were differentially abundant between treatments or not.

**Potential choices for differential abundance analysis**

1. DESeq 
2. MetagenomeSeq
3. LefSE

**LefSe is the softwared I used for differential abundance analysis**

1. Generate the lefse input file generated from mothur using ``make.lefse``

  * Design file is needed for make.lefse
  ```
  group   Treat
104B_1w_ML      Bulk_1wk
104E_1w_ML      Endo_1wk
104R_1w_ML      Rhi_1wk
105B_1w_ML      Bulk_1wk
105E_1w_ML      Endo_1wk
105R_1w_ML      Rhi_1wk
108B_1w_ML      Bulk_1wk
108E_1w_ML      Endo_1wk
108R_1w_ML      Rhi_1wk
109B_1w_ML      Bulk_1wk
109E_1w_ML      Endo_1wk
109R_1w_ML      Rhi_1wk
201B_1w_ML      Bulk_1wk
201E_1w_ML      Endo_1wk
  ```
  
  * Generate the sample list for subset the dataset
  
  ```
  104B_1w_ML-104E_1w_ML-104R_1w_ML-105B_1w_ML-105E_1w_ML-105R_1w_ML-108B_1w_ML-108E_1w_ML-108R_1w_ML-109B_1w_ML-109E_1w_ML-109R_1w_ML-201B_1w_ML-201E_1w_ML-201R_1w_ML-202B_1w_ML-202E_1w_ML-202R_1w_ML-206B_1w_ML-206E_1w_ML-206R_1w_ML-207B_1w_ML-207E_1w_ML-207R_1w_ML-209B_1w_ML-209E_1w_ML-209R_1w_ML-301B_1w_ML-301E_1w_ML-301R_1w_ML-302B_1w_ML-302E_1w_ML-302R_1w_ML-306B_1w_ML-306E_1w_ML-309B_1w_ML-309E_1w_ML-309R_1w_ML-310B_1w_ML-310E_1w_ML-310R_1w_ML-411B_1w_ML-411E_1w_ML-411R_1w_ML
  ```
  
  
  2. Now process differential abundance analysis
  
  ```
  
  lefse-format_input.py Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.subsample.tx.1.pick.1.lefse Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.subsample.tx.1.pick.1.lefse_in  -c 1 -s -1 -u 2 -o 1000000

run_lefse.py Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.subsample.tx.1.pick.1.lefse_in Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.subsample.tx.1.pick.1.lefse_res_LDA4 -l 4 -b 200 -t Week1_Compartment_effects_lefse_LDA4 -y 0
  
  ```
  
  3. Now, you can look into the Lefse analysis results - table
  
  
  ```
  Bacteria.Acidobacteria  5.47792337114   Bulk_1wk        5.13101032874   5.00095789612e-09
Bacteria.Acidobacteria.Acidobacteriia   4.8756854107    Bulk_1wk        4.5305754404    6.48505048859e-09
Bacteria.Acidobacteria.Acidobacteriia.Solibacterales    4.5870311472    Bulk_1wk        4.23993016595   9.06368652001e-09
Bacteria.Acidobacteria.Acidobacteriia.Solibacterales.Solibacteraceae__Subgroup_3        4.5870311472    Bulk_1wk        4.23993016595   9.06368652001e-09
Bacteria.Acidobacteria.Blastocatellia__Subgroup_4       5.08527083787   Bulk_1wk        4.75030705452   5.00095789612e-09
Bacteria.Acidobacteria.Blastocatellia__Subgroup_4.Pyrinomonadales       4.98680941035   Bulk_1wk        4.65964501325   5.00095789612e-09
Bacteria.Acidobacteria.Blastocatellia__Subgroup_4.Pyrinomonadales.Pyrinomonadaceae      4.98680941035   Bulk_1wk        4.65964501325   5.00095789612e-09
Bacteria.Acidobacteria.Blastocatellia__Subgroup_4.Pyrinomonadales.Pyrinomonadaceae.RB41 4.98680941035   Bulk_1wk        4.65992496258   5.00095789612e-09
Bacteria.Acidobacteria.Subgroup_6       4.79560416635   Bulk_1wk        4.42335853541   7.06002856923e-09
Bacteria.Acidobacteria.Subgroup_6.Subgroup_6_or 4.79396329667   Bulk_1wk        4.42488531191   6.48505048859e-09
Bacteria.Acidobacteria.Subgroup_6.Subgroup_6_or.Subgroup_6_fa   4.79396329667   Bulk_1wk        4.42488531191   6.48505048859e-09
Bacteria.Acidobacteria.Subgroup_6.Subgroup_6_or.Subgroup_6_fa.Subgroup_6_ge     4.79396329667   Bulk_1wk        4.42488531191   6.48505048859e-09
  ```
  
  4. Visualized the results in a barplot
  
  ```
  lefse-plot_res.py Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.subsample.tx.1.pick.1.lefse_res_LDA4 Seed.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.subsample.tx.1.pick.1.lefse_res_LDA4.png --dpi 300 --title Week1_Compartment_effects_LefSe_LDA4
  ```
  
  
