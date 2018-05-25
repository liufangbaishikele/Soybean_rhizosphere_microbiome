
1. Prepare design file and group lists for ``make.lefse``

    Meta data link [strigolactone_meta_up.txt](https://github.com/liufangbaishikele/Soybean_rhizosphere_microbiome/blob/master/Differential%20abundance%20analysis/More_detailed_Lefse/strigolactone_meta_up.txt)
    ```
    grep -E 'CT|D14_RNAi' strigolactone_meta_up.txt | grep -v 'CT_GR24' |awk '{print $1 "\t" $6}'> D14_RNAi
    echo -e 'group\tTreatment' > D14_RNAi.design && cat D14_RNAi >> D14_RNAi.design
    awk '{print $1}' D14_RNAi.design | tail -16 | paste -s -d- > D14_RNAi_list
    ```
2. Generate input file for Lefse using mothur - ``make.lefse`` command. Design file could be found [here](https://github.com/liufangbaishikele/Soybean_rhizosphere_microbiome/blob/master/Differential%20abundance%20analysis/More_detailed_Lefse/D14_RNAi.design)

    ``` make.lefse(shared=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.pick.0.03.subsample.shared,constaxonomy=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.cons.taxonomy,design=D14_RNAi.design,groups=CT_1-CT_2-CT_3-CT_4-CT_5-CT_6-CT_7-CT_8-G1_1_RNAi-G1_2_RNAi-G1_4_RNAi-G1_5_RNAi-G1_6_RNAi-G1_7_RNAi-G1_8_RNAi,scale=totalgroup) 
    ```
3. Running Lefse commands (in terms of installation and parameter details please found in the strigolactone_lefse_documentation)

    ```
    lefse-format_input.py strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.pick.0.03.subsample.0.03.lefse strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.pick.0.03.subsample.0.03.lefse_in  -c 1 -s -1 -u 2 -o 1000000

    run_lefse.py strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.pick.0.03.subsample.0.03.lefse_in strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.pick.0.03.subsample.0.03.lefse_res_LDA4 -l 4 -b 200 -t Strigolactone_Control_vs_D14_RNAi_OTU_lefse_LDA4 -y 0

    lefse-plot_res.py strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.pick.0.03.subsample.0.03.lefse_res_LDA4 strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.pick.0.03.subsample.0.03.lefse_res_LDA4.png --dpi 300 --title Strigolactone_Control_vs_D14_RNAi_LefSE_LDA4.png

    run_lefse.py strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.pick.0.03.subsample.0.03.lefse_in strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.pick.0.03.subsample.0.03.lefse_res_LDA2 -l 2 -b 200 -t Strigolactone_Control_vs_D14_RNAi_OTU_lefse_LDA2 -y 0

    lefse-plot_res.py strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.pick.0.03.subsample.0.03.lefse_res_LDA2 strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.pick.0.03.subsample.0.03.lefse_res_LDA2.png --dpi 300 --title Strigolactone_Control_vs_D14_RNAi_LefSE_LDA2.png

    ```
    
 4. Lefse output filtering based on analysis purpose
 
     ```
     grep "CT" strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.pick.0.03.subsample.0.03.lefse_res_LDA2 > Control_vs_D14_RNAi_LDA2
    grep "D14_RNAi" strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.pick.0.03.subsample.0.03.lefse_res_LDA2 >> Control_vs_D14_RNAi_LDA2
    grep 'Otu.....' Control_vs_D14_RNAi_LDA2 > Control_vs_D14_RNAi_LDA2_OTU
    grep -v 'Otu.....' Control_vs_D14_RNAi_LDA2 > Control_vs_D14_RNAi_LDA2_non_OTUs
    grep -v "unclassified" Control_vs_D14_RNAi_LDA2_non_OTUs > Control_vs_D14_RNAi_LDA2_non_OTU_classified
    sort -k1,1 -k2,2 -k3,3 -k4,4 -k5,5 -k6,6 Control_vs_D14_RNAi_LDA2_non_OTU_classified > Control_vs_D14_RNAi_LDA2_non_OTU_classified_sort
    grep -E -v 'Otu.....\.' Control_vs_D14_RNAi_LDA2 | sort -k1,1 -k2,2 -k3,3 -k4,4 -k5,5 -k6,6 > Control_vs_D14_RNAi_LDA2_up
    lefse-plot_res.py Control_vs_D14_RNAi_LDA2_non_OTU_classified_sort Control_vs_D14_RNAi_LDA2_non_OTU_classified_sort.png --title 'Control vs D14_RNAi  lefse LDA>2' --title_font_size 15 --class_legend_font_size 10 --max_feature_len 50 --dpi 300
     ```
