
1. Prepare design file and group lists for ``make.lefse``

    Meta data link[strigolactone_meta_up.txt](https://github.com/liufangbaishikele/Soybean_rhizosphere_microbiome/blob/master/Differential%20abundance%20analysis/More_detailed_Lefse/strigolactone_meta_up.txt)
    ```
    grep -E 'CT|D14_RNAi' strigolactone_meta_up.txt | grep -v 'CT_GR24' |awk '{print $1 "\t" $6}'> D14_RNAi
    echo -e 'group\tTreatment' > D14_RNAi.design && cat D14_RNAi >> D14_RNAi.design
    awk '{print $1}' D14_RNAi.design | tail -16 | paste -s -d- > D14_RNAi_list
    ```
2. Generate input file for Lefse using mothur - ``make.lefse`` command. Design file could be found [here]()

    ``` make.lefse(shared=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.pick.0.03.subsample.shared,constaxonomy=strigolactone.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.cons.taxonomy,design=D14_RNAi.design,groups=CT_1-CT_2-CT_3-CT_4-CT_5-CT_6-CT_7-CT_8-G1_1_RNAi-G1_2_RNAi-G1_4_RNAi-G1_5_RNAi-G1_6_RNAi-G1_7_RNAi-G1_8_RNAi,scale=totalgroup) 
    ```
