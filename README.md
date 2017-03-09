# Soybean rhizosphere microbiome-soil type and cultivar impracts
Preliminary analysis for 16S seq data from Ag&amp;Forest cultivar project 
## sequence decompress
'''
cd /lustre/projects/staton/projects/soybean_strigolactones/16S_raw_fastq/soybean_strigolactone_16S_03_02_2017
for file in *; do echo $file; cd $file; cp * ../../fastq_gunzip/; cd ..; done

'''
