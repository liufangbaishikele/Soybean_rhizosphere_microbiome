## ITS2 sequencing PNA design for soybean



As our endosphere and seed ITS sequencing got high percentage of Glycine max ITS contamination. we decided to design ITS2 blocker to elimimate the amplification of soybean ITS amplification

---

Principles [Lundberg](https://www.nature.com/articles/nmeth.2634#supplementary-information):
Peptide nucleic acid (PNA) design.
To identify candidate PNA oligo sequences, we fragmented in silico the full length A. thaliana plastid and mitochondrial 16S sequences into short k-mers for k of length 9, 10, 11, 12 and 13, and we queried for exact matches against the 4 February 2011 version of the Greengenes 16S training set comprising 35,430 unique, high-quality full-length bacterial sequences (Supplementary Fig. 10). A. thaliana–specific k-mers falling between the 515F and 806R 16S rRNA primers (V4 region) were considered candidates and were lengthened as necessary to increase the predicted melting temperatures and were screened for design characteristics15, 20.

A successful elongation arrest PNA clamp is generally between 13 bp and 17 bp and has an annealing temperature above that of the PCR primer whose extension it blocks and a melting temperature above that used for the extension cycle20. We designed 17-mer sequences to block the plastid and mitochondria, each with a predicted melting temperature around 80 °C (Supplementary Table 1f). Melting temperature, problematic hairpins, GC content and other design considerations were calculated using the Life Technologies PNA designer (http://www6.appliedbiosystems.com/support/pnadesigner.cfm).

The anti-mitochondrial PNA (mPNA) 5′-GGCAAGTGTTCTTCGGA-3′ and the anti-plastid PNA (pPNA) 5′-GGCTCAACCCTGGACAG-3′ (Supplementary Table 1f) were ordered from PNA Bio. Lyophilized PNA was resuspended in sterile water to a stock concentration of 100 μM. For PNA concentrations that were repeatedly tested, working stocks of 5 μM, 15 μM, 25 μM and 40 μM were prepared in water. All stocks were stored at −20 °C and heated to 65 °C before use to resolubilize any precipitate.
PNA should have higher annealing temperature

## Procedure

1. I went back to my ITS sequence pipeline and extracted non-funfi contigs and their count
2. Sorted the count table and found one dominant contigs - M04398_3_000000000-C4DLV_1_1106_23817_21172
3. Extracted the fasta of this contigs
4. BLAST in NCBI against nr database. IT hits 100% identity to Glycine max - `Glycine max 18S ribosomal RNA gene, partial sequence; internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence; and 26S ribosomal RNA gene, partial sequence`
5. BLAST this contigs agaist Glycine genomes. It hit Glycine max with 100% identity and 99% identity to Glycine soja
6. Fragement this contig to kmers with length range from 9, 10, 11, 12 and 13.
7. BLAST these kmers to ITS UNITE database. After it, Kmers that have no hit to the UNITE database will be used as the candicate for PNA design.

## PNA screening criteria

* Length between 13bp to 17bp
* Has higher anealing temperature
* A meltering temperature above that used for extension cycle20.
* Melting temperature, problematic hairpins, GC content and other design considerations were calculated using the Life Technologies PNA designer (http://www6.appliedbiosystems.com/support/pnadesigner.cfm).
* PNA will be ordered from PNA Bio company.

## Practical prosedure
* Generated Kmers with length of 9,10,11,12,and 13bp based on M04398_3_000000000-C4DLV_1_1106_23817_21172 using Jellyfish software.

```
/staton/software/jellyfish-2.2.10/jellyfish count -m 9 -s 100M -t 10 -C M04398_3_000000000-C4DLV_1_1106_23817_21172
mv mer_counts.jf 9mer_counts.jf
/staton/software/jellyfish-2.2.10/jellyfish dump 9mer_counts.jf > 9_mer_counts_dumps.fa

/staton/software/jellyfish-2.2.10/jellyfish count -m 10 -s 100M -t 10 -C M04398_3_000000000-C4DLV_1_1106_23817_21172
mv mer_counts.jf 10mer_counts.jf
/staton/software/jellyfish-2.2.10/jellyfish dump 10mer_counts.jf > 10_mer_counts_dumps.fa

/staton/software/jellyfish-2.2.10/jellyfish count -m 11 -s 100M -t 10 -C M04398_3_000000000-C4DLV_1_1106_23817_21172
mv mer_counts.jf 11mer_counts.jf
/staton/software/jellyfish-2.2.10/jellyfish dump 11mer_counts.jf > 11_mer_counts_dumps.fa

/staton/software/jellyfish-2.2.10/jellyfish count -m 12 -s 100M -t 10 -C M04398_3_000000000-C4DLV_1_1106_23817_21172
mv mer_counts.jf 12mer_counts.jf
/staton/software/jellyfish-2.2.10/jellyfish dump 12mer_counts.jf > 12_mer_counts_dumps.fa

/staton/software/jellyfish-2.2.10/jellyfish count -m 13 -s 100M -t 10 -C M04398_3_000000000-C4DLV_1_1106_23817_21172
mv mer_counts.jf 13mer_counts.jf
/staton/software/jellyfish-2.2.10/jellyfish dump 13mer_counts.jf > 13_mer_counts_dumps.fa

for file in *.fa; do echo $file; grep -v '>' $file > update_$file;done
for file in update_*; do echo $file; awk '{print ">" NR; print $0}' $file > header_$file; done
```
* Concatinate all the Kmers into one file called 9to13mer.fa and remove all header information and add new header information
```
cat *.fa > 9to13mer_counts_dumps.fa
grep -v '>' 9to13mer_counts_dumps.fa > 9to13mer_seq_only
for seq in $(cat 9to13mer_seq_only); do echo $seq >> match_out; grep -c $seq UNITEv8_sh_97_s_k__Fungi.fasta >> match_out; done
grep -B 1 '^0' match_out > Glycine_max_ITS_specific_mers
#Edit this Glycine_max_ITS_specific_mers file to remove all --lines and 0 lines
awk '{print ">" NR; print $0}' Glycine_max_ITS_specific_mers_edit > Glycine_max_ITS_specific_mers_edit_add_header
```
* Modified ``UNITEv8_sh_97_s_all.fasta`` reference to exclude non-Fungi and generated ``UNITEv8_sh_97_s_k__Fungi.fasta``

```
grep 'k__Fungi' UNITEv8_sh_97_s_all.tax > UNITEv8_sh_97_s_k__Fungi.tax
awk '{print $1}' UNITEv8_sh_97_s_k__Fungi.tax > k__Fungi_seqID
grep -Ff k__Fungi_seqID UNITEv8_sh_97_s_all.fasta > UNITEv8_sh_97_s_k__Fungi.fasta
```
* Tried to use blast to find Glycine max specific Kmers by blast against UNITEv8_sh_97_s_k__Fungi.fasta. It turned out no hit, but I am sure there are lots of exact match. The problems is blast take input sequence larger than 20bp???
* Trying to using mapping tools, bowtie. However my reference are thousands of ITS sequence. Which make this more complicated even after I got results.
* So, I came up with shell script to do exact match and extract those that do not match.

```
grep -v '>' 9to13mer_counts_dumps.fa > 9to13mer_seq_only
for seq in $(cat 9to13mer_seq_only); do echo $seq >> match_out; grep -c $seq UNITEv8_sh_97_s_k__Fungi.fasta >> match_out; done
grep -B 1 '^0' match_out > Glycine_max_ITS_specific_mers
#Edit this Glycine_max_ITS_specific_mers file to remove all --lines and 0 lines
awk '{print ">" NR; print $0}' Glycine_max_ITS_specific_mers_edit > Glycine_max_ITS_specific_mers_edit_add_header
```
* After extracted all Glycine max specific Kmers. Mapped to M04398_3_000000000-C4DLV_1_1106_23817_21172 (Glycine max ITS2). 

*After I done with the mapping using bowtie and visualization, I checked the GC content and try to design PNA based on requirements, I realized that those Glycine max ITS2 unique kmers can not match with the Glycine_max_ITS2 sequence. Then I went back the kmer file and found that the kmers were extracted with both forward and reverse direction.* 

*So, I reanalyzed/updated the data from ``Glycine_max_ITS_specific_mers_edit`` to remove those reverse complement kmers that did not match fungi ITS database. Basically, I removed all kmers that do not match with the Glycine max ITS2 sequneces. Therefore, end with ``Glycine_max_ITS_specific_mers_edit_match_seq_edit`` as the Glycine specific kmers and will be used for mapping against Glycine max ITS reference to look at the location of those kmers.*

      ```
      for seq in $(cat Glycine_max_ITS_specific_mer) ; do echo $seq; echo $seq >> Glycine_max_ITS_specific_mers_edit_match_count; grep -c $seq M04398_3_000000000-C4DLV_1_1106_23817_21172>> Glycine_max_ITS_specific_mers_edit_match_count; done
      ```
**Here to remove kmers got no hit in Glycine max ITS sequence**

      ```
      grep -B 1 '1' Glycine_max_ITS_specific_mers_edit_match_count > Glycine_max_ITS_specific_mers_edit_match_seq
      ```
**Edit this seq file to make it cleaser and add header to this seq file**

      ```
      awk '{print ">" NR; print $0}' Glycine_max_ITS_specific_mers_edit_match_seq_edit > Glycine_max_ITS_specific_mers_edit_match_seq_edit_with_header
      ```

**Doing mapping using bowtie2**

  1. Creat index file

      ```
      /staton/software/bowtie2-2.3.4.3-linux-x86_64/bowtie2-build -f M04398_3_000000000-C4DLV_1_1106_23817_21172.fasta ref
      ```

  2. Run alignment

      ```
      /staton/software/bowtie2-2.3.4.3-linux-x86_64/bowtie2  -x ref -f Glycine_max_ITS_specific_mers_edit_match_seq_edit_with_header -N 0 --end-to-end --norc -S Kmer_vs_Glycine_max_ITS2.sam
      ```

   3. convert sam file to bam file

      ```
      /staton/software/samtools-1.9/samtools view Kmer_vs_Glycine_max_ITS2.sam -b  -o Kmer_vs_Glycine_max_ITS2.bam
      ```

   4. Sort bam file inder to generate IGV readable index file and bam file

      ```
      /staton/software/samtools-1.9/samtools sort Kmer_vs_Glycine_max_ITS2.bam -o Kmer_vs_Glycine_max_ITS2_sorted.bam
      ```

   5. Generate the index file

      ```
      /staton/software/samtools-1.9/samtools index -b Kmer_vs_Glycine_max_ITS2_sorted.bam 
      ```

   6. secure copy sorted bam and sorted_bam.bai index file as well as reference fasta file to local computer. Using IGV for visualization

   7. Next step will be look into the distribution of those Kmers and found a region with length around 20bp and design PNA based on that.

   8. Extract one of the fragement cover overlapped  kmers with length of 18bp. 

   9. Double check its uniqueness against all fungi ITS database

  10. Then I need to check the annealing temperature for all the candicate PNAs
  

## Now, I need to stop here and continue my committee meeting ppt.

```
https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0413-8#Decs
https://www.ncbi.nlm.nih.gov/nuccore/NC_008483.2?report=fasta
https://www.ncbi.nlm.nih.gov/nuccore/NC_008483.2?report=fasta
Populus PNA: CGAGGGCACGTCTGCCTGG 
NCBI-BLAST using PNA against populus genome, found 100 identical hit.
https://www.nature.com/articles/nmeth.2634#supplementary-information
https://media.nature.com/original/nature-assets/nmeth/journal/v10/n10/extref/nmeth.2634-S1.pdf
https://www.pnabio.com/
https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0025715
```
