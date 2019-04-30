## Select Agrobacterium specific primer set to used to quantify Agrobacterium rhizogene strains across treatments

1. We have ``Agrobacterium tumefaciens C58`` - which all names as ``Agrobacterium fabrum str. C58``
2. We also have ``Agrobacterium tumerfaciens EHA105`` strain, which we did not find genome sequence from NCBI genome database
3. In addition, we also have ``Agrobacterium rhizogenes K599`` strain, which named as ``Agrobacterium rhizogenes strain NCPPB2659``
4.Strategy to target Agrobacterium k599 specific region:

  * Download genome sequence of both Agrobacterium tumefaciens c58 - named Af_C58.fasta and Agrobacterium rhizogenes NCPPB2659 and named as Ar_k599.fasta
  
  * Make blast database using Ar_k599.fasta
  
  * Do blastn and set 1e-5 as the threshold, format using 6
  
  * Download the blast results and sort based on query start and endo location.
  
  * Find all the large gaps from the hit results
  
  * And using our unpublished python script ``seq_substr.py``to subset the sequences based on the gap location
  
  * Then copy the subseted fasta and blast against ``nr`` databased and limited to Rhizobium/Agrobacterium taxon.
  
  * After coupel of blast, Find that the region from **431507** to **447981** have about 4000bp only hit Agrobacterium rhizogene k599 itself. 
  
  * So, based on the blast results, narrow the target region to 442507-446507, and subset the sequence to this region.
  
  * Using this new sequence, blast against all nr database. It turned out only hit Agrobacterium k599 itself.
  
5. Now upload the target sequence to primer3 plus to design primer for qPCR

  * We want the amplicon length around 70-150bp
  * GC content 40-80%
  * primer length 18-22bp
  * Melting temperature around 58-60-62 celsius
  * Max polyx=3
  
6. At the same time, we need some positive control for qPCR. Based on 16S contigs count table, I could find contigs with consistent count across control and over-expression treatment. But when I extract the sequence and blast in NCBI, it hits different unculture bacteria. So, very unconfident to design primer based on 16S V3-V4 region as it could amplify numerous of bacteria with similar V3-V4 regions. And from gel band, we could not tell the difference between sequences with same length but different allels.

7. I was thinking to add a known bacteria strain with specific and equal concentration. We can design this strain specific primer set as positive control. -- This method are called using spiked known bacteria strain as internal standard.

**Internal standard issues**

 * Some researchers using IAC (internal amplification control, but the purpose is to evaluate PCR inhibation to help confirm true negative amplification instead of PCR inhibation. OR [exogenous internal positive control (IPC)](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2672.2009.04364.x) 
 
 * There is no natural internal standard sequence that can be used in analogy to housekeeping genes widely used for pure culture or tissue analysis. More often, researchers using specific bacteria strain carrying specific construct are used as extraction standard to correct for extraction losses. Lebuhn et al, compared samples spiked and non-spiked with the target sequence to achieve this correction for their qPCR and RTqPCR data.[Reference](https://www.sciencedirect.com/science/article/pii/S0167701206003083#bib30)

 * Alternatively, we could use 16S rRNA gene as internal standard with the assumption that the recovery efficiency for each of the target gene present in a soil samples is identical. 

 * For a higher accuracy of the negative control, it is favourable to add an internal control before the DNA extraction. In this way, variations in the sample treatment (DNA extraction and PCR) will be reflected in the IPC result. This will allow further normalization of the resulting data and obtain a more accurate quantification of gene copy number in soil samples. An internal standard prepared from plasmid DNA can be added to soil samples as a defined amount of e.g. Escherichia coli harbouring the plasmid or as purified DNA (Lebuhn et al. 2004; Mumy and Findlay 2004; Park and Crowley 2005).[reference](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2672.2009.04364.x)
 

8. For my study, I have already done with the DNA extraction, so, I can't use Plasmid as internal standard to role out the extraction effeciency bias or PCR efficiency. 

9. My final descision is to use 16S universal primer as a 'housekeeping' standard to quantify the relative abundance of Agrobacterium k599 within each samples. [Reference](https://academic.oup.com/femsec/article/58/3/572/523809)

* Primer F - CGGCAACGAGCGCAACCC
* Primer R - CCATTGTAGCACGTGTGTAGCC






