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
  
  
