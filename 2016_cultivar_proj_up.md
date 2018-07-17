----
Why update the analysis of the data?
----

* Want to update the plots and make all plots using same color defination
* Want to change CV# to corresponding cultivar abbreviation easier to remember and make more sense.

---
How to update the analysis 
---

  1. First changed sampleID inside of shared file in beacon server by vim
  2. Change cultivar_meta file in local computer to match with shared file

---
OTU based analysis
---
* For OTU based analysis the OTU table were normalized at 19223 depth.

---
Genus based analysis
---

    1. The OTU table were converted to phylotypes based on taxonomy information 
    2. Any genus with sum smaller than 50 were removed before Lefse analysis

  
