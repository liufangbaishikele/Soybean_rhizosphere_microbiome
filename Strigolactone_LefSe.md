#                               2016 strigolactone LefSe analysis documentation
-----------------

## Data manipulation within mothur -
The reason that I manipulate shared file inside of mothur instead of within R is that make.lefse command need .shared file as input.

### Family level shared file -- taxlevel=2

1. Rarefy\_based normalization - Subset ``.shared`` using rarefaction method to minimum read depth along all samples 

``sub_sample``
```

```
