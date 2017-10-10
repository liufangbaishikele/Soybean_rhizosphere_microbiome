
## Following [nics data transfer documentation](https://www.nics.utk.edu/computing-resources/data-transfer)

1) Regitered globalbus account
2) In the NICS portal associate their NetID with their NICS account (see the image below) and
3) In the NICS portal setup their X.509 user certificate by associating their InCommon credential with their NICS account
4) compress all of the raw_fastq file and some big alignment file and transfer data from newton to nics

```
tar -zcvf filename.tar.gz orinigalfiles*
```


### Tried gsiscp method on command line, 

```
ssh fliu21@dtn1.newton.utk.edu 
gsiscp -c NONE -oNoneEnabled=yes -oNoneSwitch=yes /lustre/projects/staton/projects/soybean_strigolactone/folder datamover8.nics.utk.edu:/lustre/haven/gamma/staton/projects/01_first_run 
```
PASSCODE are always wrong.


### Applied Globus Web-based Data Transfer.

**Endpoint of Newton** I used is **dtn1**
**ACF Endpoint** is **nics\#ACF**
