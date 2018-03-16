## Installation of various softwares in my Mac


1. Mac installation : python 3.5.2 > anaconda > qiime2-2017.12
Install python 3.5.2 from this link
Download Anaconda package for Mac OS system from this link and double-click to install
Follow this documentation for installation of qiime2 within conda environment link
conda install wget

wget https://data.qiime2.org/distro/core/qiime2-2017.12-py35-osx-conda.yml
conda env create -n qiime2-2017.12 --file qiime2-2017.12-py35-osx-conda.yml

2. One day when I use conda to install jp (jason partialization) through anaconda, it showed me 
```
ImportError: No module named conda.cli
```
Here is the solution for it:
```
export PYTHONPATH="/Users/fangliu/anaconda/lib/python3.5/site-packages:$PYTHONPATH"
```
