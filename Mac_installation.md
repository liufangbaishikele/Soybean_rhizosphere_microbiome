## Installation of various softwares in my Mac


1. Mac installation : python 3.5.2 > anaconda > qiime2-2017.12

  * Install python 3.5.2 from this [link](https://www.python.org/ftp/python/3.5.2/python-3.5.2-macosx10.6.pkg) and double     click to install
  * Download Anaconda package for Mac OS system from this [link](https://repo.anaconda.com/archive/Anaconda3-5.1.0-MacOSX-x86_64.pkg) and double-click to install
  * Follow this documentation for installation of qiime2 within conda environment [link](https://docs.qiime2.org/2018.4/install/native/#install-qiime-2-within-a-conda-environment)
  ```
   wget https://data.qiime2.org/distro/core/qiime2-2018.4-py35-osx-conda.yml
   conda env create -n qiime2-2018.4 --file qiime2-2018.4-py35-osx-conda.yml
   # OPTIONAL CLEANUP
   rm qiime2-2018.4-py35-osx-conda.yml
  ```
  * conda install wget
  * wget https://data.qiime2.org/distro/core/qiime2-2017.12-py35-osx-conda.yml
 
  * conda env create -n qiime2-2017.12 --file qiime2-2017.12-py35-osx-conda.yml

2. One day when I use conda to install jp (jason partialization) through anaconda, it showed me 
  ```
  ImportError: No module named conda.cli
  ```
  Here is the solution for it:
  ```
  export PYTHONPATH="/Users/fangliu/anaconda/lib/python3.5/site-packages:$PYTHONPATH"
  ```
