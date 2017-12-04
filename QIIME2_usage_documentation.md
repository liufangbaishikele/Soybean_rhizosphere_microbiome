# Qiime2 
----------------

## Upload qiime to env

Module load qiime from ACF by ``module load qiime2/2017.11``
Check if qiime2 loaded to the current environment by ``module list`` (during this process anaconda3 is all loaded)
It shows that qiime2 is successfully loaded to working directory.
** -bash: qiime: command not found ** When I type in ``qiime``
Then I tried to load mothur ``module load mothur/1.39.5`` and checked using ``module list`` to make sure correctly uploaded.
When I type in ``mothur`` from my terminal, it poped out an error (mothur: symbol lookup error: /sw/cs400_centos7.3_acfsoftware/anaconda3/4.4.0/centos7.3_gnu6.3.0/anaconda3-4.4.0/lib/libreadline.so.6: undefined symbol: PC)
I thought the problem may caused by the path or environment. So I double checked the module list and found anaconda3 is uploaded. I guess
the problem may comes from this. So I unload anaconda3 by ``module unload anaconda3/4.4.0`` and type in mothur again. It worked!
BUT when I type ``qiime``, the command was still not found.

## Active conda environment

**NOW go back to qiime documentation** look carefully and found that after upload qiime, I still need to active the conda environment by
``source activate qiime2-2017.11`` 
now qiime become functional;)

## Start using qiime first by ``qiime --help``


