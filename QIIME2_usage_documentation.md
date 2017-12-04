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
Here are the output
```
QIIME is caching your current deployment for improved performance. This may take a few moments and should only happen once per deployment.
Usage: qiime [OPTIONS] COMMAND [ARGS]...

  QIIME 2 command-line interface (q2cli)
  --------------------------------------

  To get help with QIIME 2, visit https://qiime2.org.

  To enable tab completion in Bash, run the following command or add it to
  your .bashrc/.bash_profile:

      source tab-qiime

  To enable tab completion in ZSH, run the following commands or add them to
  your .zshrc:

      autoload bashcompinit && bashcompinit && source tab-qiime

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  info                Display information about current deployment.
  tools               Tools for working with QIIME 2 files.
  dev                 Utilities for developers and advanced users.
  alignment           Plugin for generating and manipulating alignments.
  composition         Plugin for compositional data analysis.
  dada2               Plugin for sequence quality control with DADA2.
  deblur              Plugin for sequence quality control with Deblur.
  demux               Plugin for demultiplexing & viewing sequence quality.
  diversity           Plugin for exploring community diversity.
  emperor             Plugin for ordination plotting with Emperor.
  feature-classifier  Plugin for taxonomic classification.
  feature-table       Plugin for working with sample by feature tables.
  gneiss              Plugin for building compositional models.
  longitudinal        Plugin for paired sample and time series analyses.
  metadata            Plugin for working with Metadata.
  phylogeny           Plugin for generating and manipulating phylogenies.
  quality-control     Plugin for quality control of feature and sequence data.
  quality-filter      Plugin for PHRED-based filtering and trimming.
  sample-classifier   Plugin for machine learning prediction of sample
                      metadata.
  taxa                Plugin for working with feature taxonomy annotations.
  vsearch             Plugin for clustering and dereplicating with vsearch.

```


## To enable tab completion of command in bash by ``source tab-qiime``





