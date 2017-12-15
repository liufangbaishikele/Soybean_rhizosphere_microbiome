# Qiime2 
----------------
## Intall and get familiar with QIIME
-----------------
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

## knock around the documentations

**``qiime --help``**
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


## To enable tab completion of command in bash by 
**``source tab-qiime``**

## Have a look at qiime information by 
**``qiime info``**

```
System versions
Python version: 3.5.4
QIIME 2 release: 2017.11
QIIME 2 version: 2017.11.0
q2cli version: 2017.11.0

Installed plugins
alignment 2017.11.0
composition 2017.11.0
dada2 2017.11.0
deblur 2017.11.0
demux 2017.11.0
diversity 2017.11.0
emperor 2017.11.0
feature-classifier 2017.11.0
feature-table 2017.11.0
gneiss 2017.11.0
longitudinal 2017.11.0
metadata 2017.11.0
phylogeny 2017.11.0
quality-control 2017.11.0
quality-filter 2017.11.0
sample-classifier 2017.11.0
taxa 2017.11.0
types 2017.11.0
vsearch 2017.11.0

Application config directory
/nics/d/home/fliu21/.config/q2cli

Getting help
To get help with QIIME 2, visit https://qiime2.org

Citing QIIME 2
If you use QIIME 2 in any published work, you should cite QIIME 2 and the plugins that you used. To display the citations for QIIME 2 and all installed plugins, run:

  qiime info --citations

```
**qiime info --citation** very useful 

## QIIME command structure

```
qiime demux emp-single \
  --i-seqs emp-single-end-sequences.qza \
  --m-barcodes-file sample-metadata.tsv \
  --m-barcodes-category BarcodeSequence \
  --o-per-sample-sequences demux.qza
```
1. First give plugin name - **demux** (in this case)
emp-single is the commands from demux plugin

2. Commands from the plugin
IF I want to know all available commands from this plugin, I could just type in ``qiime demux emp-single`` and tap enter button.
```
Usage: qiime demux [OPTIONS] COMMAND [ARGS]...
Description: This QIIME 2 plugin supports demultiplexing of single-end and
  paired-end sequence reads and visualization of sequence quality
  information.

  Plugin website: https://github.com/qiime2/q2-demux

  Getting user support: Please post to the QIIME 2 forum for help with this
  plugin: https://forum.qiime2.org

  Citing this plugin: No citation available. Cite plugin website:
  https://github.com/qiime2/q2-demux

Options:
  --help  Show this message and exit.

Commands:
  emp-paired  Demultiplex paired-end sequence data generated with the EMP
              protocol.
  emp-single  Demultiplex sequence data generated with the EMP protocol.
  summarize   Summarize counts per sample.
```
3. options for each commands- could be found from help documents.
```
qiime demux emp_single
```
Then hit enter button.

```
Usage: qiime demux emp-single [OPTIONS]

  Demultiplex sequence data (i.e., map barcode reads to sample ids) for data
  generated with the Earth Microbiome Project (EMP) amplicon sequencing
  protocol. Details about this protocol can be found at
  http://www.earthmicrobiome.org/protocols-and-standards/

Options:
  --i-seqs PATH                   Artifact: EMPPairedEndSequences |
                                  EMPSingleEndSequences | RawSequences
                                  [required]
                                  The single-end sequences to be
                                  demultiplexed.
  --m-barcodes-file PATH          Metadata file or artifact viewable as
                                  metadata. This option may be supplied
                                  multiple times to merge metadata  [required]
  --m-barcodes-category TEXT      Category from metadata file or artifact
                                  viewable as metadata  [required]
                                  The sample
                                  metadata category listing the per-sample
                                  barcodes.
  --p-rev-comp-barcodes / --p-no-rev-comp-barcodes
                                  [default: False]
                                  If provided, the barcode
                                  sequence reads will be reverse complemented
                                  prior to demultiplexing.
  --p-rev-comp-mapping-barcodes / --p-no-rev-comp-mapping-barcodes
                                  [default: False]
                                  If provided, the barcode
                                  sequences in the sample metadata will be
                                  reverse complemented prior to
                                  demultiplexing.
  --o-per-sample-sequences PATH   Artifact: SampleData[SequencesWithQuality]
                                  [required if not passing --output-dir]
                                  The
                                  resulting demultiplexed sequences.
  --output-dir DIRECTORY          Output unspecified results to a directory
  --cmd-config PATH               Use config file for command options
  --verbose                       Display verbose output to stdout and/or
                                  stderr during execution of this action.
                                  [default: False]
  --quiet                         Silence output if execution is successful
                                  (silence is golden).  [default: False]
  --help                          Show this message and exit.
```
-------------------------------
## Start analysis using QIIME
-------------------------------
1. First import sequence data (if undemultiplexed, both fastq and barcode is needed) and sample data (if undemultiplexed, barcode catogery is needed)
2. Demultiplex (if already indivisual samples, skip this step)
3. Make a summary of the sequencing depth and quality of sequences
4. Denoise using ``dada2`` OR ``deblur``

**Web based dada2 documentation**

  * First, open [qiime documentation](https://docs.qiime2.org/2017.11/)
  * Second, click the ``Plugins`` menu on the left panel
  * Third: click [dada2: Plugin for sequence quality control with DADA2](https://docs.qiime2.org/2017.11/plugins/available/dada2/)
  
  * Now, I found two function of dada2
  ``(1) demoise-paired``
  ``(2) denoise-single``
  
**Documentation using command line**

  * Functions of dada2
  
```qiime dada2
```
  IT will give below information

  ```
  Usage: qiime dada2 [OPTIONS] COMMAND [ARGS]...

    Description: This QIIME 2 plugin wraps DADA2 and supports sequence quality
    control for single-end and paired-end reads using the DADA2 R library.

    Plugin website: http://benjjneb.github.io/dada2/

    Getting user support: Please post to the QIIME 2 forum for help with this
    plugin: https://forum.qiime2.org

    Citing this plugin: DADA2: High-resolution sample inference from Illumina
    amplicon data. Benjamin J Callahan, Paul J McMurdie, Michael J Rosen,
    Andrew W Han, Amy Jo A Johnson, Susan P Holmes. Nature Methods 13, 581–583
    (2016) doi:10.1038/nmeth.3869.

  Options:
    --help  Show this message and exit.

  Commands:
    denoise-paired  Denoise and dereplicate paired-end sequences
    denoise-single  Denoise and dereplicate single-end sequences
  ```
  * Look at options of **denoise-paired** by type in ``qiime dada2 denoise-paired``
  
  















