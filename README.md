# SB Digestor
SB Digestor: a tailored driver genes identification tool for dissect heterogenous of Sleeping Beauty induced tumors
The goal of SB Digestor is to digest tumor heterogeneity to improve identify driver genes of Sleepling Beauty tumors.In our computational algorithm, we not only took statistical significance into account just like other tools did, but al so we consider the biological aspects like tumor heterogeneity, which may cause missing driver gene identification. We firstly conducted the saturation analysis individually to identify the inter-tumor heterogeneity, then based on which we calculated the tailored driver genes identification parameter for further data processing. As each tumor analyzed separately, so we could get much stable results, no matter how many sequenced reads and how large scale of samples. The application of saturation analysis could precisely reflect the correlation between the reads number and identified driver gene number, therefore it benefits the data depth cut-off threshold selection.

![Image](https://github.com/AipingZhang/SBDigestor/blob/main/Main%20pipeline.pdf)




Installation
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Download the bin directory and uncompress it.


Installation of other dependencies 
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cutadpt：https://cutadapt.readthedocs.io/en/stable/

bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

samtools: http://www.htslib.org/

seqtk: https://github.com/lh3/seqtk

For those above tools, please install or link them to the bin folder

Perl modules: 

Parallel::ForkManager;

Statistics::R.
