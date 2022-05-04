# SB Digestor
SB Digestor: a tailored driver genes identification tool for dissect heterogenous of Sleeping Beauty induced tumors
The goal of SB Digestor is to digest tumor heterogeneity to improve identify driver genes of Sleepling Beauty tumors.In our computational algorithm, we not only took statistical significance into account just like other tools did, but al so we consider the biological aspects like tumor heterogeneity, which may cause missing driver gene identification. We firstly conducted the saturation analysis individually to identify the inter-tumor heterogeneity, then based on which we calculated the tailored driver genes identification parameter for further data processing. As each tumor analyzed separately, so we could get much stable results, no matter how many sequenced reads and how large scale of samples. The application of saturation analysis could precisely reflect the correlation between the reads number and identified driver gene number, therefore it benefits the data depth cut-off threshold selection.

<img width="1651" alt="Screenshot 2022-05-04 at 17 02 45" src="https://user-images.githubusercontent.com/66343257/166651959-429910b4-0f82-4f4c-847b-73b44e8d80e3.png">

In this pipeline, the raw sequences were processed by removing SB transposon sequence, DNA barcode tags, and sequence beyond the restriction enzyme site. Then we discarded the processed sequences that shorter than 20bp since they are prone to map to multiple genomic locations. Following the filtering process, we determined the reliable transposon insertions based on both statistical enrichment and biological heterogeneity. For statistical part, in order to figure out whether a gene got more insertions than expected probability, we performed a binomial test for each gene for each sample. Thus, we firstly computed TA number for each gene and the whole genome, then calculated an expected insertion probability for each gene. Simultaneously, we used the remained sequences aligned to the reference genome(mm10) by bowtie2. Then we annotated each gene by counting reads number for each sample as the observed reads number. Followed that, we did the binomial test to get the significant insertional genes for each sample, which we regard it as the statistics library. For biological heterogeneity part, we tried to define the relationship between the sequenced reads number and the detectable gene number under the tumor heterogeneity circumstances. Here, we randomly extracted reads for 50 gradients for each sample. For each gradient of reads, we mapped and annotated them one by one. For a certain sample, if an annotated gene exists in the previous statistics library, we considered it as a reliable insertion gene. Then we fitted curves by the gradient reads number and gradient gene number for each sample. Here, we obtained the adapted formula, and we defined the relationship reads number and gene number. Therefore, we could deduce depth cutoff calculation formula for each sample (Depth=reads number/gene number). Following that is the candidate driver gene determination step through the read depth cutoff. Then the common insertion genes were summarized within all the tumor samples. Last, we predicted tumor suppressors and oncogenes by the SB insertion directions.





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
