# appreci8R - A Pipeline for PREcise variant Calling Integrating 8 tools

<p align="center">
    <img height="400" src="https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/34/24/10.1093_bioinformatics_bty518/2/bioinformatics_34_24_4205_f12.jpeg?Expires=1737544357&Signature=1Dqwvr9Z9q0LejOri7Aa3LY11uyUA~kJrD-ogwIefmVEbuVckUAwZRS9cZ~aMx2IDuZo286zhgiUW5edVjZUsmudm0Fy52cF8WfvIMx4HEHlVCt4tMoIkQAn8lXeZZ9iy8MqF3JHCJub9cyW4-O8OxMu4DUhxqwmtvICtnmxIkxYmz66U8kjzdaLSK~lOiIyixXVZJNgDHN8sF~tUo59jwIrdnVcmg-8CVZkNX5-pb1E19qxpHe8soe00NYmO-~XpCjmgt6tldAos2fV-aO4vfjrexk~1UrU~HnL26Laa4APzsHtGaX3~d-e6MTGW-oQN593mKtN2pJOv9GrW46M4g__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA">
</p>

<b>Fig. 1.</b>  Overview of the analysis performed by appreci8

<p align="center">
    <img height="400" src="https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/34/24/10.1093_bioinformatics_bty518/2/bioinformatics_34_24_4205_f13.jpeg?Expires=1737544363&Signature=tMAw449ixwI6cT0qe0ZSMrEcggX2nNMd12oJHWsQjkXlZHY~WxB-guep2NKtwlR5yaRcY9IkXqtAqDzcDZgcfNA~Aed7zxwMRiIS1EHMLEuRgSTg7H3YGnTb6jkyIZSGBJRvr0aVirHCC9othDnqZSw8Zs11FRjzx5lXOIHomONeYfdPt7MCkqzHBQpb8Y~NOnbHg01zBCpG3Awyfgf6Seom6j4WpoE8del6tWYO91Y~0WtpLLm9jW--~11pUwfRe0an3m8LVFAGgJHjzplOjaJxzk75x-FZtl516DGCLMj31ILbVV9EgrYQ35vccwDIBhPwYk-n39tuZ-cB3RtXgw__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA">
</p>

<b>Fig. 2.</b> General principle of filtration with appreci8. Calls are classified as ‘Mutations’, ‘Polymorphism’ or ‘Artifact’ on the basis of an artifact- and a polymorphism score


For the use of next-generation sequencing in clinical routine valid variant calling results are crucial. However, numerous variant calling tools are available. These tools usually differ in the variant calling algorithsms, the characteristics reported along with the varaint calls, the recommended filtration strategies for the raw calls and thus, also in the output. Especially when calling variants with a low variant allele frequency (VAF), perfect results are hard to obtain. High sensitivity is usually accompanied by low positive predictive value (PPV).

Appreci8 is a variant calling pipeline for detecting single nucleotide variants (SNVs) and short indels (up to ~30 bp) in next-generation sequencing (NGS) data. By integrating and filtering the output of eight individual variant calling tools (GATK, Platypus, VarScan, Freebayes, LoFreq, SAMtools, SNVer and VarDict; for considering your own selection of variant calling tools, please have a look at the R package 'appreci8R', which allows you to evaluate an unlimited individual selection of tools, for the user interface version aof the appreci8R it is limited to 13) on the basis of an artifact- and a polymorphism score, appreci8 succeeds in calling variants with high sensitivity and positive predictive value even at variant allele frequencies of 1%.



Important note: Currently, only hg19 is supported.

Sandmann S, Karimi M, de Graaf AO, Rohde C, Göllner S, Varghese J, Ernsting J, Walldin G, van der Reijden BA, Müller-Tidow C, Malcovati L, Hellström-Lindberg E, Jansen JH, Dugas M. appreci8: a pipeline for precise variant calling integrating 8 tools. Bioinformatics. 2018 Dec 15;34(24):4205-4212. doi: 10.1093/bioinformatics/bty518. PMID: 29945233; PMCID: PMC6289140.



## Requirements
appreci8 is a collection of bash and R scripts. It was developed on Linux Ubuntu. For running the R scripts, you need Version 4.1.0 or higher.

## Prerequisites
To run appreic8, you need to

* Download and install the 8 variant calling tools GATK, Platypus, VarScan, Freebayes, LoFreq, SAMtools, SNVer and VarDict. Select the version you like best. Consider the Config-file to set the correct paths to the tools. Carefully check the execution scripts 01 to 08 and adapt the default configuration if required.
* Download all dependencies at https://uni-muenster.sciebo.de/s/TzufkAEK1stdGHv⁠
* Prepare the data as follows (as an example, check out our Example-folder):
    * SampleNames.txt: The names of the samples you wish to analyze (without file extension, one name per line)
    * vcf_header.txt: Standard vcf file header (available in the appreci8 folder)
    * Folder alignment: Containing the bam- and bai files of the samples you wish to analyze (format: sample1.bam, sample1.bai etc.)
    * Folder snpEff_ann: Hotspots.txt: A list containing known hotspot mutations, covering Gene, Mutation (change on amino acid level, one-letter-code), Min_VAF (minimum allelic frequency at which you expect these mutations); an empty list can be passed, containing the header and three NA's (available in the appreci8 folder)
    * Folder targetRegions: targetRegions.bed: Bed file containing the target regions to be analyzed (no header, no information except for chr, start, end; 1 instead of chr1 etc.; for an example see file in the Example folder)

## Configuration
Appreci8 performs filtration of the data based on the following parameters (feel free to change dependent on your individual requirements):
* min_alt: minimum number of reads with the alternate allele (default: 20)
* min_dp: minimum depth (default: 50)
* min_vaf: minimum variant allele frequency (default: 0.01; do not choose values below 0.01)
* min_bq_alt: minimum mean base quality for reads with the alternate allele (default: 15)
* max_bq_diff: maximum difference for 'mean base quality reference' - 'mean base quality alternative' (default: 7)
* max_samples: maximum number of samples that are allowed to feature the same variant without penalizing (default: 3; if your data set contains more than 3 replicates of the same sample, it is recommended to increase this value)
* provean_del: minimum provean score to classify a variant as deleterious (default: 3)
* provean_tol: maximum provean score to classify a variant as tolerated (default: 1.5)
* primer: bed file containing primer locations is provided within the documents-folder (default: "FALSE")

Please note:
* currently only hg19 is supported
* the license (MIT) just covers the appreci8-scripts provided in this repository. The individual variant calling tools have their own licenses (special focus to GATK as e.g. GATK 3.3.0 is licensed by the Broad Institute and is made available for free to academic users)


## Detailed documentation
For detailed documentation on the performance of appreci8, please consider our publication (https://doi.org/10.1093/bioinformatics/bty518).

In case of errors or feature requests, do not hesitate to open an issue or contact Sarah Sandmann (sarah.sandmann@uni-muenster.de).


