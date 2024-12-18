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

appreci8 is a pipeline for combining and filterating the output of 8 variant calling tools: GATK, Platypus, VarScan, Freebayes, LoFreq, SAMtools, SNVer and VarDict (for considering your own selection of variant calling tools, please have a look at the R package 'appreci8R', which allows you to evaluate an unlimited individual selection of tools, for the user interface version it is limited to 13). The final output contains a list of variant calls, classified as "probably true", "polymophism" or "artifact".

Important note: Currently, only hg19 is supported.

Sandmann S, Karimi M, de Graaf AO, Rohde C, Göllner S, Varghese J, Ernsting J, Walldin G, van der Reijden BA, Müller-Tidow C, Malcovati L, Hellström-Lindberg E, Jansen JH, Dugas M. appreci8: a pipeline for precise variant calling integrating 8 tools. Bioinformatics. 2018 Dec 15;34(24):4205-4212. doi: 10.1093/bioinformatics/bty518. PMID: 29945233; PMCID: PMC6289140.


