Astroid
=======

Isoform reconstruction and quantification package
This project is actively maintained at: http://www.netlab.uky.edu/p/bioinfo/
Please cite our paper: "Piecing the Puzzle Together: a Revisit to Transcript 
Reconstruction Problem in RNA-seq." BMC Bioinformatics, in press."

===============

[Yan Huang](http://protocols.netlab.uky.edu/~yan) \(yan at netlab dot uky dot edu\)

* * *

Table of Contents
-----------------

* [Introduction](#introduction)
* [Compilation & Installation](#compilation)
* [Usage](#usage)


* * *

## <a name="introduction"></a> Introduction

The advancement of RNA sequencing (RNA-seq) has provided an unprecedented opportunity
to assess both the diversity and quantity of transcript isoforms in an mRNA transcriptome. 
In this paper, we revisit the computational problem of transcript reconstruction and quantification. 
Unlike existing methods which focus on how to explain the exons and splice variants detected by the 
reads with a set of isoforms, we aim to reconstruct transcripts directly at the level of the reads, 
by piecing the reads into individual effective transcript copies. Simultaneously, the assembled 
effective copies provide an explicit measure of the quantity of each isoform, instead of estimating 
the transcript abundance solely based on the collective read count. We have developed a novel method 
named Astroid1 that solves the problem of effective copy reconstruction on the basis of a 
flow network. The RNA-seq reads are represented as vertices in the flow network and are 
connected by weighted edges that evaluate the likelihood of two reads originating from the 
same effective copy. A maximum likelihood set of transcript copies is then reconstructed by 
solving a minimum-cost flow problem on the flow network.
Simulation studies on the human transcriptome have demonstrated the superior sensitivity and 
specificity of Astroid in transcript reconstruction as well as improved accuracy in transcript 
quantification over several existing approaches. The application of Astroid on two real RNA-seq 
datasets has further demonstrated its accuracy through high correlation between the estimated 
isoform abundance and the qRT-PCR validations.

## <a name="compilation"></a> Compilation & Installation

To compile Astroid, simply run
   
    make

To install, simply put the Astroid directory in your environment's PATH
variable.

### Prerequisites

C++ is required to be installed. 


## <a name="usage"></a> Usage

Astroid samfilename outputfilefolder

The first argument is the input sam file name.
Sample sam file format:
HWI-ST254_0000:5:1:1575:1998#0/1 0 chr1 1263644 90 100M * 0 0 CAGTGCCCTCCATGCCCTGGCTGGCAGAAACCCTCAACAGCAGTCTGGGCACTGTGGGGCTCTCCCCGCCTCTCCTGCCTTGTTTGCCCCTCAGCGTGCC ccc\ccadbdgfbcfefdeedZcdYeb^ecaef_eaeea]cb^X\QWZ[[]_BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB NM:i:0 IH:i:1 HI:i:1
HWI-ST254_0000:5:1:1878:1970#0/1 0 chrX 70782818 136 50M112N50M * 0 0 NCCACACTTTTTTTATTGGTGATCATGCTAATATGTTCCCTCACCTGAAGAAAAAAGCAGTCATCGATTTTAAGTCCAATGGGCACATTTATGACAATCG BXSWTSUWXSccUYYYUUcVZY]]Y[^^\\^^\^^\\\]\^^^^^cc_ccc_YUUVVXVPUW[PV]YYSZYZYYSWWQTU[[[[]___[UUTTWW[[[[] NM:i:1 IH:i:1 HI:i:1
HWI-ST254_0000:5:1:2112:1951#0/1 0 chr5 134260898 83 100M * 0 0 NGGAGGTAGCGATGAGAGTAATAGATAGGGCTCAGGCGTTTGTGTATGATATGTTTGCGGTTTCGATGATGTGGTCTTTGGAGTAGAAACCTGTGCGGAA BWXWT[Q[WW^\^^^][[U[][[[]cccZc_____]U\\\SXTTXVVXYVVVVVV[[\]\^X^^BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB NM:i:4 IH:i:1 HI:i:1

The second argument is the folder of all output files.
