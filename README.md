## gru 
Gru is a small pipeline for evaluation of (meanwhile obolete) 2D Nanopore reads.
it requires bwa, poretools, samtools and bamstats. 

This pipeline can be extended in the future by an assembler such as minasm or spades (to colse Illumina gaps).

Unfrotunately the 2D reads a not more supported by Nanopore.

"Brown said 1D2 read sequencing will be released to customers that are licensed as developers immediately, and more broadly around the time of its UK user meeting in May. He said the accuracy of 1D2 reads is significantly better than that of the old 2D reads, which will be discontinued in early May."
(see [this article](https://www.genomeweb.com/sequencing/oxford-nanopore-launches-gridion-x5-nanopore-sequencer-details-product-improvements) )

Usage '''/gru.py job.yml'''
