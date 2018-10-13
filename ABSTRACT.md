# RNA SV Simulator

In multicellular organisms, nearly every cell contains the same genome and thus the same genes. However, different cells show different gene expression underlying the different functions and tissues. Recent work in the cancer field, has shown that gene expression is more significant to disease outcomes than the genetic mutation [1].

Structural variants of genetic information are primary drivers of many diseases [2], making them a prime target for bioinformatic tools, but new tools require test datasets. Datasets are frequently large, ambiguous and subject to restrictive sharing permissions. However, sufficient numbers of publicly available genomes has allowed for the establishment of simulation tools of SVs.

Transcriptome sequencing is a much newer field, with a very low number of publicly available datasets — particularly in rare variants, significant to disease [3]. Here we present a tool that effectively replicates the expression profiles of rare and custom SVs across different tissue profiles to allow for the training and testing of bioinformatic tools.

[1]: https://doi.org/10.1038/nrg.2016.10
[2]: https://doi.org/10.1038/35015701
[3]: https://doi.org/10.3390%2Fijms18081652