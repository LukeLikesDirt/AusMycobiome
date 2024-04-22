# The Australian Microbiome soil ITS1 dataset

This repository contains the code associated to the paper: (TBA).

## Overview

This project is dedicated to providing a publicly accessible, well-curated version of the [Australian Microbiome Initiative](https://www.australianmicrobiome.com/) (AusMicrobiome) soil ITS dataset. The ITS sequencing procedure employed by AusMicrobiome captures the ITS1 and ITS2 regions, utilising the Illumina MiSeq short-read platform. Due to the length of the ITS regions, forward and reverse reads cannot be merged, which, along with other procedural constraints, contributes to an overestimation of diversity in the current AusMicrobiome ITS dataset. To that end, we have re-processed the entire AusMicrobiome soil ITS1 dataset using a relatively conservative approach. We aim to provide a readily available resource to assist researchers in assessing the diversity of Australian soil fungal communities. This dataset has been paired with a run-of-the-mill diversity analysis, as well as maps highlighting areas of uncertainty regarding diversity measures, to identify areas in need of sampling.

## Contents

#### Bioinformatics pipeline

#### Statistical analysys

## Data availability

## Dependencies

This project was developed using a shell computing environment and the following software:

[R](https://cran.r-project.org/) version 4.3.1 -- "Beagle Scouts"
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) version 0.36
[ITSxpress](https://github.com/USDA-ARS-GBRU/itsxpress) v2.0.0
[VSEARCH](https://github.com/torognes/vsearch) v2.22.1
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) v0.12.1
[MultiQC](https://multiqc.info/) version 1.15
[DADA2](https://benjjneb.github.io/dada2/) version 3.18
[seqinr](https://github.com/lbbe-software/seqinr) version 4.2-30
[here](https://here.r-lib.org/) version 1.0.1
[tidyverse](https://www.tidyverse.org/) version 2.0.0
[ITSx](https://microbiology.se/software/itsx/) version 1.1.3
[HMMER](http://hmmer.org/) version 3.1b2 
[Perl](https://www.perl.org/get.html) version 5.32.1
[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) version 2.14.1
[data.table](https://github.com/Rdatatable/data.table) version 1.14.10
[ggpubr](https://github.com/kassambara/ggpubr) version 0.4.0
[rnaturalearth]
[rnaturalearthdata]
[patchwork]

## Contact
