# The Australian Microbiome soil ITS1 dataset

This repository contains the code associated to the paper: The Australian Microbiome initiative for accurate detection and identification of fungi.

Authors
Luke Florence<sup>1<sup>, Sean Tomlinson<sup>2<sup>, Marc Freestone<sup>3<sup>, John Morgan<sup>1<sup>, Jen Wood<sup>4<sup>, Camille Truong<sup>3<sup>

Affiliations
1. Department of Environment and Genetics, La Trobe University, Bundoora, VIC, 3083, Australia.
2. Biodiversity and Conservation Science, Department of Biodiversity, Conservation and Attractions, Kensington, WA, 6151, Australia.
3. Royal Botanic Gardens Victoria, Melbourne, VIC, 3004, Australia.
4. Department of Department of Microbiology, Anatomy, Physiology and Pharmacology, La Trobe University, Bundoora, VIC, 3083, Australia.

Corresponding author: Luke Florence (L.Florence@latrobe.edu.au)

## Overview

DNA metabarcoding has played a pivotal role in advancing our understanding of soil-inhabiting fungi worldwide. The [Australian Microbiome Initiative](https://www.australianmicrobiome.com/) has produced an extensive soil fungal dataset that covers more than 2000 plots in a variety of bioregions and ecosystem types. The integration of Australian microbiome sequence data into platforms like the Global Biodiversity Information Facility holds immense promise to provide valuable insights into fungal biodiversity and ecology in the Southern Hemisphere. However, recent studies raised concerns about the uneven quality of the Australian Microbiome fungal dataset, highlighting overinflated diversity and erroneous species occurrences (false positives). To address these concerns, we aim to generate a robust dataset that can be leveraged by end-users interested in the conservation and ecology of fungi in soils. We reanalysed the Australian Microbiome fungal dataset using a conservative approach following best practices in metabarcoding and carefully ground-truthed our methodology by comparing our contemporary dataset to studies conducted on historical versions the dataset. This resource will significantly benefit future research on soil fungi in Australia and beyond.

The scripts used to process the raw data and generate the final dataset are available in the `code` directory. The associated raw data files can be obtained from the [Bioplatforms Australia data portal](https://data.bioplatforms.com/organization/australian-microbiome) by selecting ITS data from terrestrial soil environments. The final dataset is available on the dedicated [figshare data repository]().

## Repository contents

**Note** that shell scripts are run from the source file but R scripts are run from an R project directory (i.e. the project base directory).

**01.Extract_ITS.sh** performs five key functions:
1. Quality truncate reads with Trimmomatic
2. Extract the ITS region with ITSxpress
3. Quality filter reads with VSEARCH
4. Quality check reads quality with FastQC and MultiQC
5. Track reads across the pipeline using a library-wise approach

**02.Denoise.R** can be executed from a HPC cluster by running the **02.Denoise.sh** and performs the following tasks:
1. Denoise reads with DADA2
2. Merge sequence tables from multiple sequencing runs
3. Convert the merged sequence table to a FASTA sequence file formatted for chimera detection in VSEARCH

**03.Prepare_UNITE.sh** downloads the UNITE+INSD reference dataset used in this study and extracts the ITS1 region for use in chimera detection and taxonomic assignment.

**04.Chimera_detection.sh** performs two main tasks:
1. De novo chimera detection with VSEARCH
2. Reference-based chimera detection with VSEARCH

**05.Cluster.sh** clusters the sequences into operational taxonomic units (OTUs) at 97% similarity using an abundance-based centroid approach in VSEARCH.

**06.BLAST.sh** assigns top ten BLAST hits to OTUs using the UNITE+INSD reference dataset.

**07.Abundance_filter_OTUs**

**08.Quality_filter_BLAST.R**

**09.Annotate_OTUs.R**

**10.Collect_metadata.R**

**Other scripts:**
The **map.pl** script is used by the **04.Chimera_detection.sh** script to map dereplicated reads back to their original samples. The **functions.R** script contains custom functions used in multiple R scripts within this repository.

## Dependencies

This project was developed using a Unix computing environment and the following software:

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
[rnaturalearth]()
[rnaturalearthdata]()
[patchwork]()
