# The Australian Microbiome dataset for soil fungi

This repository contains the code associated to the paper:

**The Australian Microbiome dataset for accurate detection and ecological modelling of fungi**

Authors:

Luke Florence<sup>1</sup>, Sean Tomlinson<sup>2,3</sup>, Marc Freestone<sup>4,5</sup>, John W. Morgan<sup>1</sup>, Jen L. Wood<sup>6</sup>, Camille Truong<sup>4</sup>

Affiliations:
1. Department of Environment and Genetics, La Trobe University, Bundoora, VIC, 3083, Australia.
2. Biodiversity and Conservation Science, Department of Biodiversity, Conservation and Attractions, Kensington, WA, 6151, Australia.
3. School of Biological Sciences, University of Adelaide, Adelaide, SA, 5000, Australia.
4. Royal Botanic Gardens Victoria, Melbourne, VIC, 3004, Australia.
5. Ecology and Evolution, Research School of Biology, Australian National University, Canberra, ACT, 2600, Australia.
6. Department of Department of Microbiology, Anatomy, Physiology and Pharmacology, La Trobe University, Bundoora, VIC, 3083, Australia.

Corresponding author: Luke Florence (L.Florence@latrobe.edu.au)

## Overview

DNA metabarcoding has played a pivotal role in advancing our understanding of soil-inhabiting fungi worldwide. The [Australian Microbiome Initiative](https://www.australianmicrobiome.com/) has produced an extensive soil fungal dataset, sampling more than 2000 plots across a breadth of ecosystems in Australia and Antarctica. Recent studies have highlighted overinflated diversity and false positives in the Australian Microbiome fungal dataset, leading to erroneous species occurrences. To address this, we reanalysed the dataset using a conservative approach, following best practices in metabarcoding. We carefully ground-truthed our methodology by comparing our contemporary dataset to studies conducted on historical versions of the dataset. Our approach generated a robust dataset that can be leveraged by end-users interested in the biodiversity, distribution, and conservation of fungi in soils. We paired this dataset with over 100 predictor variables to fast-track data exploration.

The scripts used to process the raw data and generate the primary OTU and taxonomy files are available in the `code/bioinformatics` directory. Code associated with downstream data curation and technical validation are available in the `code/data_curation` directory. The associated raw data files can be obtained from the [Bioplatforms Australia data portal](https://data.bioplatforms.com/organization/australian-microbiome) by selecting ITS data from terrestrial soil environments. The primary data products from this analysis are available on the dedicated [figshare data repository](https://doi.org/10.26181/25706673).

## Repository contents

**Note** that shell scripts are run from the source file but R scripts are run from an R project directory (i.e. the project base directory).

### Bionformatics

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
1. Sample-wise abundance filter of OTUs in each sample based on relative abundance <0.1% of the total sequence count per sample
2. Library-wise abundance filter of each OTU with a relative sequence abundance <0.5% of the total OTU within a given library
3. Positive control filter of remove positive control OTUs with a relative sequence abundance <3% of the total positive control OTU within a given library
4. Dereplicate samples and remove samples with sequencing deopth <5,000 reads

**08.Quality_filter_BLAST.R**: Remove poor quality BLAST hits based on similarity and couverage thresholds of 70% with some exceptions.

**09.Annotate_OTUs.R**: Annotate OTUs based on taxa-specific thresholds and 66.67% consensus of across the top three hits. 

The **map.pl** script is used by the **04.Chimera_detection.sh** script to map dereplicated reads back to their original samples. The **functions.R** script contains custom R functions used in the bioinformatics analyais.

### Data curation

**1.Collect_metadata.R**: Calculate basic diversity metrics, estimate mould contaminatuion, and gather predictor variables.

**02.Check_mould_contamination.R**: Evaluate the impact of mould contamination on sample richness.

**03.Normalise_OTUs.R**:
1. Evaluate the relationship between seuqnecing depth and OTU richness.
2. Normalise OTU counts using the scaling with ranked subsampling method.

**04.ECM_antarctica.R**: Evaluate the occurrence of ectomycorrhizal OTUs in Antarctica.

**05.Map_cortinarius**: Map the distribution of Cortinarius species in Australia.

**06.Map_amanita**: Map the distribution of Amanita species in Australia.

**07.Compare_filtering_thresholds.R**: Compare the impact of different filtering thresholds on OTU richness, abundance and prevalence.

The **functions.R** script contains custom R functions used in the data curation process.

## Dependencies

This project was developed using a Unix computing environment and the following software:

* [R](https://cran.r-project.org/) version 4.3.1 -- "Beagle Scouts"
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) version 0.36
* [ITSxpress](https://github.com/USDA-ARS-GBRU/itsxpress) version 2.0.0
* [VSEARCH](https://github.com/torognes/vsearch) version 2.22.1
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) version 0.12.1
* [MultiQC](https://multiqc.info/) version 1.15
* [DADA2](https://benjjneb.github.io/dada2/) version 1.26
* [seqinr](https://github.com/lbbe-software/seqinr) version 4.2-30
* [here](https://here.r-lib.org/) version 1.0.1
* [tidyverse](https://www.tidyverse.org/) version 2.0.0
* [ITSx](https://microbiology.se/software/itsx/) version 1.1.3
* [HMMER](http://hmmer.org/) version 3.1b2 
* [Perl](https://www.perl.org/get.html) version 5.32.1
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) version 2.14.1
* [data.table](https://github.com/Rdatatable/data.table) version 1.14.10
* [ggpubr](https://github.com/kassambara/ggpubr) version 0.4.0
* [terra](https://rspatial.org/pkg/) version 1.7-71
* [tidyterra](https://dieghernan.github.io/tidyterra/) version 0.5.2
* [sf](https://r-spatial.github.io/sf/) version 1.0-15
* [SRS]() version 0.2.3
* [scales]() version 1.3.0
* [vegan](https://cran.r-project.org/web/packages/vegan/index.html) version 2.4-6
* [rnaturalearth](https://cran.r-project.org/web/packages/rnaturalearth/vignettes/rnaturalearth.html) version 1.0.1
* [rnaturalearthdata](https://docs.ropensci.org/rnaturalearthdata/) version 1.0.0
* [patchwork](https://patchwork.data-imaginist.com/) version 1.2.0
