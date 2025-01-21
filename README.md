# The Australian Microbiome dataset for soil fungi

This repository contains the code associated to the paper:

**A curated dataset for accurate detection of soil fungi to advance ecological research and conservation in Australia and Antarctica**

Authors:

Luke Florence<sup>1</sup>, Sean Tomlinson<sup>2,3</sup>, Marc Freestone<sup>4</sup>, John W. Morgan<sup>1</sup>, Jen L. Wood<sup>5</sup>, Camille Truong<sup>6</sup>

Affiliations:
1. Department of Environment and Genetics, La Trobe University, Bundoora, VIC, 3083, Australia.
2. Biodiversity and Conservation Science, Department of Biodiversity, Conservation and Attractions, Kensington, WA, 6151, Australia.
3. School of Biological Sciences, University of Adelaide, Adelaide, SA, 5000, Australia.
4. The Biodiversity Consultancy, Cambridge, CB2 1SJ, United Kingdom.
5. Department of Department of Microbiology, Anatomy, Physiology and Pharmacology, La Trobe University, Bundoora, VIC, 3083, Australia.
6. Royal Botanic Gardens Victoria, Melbourne, VIC, 3004, Australia.

Corresponding author: Luke Florence (L.Florence@latrobe.edu.au)

## Overview

DNA metabarcoding has played a pivotal role in advancing our understanding of the diversity and function of soil-inhabiting fungi. The [Australian Microbiome Initiative](https://www.australianmicrobiome.com/) has produced an extensive soil fungal metabarcoding dataset of more than 2000 plots across a breadth of ecosystems in Australia and Antarctica. Sequence data requires rigorous approaches for the integration of species occurrences into biodiversity platforms, addressing biases due to false positives and overinflated diversity estimates, among others. To address such biases, we conducted a rigorous analysis of the fungal dataset following best practices in fungal metabarcoding and paired this dataset with over 100 predictor variables to fast-track data exploration. We carefully validated our methodology based on studies conducted on historical versions of the dataset. Our approach generated robust information on Australian soil fungi that can be leveraged by end-users interested in biodiversity, biogeography, and conservation.

The primary data products from this analysis are available on the associated [figshare repository](https://doi.org/10.6084/m9.figshare.27938037).

## Repository contents

### Bionformatics

**Notes:** 
* The reproducible code for the bioinformatics pipeline are available in the `bioinformatics` directory and are run from source directory.
* The associated raw data files can be obtained from the [Bioplatforms Australia data portal](https://data.bioplatforms.com/organization/australian-microbiome) by using the search terms “sample_type:Soil & amplicon:ITS & depth_lower:0.1”.
* The dependencies required to reproduce this research can be installed by via mamba and the `env/dynamic_cluster.yml` file.
* This pipeline has been designed to be reproducible. However, this is a first draft of the taxonomically informed dynamic clustering and some troubleshooting will be required for replicability (i.e. to run the pipeline on new datasets).
* After organising the raw data, the scripts with numeric preifixes should be run in numeric order to reproduce the results. Scripts without numeric prefixes are auxiliary code.

**01.Extract_ITS.sh** performs five functions:
1. Quality truncate reads with Trimmomatic
2. Extract the ITS region with ITSxpress
3. Quality filter reads with VSEARCH
4. Check the quality of processed reads with FastQC and MultiQC
5. Track reads across the pipeline (library-wise approach)

**02.Denoise.R** performs the following tasks:
1. Denoise reads with DADA2
2. Merge sequence tables from multiple sequencing runs
3. Convert the merged sequence table to a FASTA sequence file formatted for chimera detection in VSEARCH

**03.Chimera_detection.sh** performs two main tasks:
1. De novo chimera detection with VSEARCH
2. Reference-based chimera detection with VSEARCH
*Note:* There are two distinct output from this and subsequent steps: (1) an `ASVs` and `OTUs` output. The `ASVs` are have ASVs that are clustered after taxonomic assignment using a taxonomically informed dynamic clustering approach and the `OTUs are clustered in this step at 97% similarity using an abundance-based centroid approach in VSEARCH. The main data output from this project uses the dynamic clusters. (ie.e the ASV files from this step) and the OTU files are intended for comparative analyses between the conventional 97% OTU approach and the dynamic clustering approach we chose to use.

**04.Abundance_filter_OTUs**
1. Sample-wise abundance filter of OTUs in each sample based on relative abundance <0.1% of the total sequence count per sample
2. Library-wise abundance filter of each OTU with a relative sequence abundance <0.5% of the total OTU within a given library
3. Positive control filter of remove positive control OTUs with a relative sequence abundance <3% of the total positive control OTU within a given library
4. Dereplicate samples and remove samples with sequencing deopth <5,000 reads

**05.Predict_cutoffs.sh** downloads the UNITE+INSD reference dataset used in this study and extracts the ITS1 region for use in chimera detection and taxonomic assignment.

**06.BLAST.sh** assigns top five BLAST hits to ASVs and OTUs using the UNITE+INSD reference dataset.

**07.Classify_OTUs.sh** Filter BLAST top five hits based on taxon-specific thresholds, coverage thresholds for genus (90%) and species (95%), and then affiliate taxonomy to ASVs and OTUs to the best hit at each rank using 66.67% consensus threshold across the all remaining hits. 

**08.Dynamic_clustering.sh** Cluster OTUs using a taxonomically informed dynamic clustering approach

### Tchnical validation

**Note:**
* Code associated with the technical validation are available in the `technical_validation` directory and are executable from within the R project.

**01.Assess_length_bias.R**: Assess sequences length distribution of the ITS region in the Australian Microbiome dataset against the UNITE+INSD reference dataset.

**02.Impact_of_ITS_extraction.R**: Assess the impact of ITS extraction on biases against taxon with long ITS regions. We processed 300 bp sequences targeting the fungal ITS1 region, removing conserved SSU and 5.8S rRNA sequences to enhance OTU clustering and taxonomic accuracy. Using ITSxpress for ITS1 extraction, we noted that the inability to merge paired-end reads in the Australian Microbiome fungal dataset biases against taxa with long ITS1 regions (>230 bp), leading to potential false negatives. This script identifies taxa with long ITS1 regions that are underrepresented in the Australian Microbiome dataset.

**03.Quantify_ECM_in_antarctica.R**: Evaluate the occurrence of ectomycorrhizal OTUs in Antarctica.

**04.Map_cortinarius.R**: Map the distribution of Cortinarius species in Australia.

**05.Map_amanita.R**: Map the distribution of Amanita species in Australia.

**06.Example_niche_analysis.R**: Test differences in the soil total nitrogen and phosphorus niches of two hypogeous (belowground) ectomycorrhizal genera from the order Pezizales: Ruhlandiella and Sphaerosoma.

**07.Assess_dataset_diversity.R**:

**Technical validation dependencies**

The technical validation was conducted using [R version 4.3.3 -- "Angel Food Cake"](https://cran.r-project.org/) and the following packages:

* [Biostrings](https://github.com/lbbe-software/seqinr) version 2.70.3
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) version 2.16.0
* [data.table](https://github.com/Rdatatable/data.table) version 1.15.4
* [emmeans](https://rvlenth.github.io/emmeans/) version 1.10.1
* [ggbeeswarm](https://cran.r-project.org/web/packages/ggbeeswarm/index.html) version 0.7.2
* [ggpubr](https://github.com/kassambara/ggpubr) version 0.6.0
* [parameters](https://easystats.github.io/parameters/) version 0.21.7
* [patchwork](https://patchwork.data-imaginist.com/) version 1.2.0
* [performance](https://easystats.github.io/performance/) version 0.11.0
* [rnaturalearth](https://cran.r-project.org/web/packages/rnaturalearth/vignettes/rnaturalearth.html) version 1.0.1
* [rnaturalearthdata](https://docs.ropensci.org/rnaturalearthdata/) version 1.0.0
* [SRS](https://github.com/vitorheidrich/SRS) version 0.2.3
* [sf](https://r-spatial.github.io/sf/) version 1.0-16
* [terra](https://rspatial.org/pkg/) version 1.7-71
* [tidyterra](https://dieghernan.github.io/tidyterra/) version 0.6.0
* [tidyverse](https://www.tidyverse.org/) version 2.0.0
* [vegan](https://cran.r-project.org/web/packages/vegan/index.html) version 2.6-8

