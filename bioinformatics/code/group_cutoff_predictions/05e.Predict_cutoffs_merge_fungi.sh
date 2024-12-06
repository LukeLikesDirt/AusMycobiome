#!/bin/bash

# Activate the conda environment
source ~/.bashrc
conda activate dynamic_clustering

################################################################################
# Prepare dnabarcoder for predictions
################################################################################

# Constants
export MIN_LEN=50
export TRUNC_LEN=500

# Clone dnabarcoder
echo "Cloning dnabarcoder..." $(date)
if [ ! -d "dnabarcoder" ]; then
  echo "Cloning dnabarcoder at: $(date)"
  git clone https://github.com/vuthuyduong/dnabarcoder.git dnabarcoder
  
  # Remove files stored in the dnabarcoder data directory
  rm -r dnabarcoder/data/*
  
  # Download the ITS1 extracted UNITE+INSD fasta file for fungi and the classification file formatted for dnabarcoder
  echo "Downloading UNITE..." $(date)
  wget -O dnabarcoder/data/unite2024ITS1.classification "https://zenodo.org/records/13336328/files/unite2024ITS1.classification?download=1"
  wget -O dnabarcoder/data/unite2024ITS1.fasta "https://zenodo.org/records/13336328/files/unite2024ITS1.fasta?download=1"
  
  # Truncate the dataset to max sequence length in ASVs (for analyses with forward reads and partial ITS reads retained)
  vsearch \
    --fastx_filter dnabarcoder/data/unite2024ITS1.fasta \
    --fastq_trunclen_keep "$TRUNC_LEN" --fastq_minlen "$MIN_LEN" \
    --fastaout dnabarcoder/data/unite2024ITS1.filtered.fasta
  
  # Create unique sequence dataset to speed up cutoff predictions
  echo "Selecting unique sequences..." $(date)
  python dnabarcoder/aidscripts/selectsequences.py -i dnabarcoder/data/unite2024ITS1.filtered.fasta -unique yes -c dnabarcoder/data/unite2024ITS1.classification -o dnabarcoder/data/unite2024ITS1.unique.fasta
  
else
  echo "dnabarcoder is already cloned..."
fi

# Work form the dnabarcoder directoy
cd dnabarcoder

# Divide kingdom Fungi into phylogenetic clusters for predictions (e.g., "dikarya_fungi", "terrestrial_basal_fungi", "zoosporic_basal_fungi")
Rscript --vanilla "../prepare_dnabarcoder.R" >> "../slurm/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.out"

################################################################################
# Predict global cutoffs for fungal phyla
################################################################################

# Predict global similarity cutoffs at the phylum level
echo 'Predicting global similarity cutoffs for phylum at:' $(date)
python dnabarcoder.py predict -i data/unite2024ITS1.unique.fasta -c data/unite2024ITS1.unique.classification -st 0.5 -et 1 -s 0.001 -rank phylum -prefix fungi -removecomplexes yes -ml 50
echo 'Finised phylum cutoff predictions at:' $(date)

################################################################################
# Compute best cutoffs across groups and merge files
################################################################################

# NOTE: I do some manual curation, changing the name of "All" in global
# predicitions to group names prior to merging. I also do post merging curation
# to remove ambiguous values. For example, group number of 10, minimum sequence
# number of 30 and max portion of 50% was mostly stable with few ambiguous cutoff
# values. But some obviously too low cutoff values with high confidance values
# still remain for a few groups that are boarder line to the dnabarcoder options
# values used.

# Compute best cutoffs for locall cutoffs:
python dnabarcoder.py best -i ../dnabarcoder_dikarya_fungi/dnabarcoder/dikarya.cutoffs.json -c ./data/dikarya_fungi.classification
python dnabarcoder.py best -i ../dnabarcoder_terrestrial_basal_fungi/dnabarcoder/terrestrial_basal_fungi.cutoffs.json -c ./data/terrestrial_basal_fungi.classification
python dnabarcoder.py best -i ../dnabarcoder_zoosporic_basal_fungi/dnabarcoder/zoosporic_basal_fungi.cutoffs.json -c ./data/zoosporic_basal_fungi.classification

# Merge global cutoffs
python dnabarcoder.py merge -i ./dnabarcoder/global_cutoffs.json,../dnabarcoder_global_fungi/dnabarcoder/dikarya_fungi.cutoffs.json,../dnabarcoder_global_fungi/dnabarcoder/terrestrial_basal_fungi.cutoffs.json,../dnabarcoder_global_fungi/dnabarcoder/zoosporic_basal_fungi.cutoffs.json -o ./dnabarcoder/global_cutoffs.json

# Merge all cutoffs
python dnabarcoder.py merge -i ./dnabarcoder/global_cutoffs.json,./dnabarcoder/dikarya.cutoffs.best.json,./dnabarcoder/terrestrial_basal_fungi.cutoffs.best.json,./zoosporic_basal_fungi.cutoffs.best.json -o ./data/aus_mic_cutoffs.json

conda deactivate