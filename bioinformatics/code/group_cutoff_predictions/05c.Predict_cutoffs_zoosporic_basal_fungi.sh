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
if [ ! -d "dnabarcoder_zoosporic_basal_fungi" ]; then
  echo "Cloning dnabarcoder at: $(date)"
  git clone https://github.com/vuthuyduong/dnabarcoder.git dnabarcoder_zoosporic_basal_fungi
  
  # Remove files stored in the dnabarcoder data directory
  rm -r dnabarcoder_zoosporic_basal_fungi/data/*
  
  # Download the ITS1 extracted UNITE+INSD fasta file for fungi and the classification file formatted for dnabarcoder
  echo "Downloading UNITE..." $(date)
  wget -O dnabarcoder_zoosporic_basal_fungi/data/unite2024ITS1.classification "https://zenodo.org/records/13336328/files/unite2024ITS1.classification?download=1"
  wget -O dnabarcoder_zoosporic_basal_fungi/data/unite2024ITS1.fasta "https://zenodo.org/records/13336328/files/unite2024ITS1.fasta?download=1"
  
  # Truncate the dataset to max sequence length in ASVs (for analyses with forward reads and partial ITS reads retained)
  vsearch \
    --fastx_filter dnabarcoder_zoosporic_basal_fungi/data/unite2024ITS1.fasta \
    --fastq_trunclen_keep "$TRUNC_LEN" --fastq_minlen "$MIN_LEN" \
    --fastaout dnabarcoder_zoosporic_basal_fungi/data/unite2024ITS1.filtered.fasta
  
  # Create unique sequence dataset to speed up cutoff predictions
  echo "Selecting unique sequences..." $(date)
  python dnabarcoder_zoosporic_basal_fungi/aidscripts/selectsequences.py -i dnabarcoder_zoosporic_basal_fungi/data/unite2024ITS1.filtered.fasta -unique yes -c dnabarcoder_zoosporic_basal_fungi/data/unite2024ITS1.classification -o dnabarcoder_zoosporic_basal_fungi/data/unite2024ITS1.unique.fasta
  
else
  echo "dnabarcoder is already cloned..."
fi

# Work form the dnabarcoder directoy
cd dnabarcoder_zoosporic_basal_fungi

# Divide kingdom Fungi into phylogenetic clusters for predictions (e.g., "dikarya_fungi", "terrestrial_basal_fungi", "zoosporic_basal_fungi")
Rscript --vanilla "../prepare_dnabarcoder.R" >> "../slurm/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.out"

################################################################################
# Predict species cutoffs
################################################################################

echo 'Starting species cutoff predictions at:' $(date)

# Select sequences for species
echo 'Selecting sequences for species at:' $(date)
python aidscripts/selectsequences.py -i data/zoosporic_basal_fungi.fasta -c data/zoosporic_basal_fungi.classification -rank species -o zoosporic_basal_fungi.species.fasta 

# Predict local similarity cutoffs at the species level
echo 'Predicting local similarity cutoffs for species at:' $(date)
python dnabarcoder.py predict -i zoosporic_basal_fungi.species.fasta -c zoosporic_basal_fungi.species.classification -st 0.8 -et 1 -s 0.001 -rank species -prefix zoosporic_basal_fungi -higherrank genus,family,order,class,phylum,kingdom -minseqno 30 -mingroupno 10 -maxproportion 0.5 -removecomplexes yes -ml 50

# Remove intermediate files
rm zoosporic_basal_fungi.species.*

echo 'Finised species cutoff predictions at:' $(date)

################################################################################
# Predict genus cutoffs
################################################################################

echo 'Starting genus cutoff predictions at:' $(date)

# Select sequences for genus
echo 'Selecting sequences for genus at:' $(date)
python aidscripts/selectsequences.py -i data/zoosporic_basal_fungi.fasta -c data/zoosporic_basal_fungi.classification -rank genus -o zoosporic_basal_fungi.genus.fasta  

# Predict local similarity cutoffs at the genus level
echo 'Predicting local similarity cutoffs for genus at:' $(date)
python dnabarcoder.py predict -i zoosporic_basal_fungi.genus.fasta -c zoosporic_basal_fungi.genus.classification -st 0.6 -et 1 -s 0.001 -prefix zoosporic_basal_fungi -rank genus -higherrank family,order,class,phylum,kingdom -minseqno 30 -mingroupno 10 -maxproportion 0.5 -removecomplexes yes -ml 50

# Remove intermediate files
rm zoosporic_basal_fungi.genus.*

echo 'Finised genus cutoff predictions at:' $(date)

################################################################################
# Predict family cutoffs
################################################################################

echo 'Starting family cutoff predictions at:' $(date)

# Select sequences for family
echo 'Selecting sequences for family at:' $(date)
python aidscripts/selectsequences.py -i data/zoosporic_basal_fungi.fasta -c data/zoosporic_basal_fungi.classification -rank family -o zoosporic_basal_fungi.family.fasta  

# Predict local similarity cutoffs at the family level
echo 'Predicting local similarity cutoffs for family at:' $(date)
python dnabarcoder.py predict -i zoosporic_basal_fungi.family.fasta -c zoosporic_basal_fungi.family.classification -st 0.5 -et 1 -s 0.001 -prefix zoosporic_basal_fungi -rank family -higherrank order,class,phylum,kingdom -minseqno 30 -mingroupno 10 -maxproportion 0.5 -removecomplexes yes -ml 50

# Remove intermediate files
rm zoosporic_basal_fungi.family.*

echo 'Finised family cutoff predictions at:' $(date)

################################################################################
# Predict order cutoffs
################################################################################

echo 'Starting order cutoff predictions at:' $(date)

# Select sequences for order
echo 'Selecting sequences for order at:' $(date)
python aidscripts/selectsequences.py -i data/zoosporic_basal_fungi.fasta -c data/zoosporic_basal_fungi.classification -rank order -o zoosporic_basal_fungi.order.fasta  

# Predict local similarity cutoffs at the order level
echo 'Predicting local similarity cutoffs for order at:' $(date)
python dnabarcoder.py predict -i zoosporic_basal_fungi.order.fasta -c zoosporic_basal_fungi.order.classification -st 0.5 -et 1 -s 0.001 -prefix zoosporic_basal_fungi -rank order -higherrank class,phylum,kingdom -minseqno 30 -mingroupno 10 -maxproportion 0.5 -removecomplexes yes -ml 50

# Remove intermediate files
rm zoosporic_basal_fungi.order.*

echo 'Finised order cutoff predictions at:' $(date)

################################################################################
# Predict class cutoffs
################################################################################

echo 'Starting class cutoff predictions at:' $(date)

# Select sequences for class
echo 'Selecting sequences for class at:' $(date)
python aidscripts/selectsequences.py -i data/zoosporic_basal_fungi.fasta -c data/zoosporic_basal_fungi.classification -rank class -o zoosporic_basal_fungi.class.fasta  

# Predict local similarity cutoffs at the class level
echo 'Predicting local similarity cutoffs for class at:' $(date)
python dnabarcoder.py predict -i zoosporic_basal_fungi.class.fasta -c zoosporic_basal_fungi.class.classification -st 0.5 -et 1 -s 0.001 -prefix zoosporic_basal_fungi -rank class -higherrank phylum,kingdom -minseqno 30 -mingroupno 10 -maxproportion 0.5 -removecomplexes yes -ml 50

# Remove intermediate files
rm zoosporic_basal_fungi.class.*

echo 'Finised class cutoff predictions at:' $(date)

# Deactivate the conda environment
conda deactivate