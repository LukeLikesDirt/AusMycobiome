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
if [ ! -d "dnabarcoder_global_fungi" ]; then
  echo "Cloning dnabarcoder at: $(date)"
  git clone https://github.com/vuthuyduong/dnabarcoder.git dnabarcoder_global_fungi
  
  # Remove files stored in the dnabarcoder data directory
  rm -r dnabarcoder_global_fungi/data/*
  
  # Download the ITS1 extracted UNITE+INSD fasta file for fungi and the classification file formatted for dnabarcoder
  echo "Downloading UNITE..." $(date)
  wget -O dnabarcoder_global_fungi/data/unite2024ITS1.classification "https://zenodo.org/records/13336328/files/unite2024ITS1.classification?download=1"
  wget -O dnabarcoder_global_fungi/data/unite2024ITS1.fasta "https://zenodo.org/records/13336328/files/unite2024ITS1.fasta?download=1"
  
  # Truncate the dataset to max sequence length in ASVs (for analyses with forward reads and partial ITS reads retained)
  vsearch \
    --fastx_filter dnabarcoder_global_fungi/data/unite2024ITS1.fasta \
    --fastq_trunclen_keep "$TRUNC_LEN" --fastq_minlen "$MIN_LEN" \
    --fastaout dnabarcoder_global_fungi/data/unite2024ITS1.filtered.fasta
  
  # Create unique sequence dataset to speed up cutoff predictions
  echo "Selecting unique sequences..." $(date)
  python dnabarcoder_global_fungi/aidscripts/selectsequences.py -i dnabarcoder_global_fungi/data/unite2024ITS1.filtered.fasta -unique yes -c dnabarcoder_global_fungi/data/unite2024ITS1.classification -o dnabarcoder_global_fungi/data/unite2024ITS1.unique.fasta
  
else
  echo "dnabarcoder is already cloned..."
fi

# Work form the dnabarcoder directoy
cd dnabarcoder_global_fungi

# Divide kingdom Fungi into phylogenetic clusters for predictions (e.g., "dikarya_fungi", "terrestrial_basal_fungi", "zoosporic_basal_fungi")
Rscript --vanilla "../prepare_dnabarcoder.R" >> "../slurm/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.out"

################################################################################
# Predict species cutoffs
################################################################################

echo 'Starting species cutoff predictions at:' $(date)

# Predict global similarity cutoffs at the species level
echo 'Predicting global similarity cutoffs for species at:' $(date)
python dnabarcoder.py predict -i data/unite2024ITS1.unique.fasta -c data/unite2024ITS1.unique.classification -st 0.8 -et 1 -s 0.001 -rank species -prefix fungi -removecomplexes yes -ml 50
python dnabarcoder.py predict -i data/dikarya_fungi.fasta -c data/dikarya_fungi.classification -st 0.8 -et 1 -s 0.001 -rank species -prefix dikarya_fungi -removecomplexes yes -ml 50
python dnabarcoder.py predict -i data/terrestrial_basal_fungi.fasta -c data/terrestrial_basal_fungi.classification -st 0.8 -et 1 -s 0.001 -rank species -prefix terrestrial_basal_fungi -removecomplexes yes -ml 50
python dnabarcoder.py predict -i data/zoosporic_basal_fungi.fasta -c data/zoosporic_basal_fungi.classification -st 0.8 -et 1 -s 0.001 -rank species -prefix zoosporic_basal_fungi -removecomplexes yes -ml 50

echo 'Finised species cutoff predictions at:' $(date)

################################################################################
# Predict genus cutoffs
################################################################################

# Predict global similarity cutoffs at the genus level
echo 'Predicting global similarity cutoffs for genus at:' $(date)
python dnabarcoder.py predict -i data/unite2024ITS1.unique.fasta -c data/unite2024ITS1.unique.classification -st 0.5 -et 1 -s 0.001 -rank genus -prefix fungi -removecomplexes yes -ml 50
python dnabarcoder.py predict -i data/dikarya_fungi.fasta -c data/dikarya_fungi.classification -st 0.5 -et 1 -s 0.001 -rank genus -prefix dikarya_fungi -removecomplexes yes -ml 50
python dnabarcoder.py predict -i data/terrestrial_basal_fungi.fasta -c data/terrestrial_basal_fungi.classification -st 0.5 -et 1 -s 0.001 -rank genus -prefix terrestrial_basal_fungi -removecomplexes yes -ml 50
python dnabarcoder.py predict -i data/zoosporic_basal_fungi.fasta -c data/zoosporic_basal_fungi.classification -st 0.5 -et 1 -s 0.001 -rank genus -prefix zoosporic_basal_fungi -removecomplexes yes -ml 50
echo 'Finised genus cutoff predictions at:' $(date)

################################################################################
# Predict family cutoffs
################################################################################

# Predict global similarity cutoffs at the family level
echo 'Predicting global similarity cutoffs for family at:' $(date)
python dnabarcoder.py predict -i data/unite2024ITS1.unique.fasta -c data/unite2024ITS1.unique.classification -st 0.5 -et 1 -s 0.001 -rank family -prefix fungi -removecomplexes yes -ml 50
python dnabarcoder.py predict -i data/dikarya_fungi.fasta -c data/dikarya_fungi.classification -st 0.5 -et 1 -s 0.001 -rank family -prefix dikarya_fungi -removecomplexes yes -ml 50
python dnabarcoder.py predict -i data/terrestrial_basal_fungi.fasta -c data/terrestrial_basal_fungi.classification -st 0.5 -et 1 -s 0.001 -rank family -prefix terrestrial_basal_fungi -removecomplexes yes -ml 50
python dnabarcoder.py predict -i data/zoosporic_basal_fungi.fasta -c data/zoosporic_basal_fungi.classification -st 0.5 -et 1 -s 0.001 -rank family -prefix zoosporic_basal_fungi -removecomplexes yes -ml 50
echo 'Finised family cutoff predictions at:' $(date)

################################################################################
# Predict order cutoffs
################################################################################

# Predict global similarity cutoffs at the order level
echo 'Predicting global similarity cutoffs for order at:' $(date)
python dnabarcoder.py predict -i data/unite2024ITS1.unique.fasta -c data/unite2024ITS1.unique.classification -st 0.5 -et 1 -s 0.001 -rank order -prefix fungi -removecomplexes yes -ml 50
python dnabarcoder.py predict -i data/dikarya_fungi.fasta -c data/dikarya_fungi.classification -st 0.5 -et 1 -s 0.001 -rank order -prefix dikarya_fungi -removecomplexes yes -ml 50
python dnabarcoder.py predict -i data/terrestrial_basal_fungi.fasta -c data/terrestrial_basal_fungi.classification -st 0.5 -et 1 -s 0.001 -rank order -prefix terrestrial_basal_fungi -removecomplexes yes -ml 50
python dnabarcoder.py predict -i data/zoosporic_basal_fungi.fasta -c data/zoosporic_basal_fungi.classification -st 0.5 -et 1 -s 0.001 -rank order -prefix zoosporic_basal_fungi -removecomplexes yes -ml 50
echo 'Finised order cutoff predictions at:' $(date)

################################################################################
# Predict class cutoffs
################################################################################

# Predict global similarity cutoffs at the class level
echo 'Predicting global similarity cutoffs for class at:' $(date)
python dnabarcoder.py predict -i data/unite2024ITS1.unique.fasta -c data/unite2024ITS1.unique.classification -st 0.5 -et 1 -s 0.001 -rank class -prefix fungi --removecomplexes yes -ml 50
python dnabarcoder.py predict -i data/dikarya_fungi.fasta -c data/dikarya_fungi.classification -st 0.5 -et 1 -s 0.001 -rank class -prefix dikarya_fungi --removecomplexes yes -ml 50
python dnabarcoder.py predict -i data/terrestrial_basal_fungi.fasta -c data/terrestrial_basal_fungi.classification -st 0.5 -et 1 -s 0.001 -rank class -prefix terrestrial_basal_fungi --removecomplexes yes -ml 50
python dnabarcoder.py predict -i data/zoosporic_basal_fungi.fasta -c data/zoosporic_basal_fungi.classification -st 0.5 -et 1 -s 0.001 -rank class -prefix zoosporic_basal_fungi -removecomplexes yes -ml 50
echo 'Finised class cutoff predictions at:' $(date)

################################################################################
# Predict phylum cutoffs
################################################################################

# Predict global similarity cutoffs at the phylum level
echo 'Predicting global similarity cutoffs for class at:' $(date)
python dnabarcoder.py predict -i data/unite2024ITS1.unique.fasta -c data/unite2024ITS1.unique.classification -st 0.5 -et 1 -s 0.001 -rank class -prefix fungi -minseqno -removecomplexes yes -ml 50
echo 'Finised phylum cutoff predictions at:' $(date)

# Deactivate the conda environment
conda deactivate

################################################################################
# Vizualise global cutoffs
################################################################################

python dnabarcoder.py predict -i data/unite2024ITS1.unique.fasta -c data/unite2024ITS1.unique.classification -rank species,genus,family,order,class,phylum
python dnabarcoder.py predict -i data/dikarya_fungi.fasta -c data/dikarya_fungi.classification -rank species,genus,family,order,class
python dnabarcoder.py predict -i data/terrestrial_basal_fungi.fasta -c data/terrestrial_basal_fungi.classification -rank species,genus,family,order,class
python dnabarcoder.py predict -i data/zoosporic_basal_fungi.fasta -c data/zoosporic_basal_fungi.classification -rank species,genus,family,order,class

conda deactivate