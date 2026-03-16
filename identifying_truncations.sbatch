#!/bin/bash

#SBATCH --job-name=missing_domains
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200gb
#SBATCH --time=12:00:00
#SBATCH --partition=short
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

#Activate Biopython environment
source /Users/anor9792/biopython/biopython/bin/activate

# Execute the Python script
python identifying_truncations.py


echo $(date +"[%b %d %H:%M:%S] Done")



