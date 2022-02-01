#!/bin/bash

#SBATCH --job-name=variantcalling
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=144:00:00
#SBATCH --mem-per-cpu=64000M
#SBATCH --output=/fast/users/altwassr_c/slurm_logs/%x.%j.out
#SBATCH --error=/fast/users/altwassr_c/slurm_logs/%x.%j.err

echo 'Start'
snakemake -r --nt --jobs 40 --use-conda -p --rerun-incomplete --conda-prefix=/fast/users/altwassr_c/work/conda-envs/ \
--skip-script-cleanup \
--drmaa " --error=/fast/users/altwassr_c/scratch/slurm_logs/preprocessing.%j.err \
--output=/fast/users/altwassr_c/scratch/slurm_logs/preprocessing.%j.out \
--time=24:00:00 \
--partition=medium \
--mem=160000 \
--mem-per-cpu=150480 \
--ntasks-per-node=4"

# --keep-going \
echo 'Finished'
