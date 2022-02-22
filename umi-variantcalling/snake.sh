#!/bin/bash

#SBATCH --job-name=variantcalling
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=144:00:00
#SBATCH --mem-per-cpu=24000M
#SBATCH --output=/fast/users/altwassr_c/scratch/slurm_logs/%x.%j.out
#SBATCH --error=/fast/users/altwassr_c/scratch/slurm_logs/%x.%j.err

echo 'Start'
snakemake -r --nt --jobs 20 --use-conda -p --rerun-incomplete --conda-prefix=/fast/users/altwassr_c/work/conda-envs/ \
--drmaa " --error=/fast/users/altwassr_c/scratch/slurm_logs/variantcalling.%j.err \
--output=/fast/users/altwassr_c/scratch/slurm_logs/variantcalling.%j.out \
--time=02:00:00 \
--partition=short \
--mem=16000 \
--mem-per-cpu=15048 \
--ntasks-per-node=4"
# --skip-script-cleanup \
# --keep-going \

# --until annovar \
echo 'Finished'
