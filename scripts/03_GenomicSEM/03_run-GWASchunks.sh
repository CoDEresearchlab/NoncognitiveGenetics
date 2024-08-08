#!/bin/bash -l
#SBATCH --output=../out/GSEM_gwas.chol.chunk.array_%A_%a.out
#SBATCH --error=../out/GSEM_gwas.chol.chunk.array_%A_%a.err
#SBATCH --time=72:00:00
#SBATCH --ntasks=15
#SBATCH --nodes=1
#SBATCH --array=1-135
#SBATCH --mem=15GB
#SBATCH --partition=brc

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

job=$(sed -n "$SLURM_ARRAY_TASK_ID"p listGwasChunks)

module load apps/R/3.6.0

Rscript --vanilla GWAS.R $job
