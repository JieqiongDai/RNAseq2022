#!/bin/bash
mkdir -p log
sbatch wraaper/sbatch_snakemake.batch 2>log/sbatch.err
