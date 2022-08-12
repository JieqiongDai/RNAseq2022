#!/bin/bash
mkdir -p log
qsub wrapper/qsub_snakemake.batch 2>log/qsub.err
