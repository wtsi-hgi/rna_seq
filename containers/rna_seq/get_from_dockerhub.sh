#!/usr/bin/env bash

mkdir -p cache_dir
export SINGULARITY_CACHEDIR=$PWD/cache_dir

mkdir -p tmp_dir
export TMPDIR=$PWD/tmp_dir

# see https://github.com/wtsi-hgi/nextflow_rna_seq_container
rm -f rna_seq_1.0.sif
singularity pull docker://wtsihgi/rna_seq:1.0

rm -r cache_dir
rm -r tmp_dir
