#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate /lustre/scratch115/humgen/resources/conda/star

salmon index \
-t Homo_sapiens.GRCh38.ncrna_plus_cdna.fa \
-i salmon_index_Homo_sapiens.GRCh38.cdna_plus_ncrna.all \
-k 31

