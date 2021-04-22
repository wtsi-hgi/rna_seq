echo starting bsub
rm -f hs_err_pid* \
    && rm -f timeline* \
    && rm -f trace* \
    && rm -rf report* \
    && rm -f nohup.out \
    && rm -f bsub.o \
    && rm -f bsub.e \
    && rm -f .nextflow.log && \
    bsub -G hgi \
	 -R'select[mem>100000] rusage[mem=100000] span[hosts=1]' \
	 -M 100000 \
	 -n 4 -o bsub.o -e bsub.e \
	 -q yesterday bash ./build_star_index.sh > bjob.id
echo finished bsub

cat bjob.id
