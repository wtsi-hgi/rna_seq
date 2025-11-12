# rna_seq_container

Dockerfiles for images used in the HGI RNA-seq pipeline.

## Version history
* The current version of the container used by HGI `rna-seq` pipeline is `rna_seq_1.2`.
* The `rna_seq` Dockerfile was taken from https://github.com/wtsi-hgi/nextflow_rna_seq_container. 
  It should match the latest container version. However, this has not been verified. 
* The `tximport` Dockerfile created form scratch because TxImport script form the container environment 
  doesn't work with the latest Ensembl.

## Availability

All containers are available at Dockerhub and can be pulled to farm:
```bash
singularity pull docker://wtsihgi/rna_seq:1.2
```


## Software versions  
The versions are saved during docker build in container file `/usr/conda_software_versions.txt`:  
(`docker run wtsihgi/rna_seq:1.2 cat /usr/conda_software_versions.txt`)

docker tag **1.2**:
No info recorded, TBD.

docker tag **1.1** adds fastqc and updates Python:

```
FastQC v0.11.9
Python 3.9.4
```

For docker tag **1.0**:
```
featureCounts v2.0.2
salmon 1.4.0
STAR 2.7.8a
multiqc, version 1.10.1
QTLtools_1.2_Ubuntu16.04_x86_64
R /opt/conda/envs/conda_rna_seq/bin/R
R version 4.0.3 with libpath "/opt/conda/envs/conda_rna_seq/lib/R/library"
Python 3.9.2 (/opt/conda/envs/conda_rna_seq/bin/python)
```

## Building images

### Dockerhub auto-build:  
see https://hub.docker.com/repository/docker/wtsihgi/rna_seq

### RNASeq container build on Sanger farm

The docker containers building on the Farm don't work.
Therefore, to build containers, you should use an OpenStack VM
with Singularity and Docker installed.

```bash
docker build --tag local/rna_seq_tximport:1.2 .
singularity build rna_seq_tximport_1.2.sif docker-daemon://local/rna_seq_tximport:1.2
```


### TxImport manual build and push

```bash
docker build --tag wtsihgi/rna_seq_tximport:1.2 .
docker image ls
#docker image tag e682162f5052 wtsihgi/rna_seq:1.1
```

### rna_seq Docker manual build:

```bash
docker build . # or docker build --tag wtsihgi/rna_seq:1.1 .
docker image ls
docker image tag e682162f5052 wtsihgi/rna_seq:1.1
docker login
docker image push wtsihgi/rna_seq:1.1

# check that conda env is loaded on 'run":
docker run wtsihgi/rna_seq:1.1 conda env list

# check software versions:
docker run wtsihgi/rna_seq:1.1 cat /usr/conda_software_versions.txt

# check that path has conda env bin dir first:
docker run wtsihgi/rna_seq:1.1 printenv
```

### Convert docker image to singularity:

```
## option 1 (pull first):
singularity pull docker://wtsihgi/rna_seq:1.1
# this creates image file rna_seq_1.1.sif in current dir.

# check conda env is loaded by default (requires --containall):
singularity exec --containall rna_seq_1.1.sif conda env list

## option 2 (single command, exec directly from Dockerhub):
singularity exec --containall docker://wtsihgi/rna_seq:1.1 conda env list

## option 3 (user Docker to create singularity image):
export IMAGE=wtsihgi/rna_seq:1.1
mkdir -p ~/singu &&  rm -rf singu/*.sif
docker run -v /var/run/docker.sock:/var/run/docker.sock -v ~/singu:/output --privileged -t --rm quay.io/singularity/docker2singularity $IMAGE
# check image:
singularity shell --containall singu/wtsihgi_rna_seq_1.1.sif  conda env list
```
