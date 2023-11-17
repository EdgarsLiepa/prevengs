# Prevengs

Pediatric oncogenomics data analysis pipeline for RTU Innovation Health Hub med tech development program.

THIS IS A DEVELOPMENT VERSION OF THE PIPELINE. IT IS NOT READY FOR PRODUCTION USE.

## Description

dependency list.

- Python
    - Numpy
    - MatplotLib
	- R
    - Deseq2

Input Files:

- Sample Counts
- Reference Lengths

Output:

- Sample Reference table
- PCA
- Top10 genes plot.
- TPM list if calculated.
- QQplots.

## About The Project

## Getting Started

add dependencies to python dependency list to requirements.txt

add dependencies to R library install code to requirements.R which is executed while image is built.

Build the docker image

```sh
docker build -t prevengs .
```

Run docker image with full pipeline

```sh
TODO add command
```

Run the python script through the docker image.

```sh
docker run -v "$PWD":/usr/src/app -it --rm prevengs python src/script.py 'data/RNS_FLT3_156.F.fastq.genome.htseq_counts.txt' 'data/gencode.v31.chr_patch_hapl_scaff.annotation.gtf'
```

**$PWD** is the current working directory. It is mounted to the docker container as /usr/src/app

Run docker container in interactive mode

```sh
docker run -v "$PWD":/usr/src/app -it --rm prevengs
```



This code processes transcriptome data and generates a plot of the results.
The output plot shows the expression levels of different genes.
![Transkriptoma_datu_plusma.jpg](doc/Transkriptoma_datu_plusma.jpg)

## Nextflow

## TODO

- [ ] Add output parameters

## Authors

- Edgars Liepa
- Ņikita Fomins
- Pauls Daugulis
- Agate Jarmakovica
- Aivija Stugle
