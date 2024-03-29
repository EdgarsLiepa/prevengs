# PreveNGS

Pediatric oncogenomics data analysis pipeline for RTU Innovation Health Hub med tech development program.

THIS IS A DEVELOPMENT VERSION OF THE PIPELINE. IT IS NOT READY FOR PRODUCTION USE.

## Description

dependencies list.

- Docker.
- [outsingle](https://github.com/esalkovic/outsingle)
- N-of-1 [gene_outlier_detection](https://github.com/jvivian/gene-outlier-detection)
- Desq2 [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- testthat

Input Files:

- HT-Seq Sample Counts - SampleName.htseq_counts.txt
    - HT-Seq output file with counts for each gene.
    - 2 collumns
        - Gene ID
        - Count
- Reference Lengths
- metadata file.

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

clone this repository and update submodules

```sh
git clone https://github.com/EdgarsLiepa/prevengs.git

# Update repository submodules.
git submodule update --init --recursive
```

### Download docker image from Docker Hub

```sh
docker pull edgarsliepa/prevengs:latest 
```

### Get Human Release 31 (GRCh38.p12) reference gene anotations

``` bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.basic.annotation.gtf.gz 
```

### Run image from Docker Hub


Example code to run Rscript using edgarsliepa/prevengs:latest image.

``` bash
docker run -v "$PWD":/usr/src/app edgarsliepa/prevengs:latest  Rscript src/pipeline.R data/BKUS_SAMPLES data/gencode.v31.chr_patch_hapl_scaff.annotation.gtf ./rez data/metadata_BKUS.tsv
```

## Build the docker image

Docker image can be built from Docker file provided in the repository.
From prevengs project.

```sh
docker build -t prevengs .
```

corscompile docker image with buildx

```sh
docker buildx build --platform linux/amd64,linux/arm64,linux/arm/v7 -t <username>/<image>:latest --push .
```

### Run built docker image

Run the R pipeline script through the docker image.

```sh
# Usage: Rscript pipeline.R <input_directory> <output_folder>

docker run -v "$PWD":/usr/src/app -it --rm prevengs Rscript src/pipeline.R data/BKUS_SAMPLES data/gencode.v31.chr_patch_hapl_scaff.annotation.gtf ./rez data/metadata_BKUS.tsv
```

Change **$PWD** to directory path with HTseqfiles and scripts. Local directory is mounted to the docker container at /usr/src/app.  
This needed to access the files from the docker container.

Currently the pipeline script needs to mounted as well.

### Run the python script through the docker image  

```bash
# Process transcriptome featureCounts.

# positional arguments:
#   counts_file  The featureCounts file to process.
#   gtf_file     The GTF file to use for gene length calculation.

# options:
#   -h, --help   show this help message and exit

docker run -v "$PWD":/usr/src/app -it --rm prevengs python3 src/script.py 'data/RNS_FLT3_156.F.fastq.genome.htseq_counts.txt' 'data/gencode.v31.chr_patch_hapl_scaff.annotation.gtf'
```

## Differential expression analysis

```sh
docker run -v "$PWD":/usr/src/app -it --rm prevengs Rscript src/dge.R path/to/counts/file.tsv path/to/metadata/file.tsv
```


### Run the R script through the docker image

Run the R top5_boxplot script through the docker image.

```sh

docker run -v "$PWD":/usr/src/app -it --rm prevengs Rscript src/top5_boxplot.R data/ ./
```

Run the R PCA script through the docker image.

Usage: Rscript PCA_for_all_genes.R <input_directory> <output_file>

```sh
docker run -v "$PWD":/usr/src/app -it --rm prevengs Rscript src/PCA_for_all_genes.R data/ ./
```


Run docker container in interactive mode

```sh
docker run -v "$PWD":/usr/src/app -it --rm prevengs
```

This code processes transcriptome data and generates a plot of the results.
The output plot shows the expression levels of different genes.
![Transkriptoma_datu_plusma.jpg](doc/Transkriptoma_datu_plusma.jpg)

## Run tests

To run the tests, run the following command:

```sh
docker run -v "$PWD":/usr/src/app -it --rm prevengs Rscript tests/dge_test.R
```

the test that package is used to run the tests.

## TODO

- [ ] Add exception when file in HT seq folder is not in right forma (create_reference_table.R)
    - [ ] Define HT seq file format.
    - [ ] BUG: breaks if htseq files are not named with *counts.txt
    - [ ] Add maybe patern match for file
- [ ] Fix Submodules.
- [ ] Definēt kādu references genome tabulu izmantot.

## Authors

- Edgars Liepa
- Ņikita Fomins
- Pauls Daugulis
- Agate Jarmakovica
- Aivija Stugle
