
# requirements.txt

numpy
matplotlib
statsmodels
scipy

# Dockerfile

# Use an official Python runtime as a parent image
FROM python:3.9-slim

# Set the working directory in the container
WORKDIR /usr/src/app

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base \
    r-base-dev \
    && rm -rf /var/lib/apt/lists/*

# First, copy just the files defining dependencies.
# It ensures that the Docker cache doesn't get invalidated unless dependencies change
COPY requirements.txt .
COPY install_packages.R .

# Install any needed Python packages specified in requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Install any needed R packages
RUN Rscript install_packages.R

# Now copy the rest of your application's code.
# This way, the earlier layers (dependencies) are cached and won't be rebuilt unless the requirements files change.
# COPY src/ .

# execute ls (Not usually needed in production Dockerfile, used here just for verifying contents)
# Commenting out, as it's not typically needed for production images
# RUN ls -la

# Copy the entrypoint script into the container
COPY entrypoint.sh .

# Adjust file permissions to ensure the script is executable
RUN chmod +x entrypoint.sh

# Use the entrypoint script to run commands
# Make sure to remove the comment to enable the ENTRYPOINT
# ENTRYPOINT ["./entrypoint.sh"]

# The CMD should provide default parameters for the ENTRYPOINT
CMD ["bash"]



# README.md

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


# install_packages.R

install.packages(c("dplyr","data.table", "factoextra", "funrar", "tidyverse"), repos = "https://cloud.r-project.org/")


# doc

