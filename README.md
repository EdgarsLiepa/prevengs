# prevengs

Pediatric oncogenomics data analysis pipeline for RTU Innovation Health Hub med tech development programm

dependency list.

- Python
    - Numpy
    - MatplotLib
	- R
    - Deseq2

Input Faili:

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

Run the docker image
```sh
docker run -it --rm prevengs '/Users/edgars/Projects/prevengs/data/RNS_FLT3_156.F.fastq.genome.htseq_counts.txt'
```

This code processes transcriptome data and generates a plot of the results.
The output plot shows the expression levels of different genes.
![Transkriptoma_datu_plusma.jpg](doc/Transkriptoma_datu_plusma.jpg)

## Authors

Edgars Liepa
Ņikita Fomins
Pauls Daugulis
Agate Jarmakovica
Aivija Stugle