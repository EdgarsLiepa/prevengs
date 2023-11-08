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
docker run -it --rm prevengs 'asd'
```
