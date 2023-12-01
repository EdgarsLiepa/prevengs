# Use an official Python runtime as a parent image
FROM ghcr.io/rocker-org/r-ver:latest

# Set the working directory in the container
WORKDIR /usr/src/app


COPY install_tidyverse.sh .
RUN ./install_tidyverse.sh



# The CMD should provide default parameters for the ENTRYPOINT
CMD ["bash"]

