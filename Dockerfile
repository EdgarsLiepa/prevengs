# Use an official Python runtime as a parent image
FROM ghcr.io/rocker-org/r-ver:latest

# Set the working directory in the container
WORKDIR /usr/src/app



# First, copy just the files defining dependencies.
# It ensures that the Docker cache doesn't get invalidated unless dependencies change
COPY install_tidyverse.sh .
COPY requirements.txt .

# Install any needed Python packages specified in requirements.txt
RUN ./install_tidyverse.sh

# install pip
RUN apt-get update && apt-get install -y python3-pip
RUN pip install --no-cache-dir -r requirements.txt

# The CMD should provide default parameters for the ENTRYPOINT
CMD ["bash"]

