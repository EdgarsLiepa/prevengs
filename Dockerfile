# Use an official Python runtime as a parent image
FROM ghcr.io/rocker-org/r-ver:latest

# Set the working directory in the container
WORKDIR /usr/src/app

ENV SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True

# First, copy just the files defining dependencies.
# It ensures that the Docker cache doesn't get invalidated unless dependencies change
COPY install_tidyverse.sh .
COPY requirements.txt .
COPY util/outpyr util/outpyr
COPY util/outsingle/requirements.txt util/outsingle/requirements.txt

# Install any needed Python packages specified in requirements.txt
RUN ./install_tidyverse.sh

# install pip
RUN apt-get update && apt-get install -y python3-pip
RUN pip install --no-cache-dir -r requirements.txt


# install R package to run python
RUN pip install util/outpyr
RUN R -e "install.packages('reticulate')"
RUN R -e "BiocManager::install(\"DESeq2\")"
RUN pip install Cython
RUN pip install -r util/outsingle/requirements.txt
RUN pip install scikit-learn 
# RUN pip install gene-outlier-detection


# Install any needed R packages
# COPY install_packages.R .
# RUN Rscript install_packages.R

# Now copy the rest of your application's code.
# This way, the earlier layers (dependencies) are cached and won't be rebuilt unless the requirements files change.
# COPY src/ .

# execute ls (Not usually needed in production Dockerfile, used here just for verifying contents)
# Commenting out, as it's not typically needed for production images
# RUN ls -la

# Copy the entrypoint script into the container
# COPY entrypoint.sh .

# Adjust file permissions to ensure the script is executable    
# RUN chmod +x entrypoint.sh

# Use the entrypoint script to run commands
# Make sure to remove the comment to enable the ENTRYPOINT
# ENTRYPOINT ["./entrypoint.sh"]

# The CMD should provide default parameters for the ENTRYPOINT
CMD ["bash"]

