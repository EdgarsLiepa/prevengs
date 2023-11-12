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

