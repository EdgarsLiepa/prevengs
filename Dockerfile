# Use an official Python runtime as a parent image
FROM python:3.9-slim

# Set the working directory in the container
WORKDIR /usr/src/app

# Install R and some dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base \
    r-base-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy the current directory contents into the container at /usr/src/app
COPY src/ /usr/src/app

# Install any needed Python packages by specifying requirements.txt
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Install any needed R packages
# You can automate the installation of R packages by listing them in a R script
# Example install_packages.R content: 
# install.packages(c("ggplot2", "dplyr"), repos = "https://cloud.r-project.org/")
COPY install_packages.R ./
RUN Rscript install_packages.R


# execute ls
RUN ls -la

# Run Python script (make sure you have a script named `script.py` in your directory)
COPY entrypoint.sh /usr/src/app/
# ENTRYPOINT ["./entrypoint.sh"]
# CMD ["python", "/src/script.py"]

# If you also want to run an R script, you might need to modify the CMD to execute a shell script that runs both, or use multiple CMDs in an entrypoint script.
