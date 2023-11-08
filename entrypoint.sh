#!/bin/bash
# entrypoint.sh

# Run a Python script
python script.py "$@"

# Run an R script
Rscript script.R

# Keep the container running (if you need it to stay up after scripts)
# tail -f /dev/null
