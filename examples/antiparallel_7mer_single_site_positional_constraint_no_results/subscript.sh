#!/bin/bash

# Choose ONE of the following methods to run PepAD.

# ============================================================
# Method 1: Run PepAD using the Apptainer image
# ============================================================

module load apptainer                                   # If required on the local system
PEPAD_SIF=/path/to/PepAD_container/PepAD_package.sif    # Required: path to PepAD_package.sif
RUN_DIR="$(pwd -P)"                                     # Absolute path of the current working directory
apptainer run --bind "${RUN_DIR}:/work" --pwd /work "$PEPAD_SIF"

# ============================================================
# Method 2: Run PepAD after manual compilation
# ============================================================
# module load PrgEnv-intel                                # If required on the local system
# export PATH="../../src/PepAD:$PATH"                     # Directory containing the PepAD executable
# PepAD