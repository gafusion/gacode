#!/bin/bash

#
# This is a helper script for installing the NVIDIA HPC SDK 
# Note: The script currently assumes Linux_x86_64 platform.
#

# Install the NVIDIA HPC SDK

# This link may need to be updated, as new compiler versions are released
# Note: Verified that it works with v25.1
NV_URL=https://developer.download.nvidia.com/hpc-sdk/25.1/nvhpc_2025_251_Linux_x86_64_cuda_12.6.tar.gz
NV_VER=25.1

echo "Downloading the NVIDIA HPC SDK"
curl -s "${NV_URL}" | tar xpzf -

echo "Installing NVIDIA HPC SDK"

export NVHPC_INSTALL_DIR=$PWD/hpc_sdk
export NVHPC_SILENT=true

(cd nvhpc_*; ./install)

# create helper scripts
mkdir setup_scripts
cat > setup_scripts/setup_nv_hpc_bins.sh << EOF
#module add ./hpc_sdk/modulefiles/nvhpc/*
export PATH=$NVHPC_INSTALL_DIR/Linux_x86_64/$NV_VER/compilers/bin:\$PATH
export PATH=$NVHPC_INSTALL_DIR/Linux_x86_64/$NV_VER/comm_libs/mpi/bin:\$PATH
export MANPATH=$NVHPC_INSTALL_DIR/Linux_x86_64/$NV_VER/compilers/man:\$MANPATH

EOF

echo "Setup script avaiabile in $PWD/setup_scripts/setup_nv_hpc_bins.sh"
