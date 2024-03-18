#!/bin/bash

#
# This is a helper script for installing the NVIDIA HPC SDK 
# Note: The script currently assumes Linux_x86_64 platform.
#

# Install the NVIDIA HPC SDK

# This link may need to be updated, as new compiler versions are released
# Note: Verified that it works with v24.1
if [ "x${NV_URL}" == "x" ]; then
  NV_URL=https://developer.download.nvidia.com/hpc-sdk/24.1/nvhpc_2024_241_Linux_x86_64_cuda_multi.tar.gz
fi

echo "Downloading the NVIDIA HPC SDK"
curl -s "${NV_URL}" | tar xpzf -

echo "Installing NVIDIA HPC SDK"

export NVHPC_INSTALL_DIR=$PWD/hpc_sdk
export NVHPC_SILENT=true

(cd nvhpc_*; ./install)

# create helper scripts
mkdir setup_scripts
cat > setup_scripts/setup_nv_hpc_bins.sh << EOF
module add ./hpc_sdk/modulefiles/nvhpc-openmpi3/*

EOF

echo "Setup script avaiabile in $PWD/setup_scripts/setup_nv_hpc_bins.sh"
