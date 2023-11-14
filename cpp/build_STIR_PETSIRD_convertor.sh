#! /bin/bash
# FILEPATH: /workspaces/STIR2PETSIRD/cpp/build.sh

# Utility script to build STIR_PETSIRD_convertor in the cpp dir
STIR_DIR="/workspaces/STIR2PETSIRD/stir/install"
PETSIRD_DIR="/workspaces/STIR2PETSIRD/cpp"
cpp_dir="/workspaces/STIR2PETSIRD/cpp"

# rm -rf ${cpp_dir}/build && mkdir -p ${cpp_dir}/build
cd ${cpp_dir}/build
cmake -G Ninja -S .. -DHDF5_ROOT=${CONDA_PREFIX} -DSTIR_DIR=${STIR_DIR} -DCMAKE_PREFIX_PATH=${PETSIRD_DIR}
ninja
