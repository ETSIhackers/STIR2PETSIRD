#! /bin/bash
# FILEPATH: /workspaces/STIR2PETSIRD/cpp/build.sh

# Utility script to build STIR_PETSIRD_convertor in the cpp dir
STIR_DIR="/workspaces/STIR2PETSIRD/stir/install"
PETSIRD_DIR="/workspaces/STIR2PETSIRD/cpp"
cpp_dir="/workspaces/STIR2PETSIRD/cpp"

# rm -rf ${cpp_dir}/build
if [ ! -d "${cpp_dir}/build" ]; then
    mkdir -p "${cpp_dir}/build"
fi

cd ${cpp_dir}/build
cmake -G Ninja -S ${cpp_dir} -DHDF5_ROOT=${CONDA_PREFIX} -DSTIR_DIR=${STIR_DIR} -DCMAKE_PREFIX_PATH=${PETSIRD_DIR}
ninja
