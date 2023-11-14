#! /bin/bash
# FILEPATH: /workspaces/STIR2PETSIRD/cpp/build.sh

# Utility script to help build STIR_PETSIRD_convertor in the cpp dir
STIR2PETSIRD_DIR="/workspaces/STIR2PETSIRD"
STIR_DIR="${STIR2PETSIRD_DIR}/stir/install/lib/cmake/"
PETSIRD_DIR="${STIR2PETSIRD_DIR}/PETSIRD/cpp/generated"
cpp_dir="${STIR2PETSIRD_DIR}/cpp"
BUILD_DIR="${cpp_dir}/build"

# rm -rf ${BUILD_DIR}  # Optional
if [ ! -d "${BUILD_DIR}" ]; then
    mkdir -p "${BUILD_DIR}"
fi

cd ${BUILD_DIR}
cmake -G Ninja -S ${cpp_dir} -DHDF5_ROOT=${CONDA_PREFIX} -DSTIR_DIR=${STIR_DIR} -DCMAKE_PREFIX_PATH=${PETSIRD_DIR}
ninja
