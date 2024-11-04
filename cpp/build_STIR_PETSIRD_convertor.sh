#! /usr/bin/env bash
# Utility script to help build STIR_PETSIRD_convertor in the cpp dir

# get current directory, see https://stackoverflow.com/a/4774063
FILEPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

# use this to point to your local installation of CMake doesn't find it
#CMAKE_STIR_OPTS="-DSTIR_DIR=\"${HOME}/stir/install/lib/cmake/STIR-6.2\""
cpp_dir="${FILEPATH}"
BUILD_DIR="${cpp_dir}/build"

pushd "${cpp_dir}/../PETSIRD/model"
yardl generate
popd

# rm -rf ${BUILD_DIR}  # Optional
mkdir -p "${BUILD_DIR}"

cd ${BUILD_DIR}
cmake -G Ninja -S ${cpp_dir} -DHDF5_ROOT=${CONDA_PREFIX} -DCMAKE_PREFIX_PATH=${CONDA_PREFIX} ${CMAKE_STIR_OPTS}
ninja
