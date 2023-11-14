#!/bin/bash
# FILEPATH: /workspaces/STIR2PETSIRD/stir/install_stir.sh

conda install boost

STIR_DIR=/workspaces/STIR2PETSIRD/stir

if [ ! -d "$STIR_DIR" ]; then
    mkdir -p "$STIR_DIR"
fi

cd ${STIR_DIR}
git clone https://github.com/UCL/STIR.git
mkdir -p build install
cd build

CMAKE_OPTIONS="-DCMAKE_INSTALL_PREFIX=${STIR_DIR}/install -DBUILD_SWIG_PYTHON=OFF"
cmake ${STIR_DIR}/STIR/ ${CMAKE_OPTIONS}
ninja install
