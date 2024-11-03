# PETSIRD basic C++ example

This directory provides a C++ example code to read/write PETSIRD data from STIR. You need to `yardl generate` in the `model` directory first.

## Requirements
   - STIR
   - PETSIRD

## Compile
```
cd cpp
mkdir -p build && cd build`
cmake -G Ninja -S ${cpp_dir} -DHDF5_ROOT=${CONDA_PREFIX} -DSTIR_DIR=${STIR_DIR} -DCMAKE_PREFIX_PATH=${PETSIRD_DIR}
ninja
```
If you did not use `conda` to install HDF5, do not add the `-DHDF5_ROOT=$CONDA_PREFIX` part of the `cmake` line.

## Usage
```
STIR_PETSIRD_convertor <output filename> <input filename>
```
where `output filename` is the HDF5 filename to write to and `input filename` is the STIR readable list mode file.
