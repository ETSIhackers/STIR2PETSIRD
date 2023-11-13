# PETSIRD basic C++ example

This directory provides a skeleton for C++ example code to read/write PETSIRD data. You need to `yardl generate` in the `model` directory first.

Check the example in the PRDdefinition repo first.

1. Compile the code:
   ```sh
   cd cpp
   mkdir -p build && cd build`
   cmake -G Ninja -S .. -DHDF5_ROOT=$CONDA_PREFIX
   ninja`
   ```
   If you did not use `conda` to install HDF5, do not add the `-DHDF5_ROOT=$CONDA_PREFIX` part of the `cmake` line.

