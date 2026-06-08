set shell := ['bash', '-ceuo', 'pipefail']

cmake_install_prefix := "$CONDA_PREFIX"
cmake_build_type := "Release"
cmake_build_dir := "cpp/build"

@default: build

@help:
    echo 'Usage: "just [cmake_install_prefix=<CMAKE_INSTALL_PREFIX>] [cmake_build_dir=<BUILD_DIR>] [cmake_build_type=Release]"'
    echo 'The variables will be passed to CMake for the C++ build.'
    echo 'Defaults:'
    echo "  cmake_install_prefix=\$CONDA_PREFIX (which is currently set to \"$CONDA_PREFIX\")"
    echo '  cmake_build_type=Release'
    echo '  cmake_build_dir=cpp/build (this is relative to the PETSIRD folder)'
    echo 'Run "just --summary" for possible recipes (default recipe is "build")'

@petsird:
    cd PETSIRD; \
    just cmake_install_prefix={{cmake_install_prefix}} build

@configure: petsird
   cmake -GNinja -S cpp -B {{cmake_build_dir}} \
      -DCMAKE_BUILD_TYPE:BOOL={{cmake_build_type}} \
      -DCMAKE_INSTALL_PREFIX:PATH={{cmake_install_prefix}}

@build: configure
   cd {{cmake_build_dir}} && \
    cmake --build . --config {{cmake_build_type}} && \
    cmake --install .


