set shell := ['bash', '-ceuo', 'pipefail']

install_prefix := "$CONDA_PREFIX"

default: (build install_prefix)

@help:
    echo Usage: "just build [<CMAKE_INSTALL_PREFIX>]"
    echo The third argument is optional and will default to \$CONDA_PREFIX
    echo "which is currently set to $CONDA_PREFIX"

@ensure-build-dir:
    mkdir -p cpp/build

@generate:
    cd PETSIRD/model; \
    yardl generate

@configure install_prefix: generate ensure-build-dir
    cd cpp; \
    cmake -GNinja -S . -B build/ -DCMAKE_INSTALL_PREFIX:PATH={{install_prefix}}

@build install_prefix: generate (configure install_prefix)
    cd cpp/build; \
    ninja install


