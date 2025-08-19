set shell := ['bash', '-ceuo', 'pipefail']

default: build

@ensure-build-dir:
    mkdir -p cpp/build

@generate:
    cd PETSIRD/model; \
    yardl generate

@configure: generate ensure-build-dir
    cd cpp/build; \
    cmake -GNinja ..

@build: generate configure
    cd cpp/build; \
    ninja


