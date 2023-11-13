# STIR2PETSIRD converter

The purpose of this repo is to provide a starting point for converting list mode PET data read by STIR to PETSIRD.

## Background
The [Emission Tomography Standardization Initiative (ETSI)](https://etsinitiative.org/)
is working towards establishing a standard for PET Raw Data, called PETSIRD ("PET ETSI Raw Data").

The specification uses the [yardl](https://aka.ms/yardl) tool to define the model. `yardl` can be
used to read the specification (in the `model` directory) and generate an PI for both C++ and API to read/write PETSIRD data.

Currently, a draft model is defined in https://github.com/ETSInitiative/PRDdefinition.

CAVEAT: the draft model generates code in the `prd` namespace. Nevertheless, we have used the name PETSIRD here
in most places (except where needed).

### Using your repo

1. Open ***your*** repo in [GitHub Codespaces](https://code.visualstudio.com/docs/remote/codespaces) or
in a [VS Code devcontainer](https://code.visualstudio.com/docs/devcontainers/containers).
This codespace/container will contain all necessary tools, including `yardl` itself, as well as your repository.
2. Use `yardl` to generate C++ and Python code for the model:
  ```sh
  cd YourRepoName
  cd PETSIRD/model
  yardl generate
  cd ../..
  ```
3. Start working in either the [`cpp`](cpp/README.md) and/or [`python`](python/README.md) directories.
