import sys

# Currently hard-wire location of generated files. This will need to change!
STIR2PETSIRD_PATH = "/workspaces/STIR2PETSIRD"
import stir
import petsird

data_path = f"{STIR2PETSIRD_PATH}/data/root_header_test1.hroot"

lm_data = stir.ListModeData.read_from_file(data_path)

print("Done")
