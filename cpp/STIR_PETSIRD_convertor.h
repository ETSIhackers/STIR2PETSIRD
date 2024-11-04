// STIR_PETSIRD_convertor.h
#ifndef STIRPETSIRDCONVERTOR_H
#define STIRPETSIRDCONVERTOR_H

#include <string>
#include "stir/listmode/ListModeData.h"
#include "stir/listmode/ListRecord.h"

petsird::BoxSolidVolume get_crystal_template(const std::array<float, 3>& crystal_length);

class STIRPETSIRDConvertor {
public:
    STIRPETSIRDConvertor(const std::string& out_filename, const std::string& in_filename);
    void process_data();

private:
    // Output file name to write to, will delete previous file if it exists
    std::string out_filename;
    // Input file name to read from
    std::string in_filename;

    // Pointer to listmode data
    std::unique_ptr<stir::ListModeData> lm_data_ptr;
};

#endif // STIRPETSIRDCONVERTOR_H
