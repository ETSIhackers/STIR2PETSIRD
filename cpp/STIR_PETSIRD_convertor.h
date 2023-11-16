// STIR_PETSIRD_convertor.h
#ifndef STIRPETSIRDCONVERTOR_H
#define STIRPETSIRDCONVERTOR_H

#include <string>


class STIRPETSIRDConvertor {
public:
    STIRPETSIRDConvertor(const std::string& out_filename, const std::string& in_filename) : out_filename(out_filename), in_filename(in_filename) {}
    void process_data();

private:
    // Output file name to write to, will delete previous file if it exists
    std::string out_filename;
    // Input file name to read from
    std::string in_filename;
};

#endif // STIRPETSIRDCONVERTOR_H
