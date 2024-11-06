// STIR_PETSIRD_convertor.h
#ifndef STIRPETSIRDCONVERTOR_H
#define STIRPETSIRDCONVERTOR_H

#include <string>
#include <memory>

namespace stir {
  class ListModeData;
  class BinNormalisation;
}

class STIRPETSIRDConvertor {
public:
    STIRPETSIRDConvertor(const std::string& out_filename, const std::string& in_filename);
    void process_data();
    void set_normalisation_sptr(const std::shared_ptr<stir::BinNormalisation> norm_sptr)
    {
      this->normalisation_sptr = norm_sptr;
    }
private:
    // Output file name to write to, will delete previous file if it exists
    std::string out_filename;
    // Input file name to read from
    std::string in_filename;

    // Pointer to listmode data
    std::unique_ptr<stir::ListModeData> lm_data_ptr;
    std::shared_ptr<stir::BinNormalisation> normalisation_sptr;
};

#endif // STIRPETSIRDCONVERTOR_H
