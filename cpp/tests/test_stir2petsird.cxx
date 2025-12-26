
#include "stir/RunTests.h"
#include "stir/ProjDataInMemory.h"
#include "stir/Scanner.h"
#include "stir/PETSIRDInfo.h"

#include "convertor_helpers.h"
#include "petsird_helpers/create.h"

class STIR2PETSIRD_Tests : public stir::RunTests
{
public:
    void run_tests() override;

private: 
    void test_stir_to_petsird_conversion();
};

using namespace stir;
using namespace std; 

void STIR2PETSIRD_Tests::run_tests()
{
    test_stir_to_petsird_conversion();
}

void 
STIR2PETSIRD_Tests::test_stir_to_petsird_conversion()
{
    // Create small STIR projdata In memory 
    auto stir_scanner_sptr = make_shared<Scanner>(Scanner::E953);
    shared_ptr<ProjDataInfo> proj_data_info_sptr(ProjDataInfo::ProjDataInfoCTI(stir_scanner_sptr,
                                                                             /*span*/ 1,
                                                                             2,
                                                                             stir_scanner_sptr->get_num_detectors_per_ring() / 2,
                                                                             stir_scanner_sptr->get_max_num_non_arccorrected_bins(),
                                                                             /*arc_corrected*/ false));
    std::cout << "Start check" << std::endl;
    check_id_conversion(stir_scanner_sptr.get());
    std::cout << "Done" << std::endl;

    auto exam_info_sptr = std::make_shared<ExamInfo>();
    exam_info_sptr->set_calibration_factor(5.0F);
    exam_info_sptr->set_low_energy_thres(350.0F);
    exam_info_sptr->set_high_energy_thres(650.0F);
    auto proj_data_sptr = std::make_shared<ProjDataInMemory>(exam_info_sptr, proj_data_info_sptr);

    for (int seg_num = proj_data_info_sptr->get_min_segment_num(); seg_num <= proj_data_info_sptr->get_max_segment_num(); ++seg_num)
    {
        auto seg = proj_data_sptr->get_segment_by_view(seg_num);
        seg.fill(static_cast<float>(seg_num));
        proj_data_sptr->set_segment(seg);
    }


   petsird::Header header = get_dummy_header();
   auto& scanner_info = header.scanner;

   // force a single module type because STIR does not support multiple module types
   const auto num_types_of_modules = 1;
  petsird_helpers::create::initialize_scanner_information_dimensions(scanner_info, num_types_of_modules,false, false);
  set_scanner_geometry(scanner_info, *proj_data_info_sptr, *exam_info_sptr);

  std::cerr << "STIR views: " << proj_data_info_sptr->get_num_views() << '\n';
// std::cerr << "PETSIRD views: " << petsird_pdi.num_views() << '\n';
   // Now scanner_info has the full PETSIRD conversion. Let's read it in STIR
   std::cerr << "Now scanner_info has the full PETSIRD conversion. Let's read it in STIR" << std::endl;
  auto petsird_info_sptr = std::make_shared<PETSIRDInfo>(header);
  std::cerr << "Reading back scanner info from PETSIRD header" << std::endl;
 
}

USING_NAMESPACE_STIR
int main()
{
    STIR2PETSIRD_Tests tests;
    tests.run_tests();
    return tests.main_return_value();
}