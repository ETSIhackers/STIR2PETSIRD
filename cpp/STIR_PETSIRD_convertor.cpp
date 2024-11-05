#include <iostream>
#include <string>
#include <fstream>
#include <memory>

#include "stir/listmode/ListModeData.h"
#include "stir/listmode/ListEvent.h"
#include "stir/listmode/ListRecord.h"
#include "stir/listmode/CListEventCylindricalScannerWithDiscreteDetectors.h"
#include "stir/LORCoordinates.h"
#include "stir/Succeeded.h"
#include "stir/IO/read_from_file.h"
#include "stir/error.h"
#include "binary/protocols.h"
#include "petsird_helpers.h"

#include "stir/ProjDataInfoBlocksOnCylindrical.h"
#include "stir/ProjDataInfoCylindrical.h"

#include "STIR_PETSIRD_convertor.h"


//! return a cuboid volume
petsird::BoxSolidVolume get_crystal_template
(const std::array<float, 3> &  crystal_length)
{
  using petsird::Coordinate;
  petsird::BoxShape crystal_shape{ Coordinate{ { 0, 0, 0 } },
                                   Coordinate{ { 0, 0, crystal_length[2] } },
                                   Coordinate{ { 0, crystal_length[1], crystal_length[2] } },
                                   Coordinate{ { 0, crystal_length[1], 0 } },
                                   Coordinate{ { crystal_length[0], 0, 0 } },
                                   Coordinate{ { crystal_length[0], 0, crystal_length[2] } },
                                   Coordinate{ { crystal_length[0], crystal_length[1], crystal_length[2] } },
                                   Coordinate{ { crystal_length[0], crystal_length[1], 0 } } };

  petsird::BoxSolidVolume crystal{ crystal_shape, /* material_id */ 1 };
  return crystal;
}

//! return a module of NUM_CRYSTALS_PER_MODULE cuboids
petsird::DetectorModule get_detector_module_tmpl(const std::array<float, 3> &  crystal_length, 
                                                const std::array< int, 3> & NUM_CRYSTALS_PER_MODULE, 
                                                const float RADIUS)
{
  petsird::ReplicatedBoxSolidVolume rep_volume;
  {
    rep_volume.object = get_crystal_template(crystal_length);
    const auto N0 = NUM_CRYSTALS_PER_MODULE[0];
    const auto N1 = NUM_CRYSTALS_PER_MODULE[1];
    const auto N2 = NUM_CRYSTALS_PER_MODULE[2];
    for (int rep0 = 0; rep0 < N0; ++rep0)
      for (int rep1 = 0; rep1 < N1; ++rep1)
        for (int rep2 = 0; rep2 < N2; ++rep2)
          {
            petsird::RigidTransformation transform{ { { 1.0, 0.0, 0.0, RADIUS + rep0 * crystal_length[0] },
                                                      { 0.0, 1.0, 0.0, (rep1 - N1 / 2) * crystal_length[1] },
                                                      { 0.0, 0.0, 1.0, (rep2 - N2 / 2) * crystal_length[2] } } };
            rep_volume.transforms.push_back(transform);
            rep_volume.ids.push_back(rep0 + N0 * (rep1 + N1 * rep2));
          }
  }

  petsird::DetectorModule detector_module;
  detector_module.detecting_elements.push_back(rep_volume);
  detector_module.detecting_element_ids.push_back(0);

  return detector_module;
}

// Convert from STIR scanner to petsird scanner info (for now, just cylindrical non-TOF scanners)
petsird::ScannerInformation
get_scanner_geometry(const stir::ProjDataInfo& stir_proj_data_info, 
const stir::ExamInfo& stir_exam_info)
{
  const stir::Scanner* stir_scanner = stir_proj_data_info.get_scanner_ptr(); 
  const float radius = stir_scanner->get_inner_ring_radius();
  petsird::ReplicatedDetectorModule rep_module;

if (!stir::is_null_ptr(dynamic_cast<const stir::ProjDataInfoBlocksOnCylindrical * >(&stir_proj_data_info))) {
  //TODO
  //rep_module.ids.push_back(module_id++);
  //rep_module.transforms.push_back(transform);
} else if (!stir::is_null_ptr(dynamic_cast<const stir::ProjDataInfoCylindrical *>(&stir_proj_data_info))) {
    
    const std::array< int, 3> NUM_CRYSTALS_PER_MODULE{ stir_scanner->get_num_detector_layers(),
        stir_scanner->get_num_transaxial_crystals_per_block(),
        stir_scanner->get_num_axial_crystals_per_block()
        };
    
    const std::array<float, 3> crystal_dims{stir_scanner->get_average_depth_of_interaction(), 
        2*M_PI*radius / stir_scanner->get_num_detectors_per_ring(),
        stir_scanner->get_ring_spacing(),
        };
    auto box = get_crystal_template(crystal_dims);

    {
      rep_module.object = get_detector_module_tmpl(crystal_dims, NUM_CRYSTALS_PER_MODULE, radius);
      int module_id = 0;
      std::vector<float> angles;
      for (unsigned int i = 0; i < stir_scanner->get_num_transaxial_blocks(); ++i)
      {
        angles.push_back(static_cast<float>((2 * M_PI * i) / stir_scanner->get_num_transaxial_blocks()));
      }

      float MODULE_AXIS_SPACING = stir_scanner->get_num_rings() * stir_scanner->get_ring_spacing() / stir_scanner->get_num_axial_blocks(); 
      
      for (auto angle : angles)
        for (unsigned ax_mod = 0; ax_mod < stir_scanner->get_num_axial_blocks(); ++ax_mod)
        {
          petsird::RigidTransformation transform{ { { std::cos(angle), std::sin(angle), 0.F, 0.F },
                                                    { -std::sin(angle), std::cos(angle), 0.F, 0.F },
                                                    { 0.F, 0.F, 1.F, MODULE_AXIS_SPACING * ax_mod } } };
          rep_module.ids.push_back(module_id++);
          rep_module.transforms.push_back(transform);
        }
    }

} else {
    std::cout << "This should never happen! Abort" << std::endl;
}

  petsird::ScannerGeometry scanner_geometry;
  scanner_geometry.replicated_modules.push_back(rep_module);
  scanner_geometry.ids.push_back(0);

  typedef yardl::NDArray<float, 1> FArray1D;
  // TOF info (in mm)
  // Variables in capitals are to do to get from scanner.
  long unsigned int NUMBER_OF_TOF_BINS = stir_proj_data_info.get_num_tof_poss();
  long unsigned int NUMBER_OF_ENERGY_BINS = 1; //stir_scanner.get_num_energy_bins();
  float TOF_RESOLUTION = stir_scanner->get_timing_resolution();
  FArray1D::shape_type tof_bin_edges_shape = { NUMBER_OF_TOF_BINS + 1 };
  FArray1D tof_bin_edges(tof_bin_edges_shape);
  for (std::size_t i = 0; i < tof_bin_edges.size(); ++i)
    tof_bin_edges[i] = (i - NUMBER_OF_TOF_BINS / 2.F) / NUMBER_OF_TOF_BINS * 2 * radius;

  FArray1D::shape_type energy_bin_edges_shape = { NUMBER_OF_ENERGY_BINS + 1 };
  FArray1D energy_bin_edges(energy_bin_edges_shape);
  for (std::size_t i = 0; i < energy_bin_edges.size(); ++i)
    energy_bin_edges[i] = stir_exam_info.get_low_energy_thres() + 
    i * (stir_exam_info.get_high_energy_thres() - stir_exam_info.get_low_energy_thres()) / NUMBER_OF_ENERGY_BINS;
    
  petsird::ScannerInformation scanner_info;
  scanner_info.scanner_geometry = scanner_geometry;
  // TODO scanner_info.detectors = detectors;
  scanner_info.tof_bin_edges = tof_bin_edges;
  scanner_info.tof_resolution = TOF_RESOLUTION; // in mm
  scanner_info.energy_bin_edges = energy_bin_edges;
  scanner_info.energy_resolution_at_511 = stir_scanner->get_energy_resolution();    // as fraction of 511
  scanner_info.event_time_block_duration = 1.F; // ms
  return scanner_info;
}

petsird::Header
get_header()
{
  petsird::Subject subject;
  subject.id = "123456";
  petsird::Institution institution;
  institution.name = "ESTI Hackathon";
  institution.address = "Tampa, FL, USA";
  petsird::ExamInformation exam_info;
  exam_info.subject = subject;
  exam_info.institution = institution;
  petsird::Header header;
  header.exam = exam_info;
  return header;
}


STIRPETSIRDConvertor::STIRPETSIRDConvertor(const std::string& out_filename, const std::string& in_filename)
  : out_filename(out_filename), in_filename(in_filename)
{
  this->lm_data_ptr = stir::read_from_file<stir::ListModeData>(this->in_filename);
}


void
STIRPETSIRDConvertor::process_data()
{
  std::cout << "Converting STIR listmode data to petsird format...\n"
        << "\t- Input file: " << this->in_filename << "\n"
        << "\t- Output file: " << this->out_filename << "\n" << std::endl;

  using namespace stir;
  const auto& stir_proj_data_info_sptr = *lm_data_ptr->get_proj_data_info_sptr();
  const auto& stir_exam_info_sptr = *lm_data_ptr->get_exam_info_sptr();

  // Setup stir record, petsird time blocks and timing info
  auto record_sptr = lm_data_ptr->get_empty_record_sptr();
  auto& record = *record_sptr;

  petsird::EventTimeBlock  event_time_blk;
  petsird::ExternalSignalTimeBlock signal_time_blk;
  petsird::BedMovementTimeBlock bed_movement_time_blk;
  petsird::GantryMovementTimeBlock gantry_movement_time_blk;

  std::vector<petsird::CoincidenceEvent> prompts_this_blk;

  double current_time = 0.0;
  unsigned long num_events = 0;

  // Setup the petsird header info
  petsird::Header header_info = get_header();
  header_info.scanner = get_scanner_geometry(stir_proj_data_info_sptr, 
                                              stir_exam_info_sptr);

  petsird::binary::PETSIRDWriter writer(this->out_filename);
  writer.WriteHeader(header_info);
  std::cout << "I am here " << std::endl; 
  while (true)
    {

      if (lm_data_ptr->get_next_record(record) == Succeeded::no)
        {
          // no more events in file for some reason
          break; // get out of while loop
        }
      if (record.is_time())
        {
          current_time = record.time().get_time_in_millisecs();
          event_time_blk.start = current_time;
          event_time_blk.prompt_events = prompts_this_blk;
          writer.WriteTimeBlocks(event_time_blk);
          prompts_this_blk.clear();
        }
      if (record.is_event())
        {
          // assume it's a cylindrical scanner for now. will need to change later.
          auto& event = dynamic_cast<CListEventCylindricalScannerWithDiscreteDetectors const&>(record.event());

          DetectionPositionPair<> dp_pair;
          event.get_detection_position(dp_pair);

          petsird::CoincidenceEvent e;
          e.detector_ids[0]
              = dp_pair.pos1().tangential_coord() + dp_pair.pos1().axial_coord() * stir_proj_data_info_sptr.get_scanner_ptr()->get_num_detectors_per_ring();
          e.detector_ids[1]
              = dp_pair.pos2().tangential_coord() + dp_pair.pos2().axial_coord() * stir_proj_data_info_sptr.get_scanner_ptr()->get_num_detectors_per_ring();
          e.energy_indices[0] = 0;
          e.energy_indices[1] = 0;
          e.tof_idx = 0;
          prompts_this_blk.push_back(e);
          ++num_events;
        } // end of spatial event processing
    }     // end of while loop over all events
    writer.EndTimeBlocks();
    writer.Close();
    std::cout << "Done! Processed " << num_events << " events." << std::endl;
}

int
main(int argc, char* argv[])
{
  std::cout << "Converting STIR to PETSIRD" << std::endl; 
  if (argc != 3)
    {
      std::cout << "Converts list mode data from STIR to petsird format.\n"
                   "Usage: " << argv[0] << " <output_filename> <input_filename>\n";
      return 1;
    }

  STIRPETSIRDConvertor my_class(argv[1], argv[2]);
  my_class.process_data();

  return 0;
}
