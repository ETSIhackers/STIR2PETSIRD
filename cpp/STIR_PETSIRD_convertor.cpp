#include <iostream>
#include <string>
#include <fstream>
#include <memory>
#include <filesystem>

#include "stir/listmode/ListModeData.h"
#include "stir/listmode/ListEvent.h"
#include "stir/listmode/ListRecord.h"
#include "stir/listmode/CListEventCylindricalScannerWithDiscreteDetectors.h"
#include "stir/LORCoordinates.h"
#include "stir/Succeeded.h"
#include "stir/IO/read_from_file.h"
#include "stir/error.h"
#include "binary/protocols.h"

#include "STIR_PETSIRD_convertor.h"

// Convert from STIR scanner to PRD scanner info (for now, just cylindrical non-TOF scanners)
prd::ScannerInformation
get_scanner_info(const stir::Scanner& stir_scanner)
{
  float radius = stir_scanner.get_inner_ring_radius();
  std::vector<float> angles;
  for (int i = 0; i < stir_scanner.get_num_detectors_per_ring(); ++i)
      angles.push_back(static_cast<float>(2 * M_PI * i / stir_scanner.get_num_detectors_per_ring()));

  std::vector<prd::Detector> detectors;
  int detector_id = 0;
  int num_rings = stir_scanner.get_num_rings();
  for (int ring = 0; ring < num_rings; ++ring)
  {
    for (auto angle : angles)
    {
      // Create a new detector
      prd::Detector d;
      d.x = radius * std::sin(angle);
      d.y = radius * std::cos(angle);
      d.z = (float)ring;
      d.id = detector_id++;
      detectors.push_back(d);
    }
  }

  typedef yardl::NDArray<float, 1> FArray1D;
  // TOF info (in mm)
  // Variables in capitals are to do to get from scanner.
  long unsigned int NUMBER_OF_TOF_BINS = 1;
  long unsigned int NUMBER_OF_ENERGY_BINS = 1;
  float TOF_RESOLUTION = -1.0;
  FArray1D::shape_type tof_bin_edges_shape = { NUMBER_OF_TOF_BINS + 1 };
  FArray1D tof_bin_edges(tof_bin_edges_shape);
  for (std::size_t i = 0; i < tof_bin_edges.size(); ++i)
    tof_bin_edges[i] = (i - NUMBER_OF_TOF_BINS / 2.F) / NUMBER_OF_TOF_BINS * 2 * radius;
  FArray1D::shape_type energy_bin_edges_shape = { NUMBER_OF_ENERGY_BINS + 1 };
  FArray1D energy_bin_edges(energy_bin_edges_shape);
  for (std::size_t i = 0; i < energy_bin_edges.size(); ++i)
    energy_bin_edges[i] = 430.F + i * (650.F - 430.F) / NUMBER_OF_ENERGY_BINS;
  prd::ScannerInformation scanner_info;
  scanner_info.detectors = detectors;
  scanner_info.tof_bin_edges = tof_bin_edges;
  scanner_info.tof_resolution = TOF_RESOLUTION; // in mm
  scanner_info.energy_bin_edges = energy_bin_edges;
  scanner_info.energy_resolution_at_511 = stir_scanner.get_energy_resolution();    // as fraction of 511
  scanner_info.listmode_time_block_duration = 1.F; // ms
  return scanner_info;
}

prd::Header
get_header()
{
  prd::Subject subject;
  subject.id = "123456";
  prd::Institution institution;
  institution.name = "ESTI Hackathon";
  institution.address = "Vancouver, Canada";
  prd::ExamInformation exam_info;
  exam_info.subject = subject;
  exam_info.institution = institution;
  prd::Header header;
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
  std::cout << "Converting STIR listmode data to PRD format...\n"
        << "\t- Input file: " << this->in_filename << "\n"
        << "\t- Output file: " << this->out_filename << "\n" << std::endl;

  using namespace stir;
  const auto& stir_scanner = *lm_data_ptr->get_proj_data_info_sptr()->get_scanner_ptr();

  // Setup stir record, prd time blocks and timing info
  auto record_sptr = lm_data_ptr->get_empty_record_sptr();
  auto& record = *record_sptr;
  prd::TimeBlock time_block;
  std::vector<prd::CoincidenceEvent> prompts_this_block;
  double current_time = 0.0;

  // Setup the prd header info
  prd::Header header_info = get_header();
  header_info.scanner = get_scanner_info(stir_scanner);

  // Remove previous out_filename if it exists
  if (std::filesystem::exists(this->out_filename))
      std::filesystem::remove(this->out_filename);

  prd::binary::PrdExperimentWriter writer(this->out_filename);
  writer.WriteHeader(header_info);

  while (true)
    {
      if (lm_data_ptr->get_next_record(record) == Succeeded::no)
        {
          // no more events in file for some reason
          break; // get out of while loop
        }
      if (record.is_time())
        {
          current_time = record.time().get_time_in_secs();
          time_block.id = current_time;
          time_block.prompt_events = prompts_this_block;
          writer.WriteTimeBlocks(time_block);
          prompts_this_block.clear();
        }
      if (record.is_event())
        {
          // assume it's a cylindrical scanner for now. will need to change later.
          auto& event = dynamic_cast<CListEventCylindricalScannerWithDiscreteDetectors const&>(record.event());

          DetectionPositionPair<> dp_pair;
          event.get_detection_position(dp_pair);

          prd::CoincidenceEvent e;
          e.detector_1_id
              = dp_pair.pos1().tangential_coord() + dp_pair.pos1().axial_coord() * stir_scanner.get_num_detectors_per_ring();
          e.detector_2_id
              = dp_pair.pos2().tangential_coord() + dp_pair.pos2().axial_coord() * stir_scanner.get_num_detectors_per_ring();
          e.energy_1_idx = 0;
          e.energy_2_idx = 0;
          e.tof_idx = 0;
          prompts_this_block.push_back(e);
        } // end of spatial event processing
    }     // end of while loop over all events
    writer.EndTimeBlocks();
    writer.Close();
}

int
main(int argc, char* argv[])
{
  if (argc != 3)
    {
      std::cout << "Converts list mode data from STIR to PRD format.\n"
                   "Usage: " << argv[0] << " <output_filename> <input_filename>\n";
      return 1;
    }

  STIRPETSIRDConvertor my_class(argv[1], argv[2]);
  my_class.process_data();

  return 0;
}
