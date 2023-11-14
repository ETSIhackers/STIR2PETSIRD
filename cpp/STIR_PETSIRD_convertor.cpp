#include <iostream>
#include <string>
#include <fstream>
#include <memory>

#include "STIR_PETSIRD_convertor.h"
#include "stir/listmode/ListModeData.h"
#include "stir/listmode/ListEvent.h"
#include "stir/listmode/ListRecord.h"
#include "stir/listmode/CListEventCylindricalScannerWithDiscreteDetectors.h"
#include "stir/LORCoordinates.h"
#include "stir/Succeeded.h"
#include "stir/IO/read_from_file.h"
#include "stir/error.h"
#include "hdf5/protocols.h"

#include "STIR_PETSIRD_convertor.h"

void
MyClass::process_data()
{
  using namespace stir;
  shared_ptr<ListModeData> lm_data_ptr(read_from_file<ListModeData>(this->in_filename));
  const auto& scanner = *lm_data_ptr->get_proj_data_info_sptr()->get_scanner_sptr();
  shared_ptr<ListRecord> record_sptr = lm_data_ptr->get_empty_record_sptr();
  ListRecord& record = *record_sptr;

  double current_time = 0.0;

  prd::hdf5::PrdExperimentWriter writer(this->out_filename);
  prd::TimeBlock time_block;
  std::vector<prd::CoincidenceEvent> prompts_this_block;

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
          /*
            auto lor = record.event().get_LOR();
            lor.p1().x();
            lor.p1().y();
            lor.p2();
          */
          DetectionPositionPair<> dp_pair;
          event.get_detection_position(dp_pair);

          prd::CoincidenceEvent e;
          e.detector_1_id
              = dp_pair.pos1().tangential_coord() + dp_pair.pos1().axial_coord() * scanner.get_num_detectors_per_ring();
          e.detector_2_id
              = dp_pair.pos2().tangential_coord() + dp_pair.pos2().axial_coord() * scanner.get_num_detectors_per_ring();
          e.energy_1_idx = 0;
          e.energy_2_idx = 0;
          e.tof_idx = 0;
          prompts_this_block.push_back(e);
        } // end of spatial event processing
    }     // end of while loop over all events
}

int
main(int argc, char* argv[])
{
  if (argc != 3)
    {
      std::cerr << "Usage: " << argv[0] << " <output_filename> <input_filename>\n";
      return 1;
    }

  MyClass my_class(argv[1], argv[2]);
  my_class.process_data();

  return 0;
}
