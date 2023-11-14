
#include <iostream>
#include "stir/listmode/ListModeData.h"
#include "stir/listmode/ListEvent.h"
#include "stir/listmode/ListRecord.h"
#include "stir/LORCoordinates.h"
#include "stir/Succeeded.h"
#include "stir/IO/read_from_file.h"
#include "/workspaces/STIR2PETSIRD/PETSIRD/cpp/generated/hdf5/protocols.h"

#include <string>
#include <fstream>
using namespace stir;

class MyClass {
public:
    MyClass(const std::string& out_filename, std::string in_filename)  {
        this->out_filename = out_filename;
        this->in_filename = in_filename;
    }



    void process_data()
    {
        shared_ptr<ListModeData> lm_data_ptr(read_from_file<ListModeData>(this->in_filename));
        shared_ptr <ListRecord> record_sptr = lm_data_ptr->get_empty_record_sptr();
        ListRecord& record = *record_sptr;

        double current_time = 0.0;

        prd::hdf5::PrdExperimentWriter writer(this->out_filename);
        prd::TimeBlock time_block;
        std::vector<prd::CoincidenceEvent> prompts_this_block;


        while(true) {
	    if (lm_data_ptr->get_next_record(record) == Succeeded::no) {
		    // no more events in file for some reason
		    break; //get out of while loop
		}
        if (record.is_time()) {
		    current_time = record.time().get_time_in_secs();
            time_block.id = current_time;
            time_block.prompt_events = prompts_this_block;
            writer.WriteTimeBlocks(time_block);
            prompts_this_block.clear();
        }
		if (record.is_event()) {
            auto lor = record.event().get_LOR();
            lor.p1().x();
            lor.p1().y();
            lor.p2();

            // prd::CoincidenceEvent e;
            // e.detector_1_id = detectors.first;
            // e.detector_2_id = detectors.second;
            // e.energy_1_idx = get_random_energy_value();
            // e.energy_2_idx = get_random_energy_value();
            // e.tof_idx = get_random_tof_value();
            // prompts_this_block.push_back(e);

        } // end of spatial event processing
    } // end of while loop over all events
}



private:
    std::string out_filename;
    std::string in_filename;
};





int
main(int argc, char* argv[])
{
  // Check if the user has provided a file
  if (argc < 3)
    {
      std::cerr << "Please provide a filename to write to" << std::endl;
      return 1;
    }

    MyClass my_class(argv[1], argv[2]);
  return 0;
}
