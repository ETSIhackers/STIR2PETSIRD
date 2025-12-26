/*!
  \file
  \brief Implementation of class STIRPETSIRDConvertor

  \author Nikos Efthimiou
  \author Eve Lennie
  \author Robert Twyman Skelly
  \author Kris Thielemans

*/
/*
    Copyright (C) 2023, Prescient Imaging
    Copyright (C) 2023, University of Sheffield
    Copyright (C) 2024, MGH
    Copyright (C) 2024, University College of London
    This file is part of STIR.

    SPDX-License-Identifier: Apache-2.0
*/

#include <iostream>
#include <string>
#include <fstream>
#include <memory>

#include "STIR_PETSIRD_convertor.h"
#include "convertor_helpers.h"

#include "stir/KeyParser.h"
#include "stir/IO/read_from_file.h"


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
  const auto stir_proj_data_info_sptr = lm_data_ptr->get_proj_data_info_sptr();
  const auto stir_exam_info_sptr = lm_data_ptr->get_exam_info_sptr();
  const auto stir_scanner = stir_proj_data_info_sptr->get_scanner_ptr(); 

  std::cout << "Start check" << std::endl;
  check_id_conversion(stir_proj_data_info_sptr->get_scanner_ptr());
  std::cout << "Done" << std::endl;

  // Setup stir record, petsird time blocks and timing info
  auto record_sptr = lm_data_ptr->get_empty_record_sptr();
  auto& record = *record_sptr;

  petsird::EventTimeBlock  event_time_blk;
  event_time_blk.time_interval.start = 0;
  event_time_blk.prompt_events.resize(1);
  event_time_blk.prompt_events[0].resize(1);
  event_time_blk.delayed_events = std::vector<std::vector<petsird::ListOfCoincidenceEvents>>();
  event_time_blk.delayed_events.resize(1);
  event_time_blk.delayed_events[0].resize(1);


  petsird::ExternalSignalTimeBlock signal_time_blk;
  petsird::BedMovementTimeBlock bed_movement_time_blk;
  petsird::GantryMovementTimeBlock gantry_movement_time_blk;

  std::vector<petsird::CoincidenceEvent> prompts_this_blk;
  std::vector<petsird::CoincidenceEvent> delayeds_this_blk;

  double current_time = 0.0;
  unsigned long num_events = 0;

  // Setup the petsird header info
  petsird::Header header_info = get_dummy_header();
  auto& scanner_info = header_info.scanner;
  const auto num_types_of_modules = 1;
  // Pre-allocate various structures to have the correct size for num_types_of_modules
  // (We will still have to set descent values into each of these.)
  petsird_helpers::create::initialize_scanner_information_dimensions(scanner_info, num_types_of_modules,
                                                                     /* allocate_detection_bin_efficiencies = */ false,
                                                                     /* allocate_module_pair_efficiencies = */ this->normalisation_sptr != nullptr);
  set_scanner_geometry(scanner_info, *stir_proj_data_info_sptr,
                       *stir_exam_info_sptr);

  if (this->normalisation_sptr)
    {
      if (!this->normalisation_sptr->set_up(stir_exam_info_sptr, stir_proj_data_info_sptr).succeeded())
        stir::error("Error setting up norm");

      get_detection_efficiencies(header_info.scanner, *stir_proj_data_info_sptr, *stir_exam_info_sptr, *this->normalisation_sptr);
    }
  petsird::binary::PETSIRDWriter writer(this->out_filename);
  writer.WriteHeader(header_info);

  // unfortunately we need to check cylindrical vs generic ATM
  auto stir_proj_data_info_generic_noarc_sptr = std::dynamic_pointer_cast<stir::ProjDataInfoGenericNoArcCorr const>(stir_proj_data_info_sptr);
  auto stir_proj_data_info_cylindrical_noarc_sptr = std::dynamic_pointer_cast<stir::ProjDataInfoCylindricalNoArcCorr const>(stir_proj_data_info_sptr);
  if (!stir_proj_data_info_generic_noarc_sptr && !stir_proj_data_info_cylindrical_noarc_sptr)
    stir::error("STIR data has to be not arccorrected due to code limitations");

  long num_events_to_process = -1; // set to -1 to process all

  while (num_events_to_process)
    {

      if (lm_data_ptr->get_next_record(record) == Succeeded::no)
        {
          // no more events in file for some reason
          break; // get out of while loop
        }
      if (record.is_time())
        {
          current_time = record.time().get_time_in_millisecs();
          event_time_blk.time_interval.stop = current_time;
          event_time_blk.prompt_events[0][0] = prompts_this_blk;
          event_time_blk.delayed_events[0][0] = delayeds_this_blk;
          writer.WriteTimeBlocks(event_time_blk);
          event_time_blk.time_interval.start = current_time;
          prompts_this_blk.clear();
          delayeds_this_blk.clear();
        }
      if (record.is_event())
        {
          if (num_events_to_process > 0) // we are using -1 when processing all
            num_events_to_process--;

          DetectionPositionPair<> dp_pair;
          if (stir_proj_data_info_cylindrical_noarc_sptr)
            dp_pair = get_det_pos_pair_help(record.event(), *stir_proj_data_info_cylindrical_noarc_sptr);
          else
            dp_pair = get_det_pos_pair_help(record.event(), *stir_proj_data_info_generic_noarc_sptr);

          petsird::CoincidenceEvent e;
          e.detection_bins[0] = get_PETSIRD_id_from_stir_det_pos(dp_pair.pos1(), stir_scanner);
          e.detection_bins[1] = get_PETSIRD_id_from_stir_det_pos(dp_pair.pos2(), stir_scanner);
          e.tof_idx = dp_pair.timing_pos() - stir_proj_data_info_sptr->get_min_tof_pos_num();
#if 0 // redundant test, so commented out
          {
            const petsird::TypeOfModule type_of_module{0};
            const auto expanded_detection_bin0
              = petsird_helpers::expand_detection_bin(scanner_info, type_of_module, e.detection_bins[0]);
            const auto expanded_detection_bin1
              = petsird_helpers::expand_detection_bin(scanner_info, type_of_module, e.detection_bins[1]);
            const auto box_shape0 = petsird_helpers::geometry::get_detecting_box(scanner_info, type_of_module, expanded_detection_bin0);
            const auto mean_pos0 = mean_position(box_shape0);
            const auto box_shape1 = petsird_helpers::geometry::get_detecting_box(scanner_info, type_of_module, expanded_detection_bin1);
            const auto mean_pos1 = mean_position(box_shape1);
            const auto p1 = make_coordinate(mean_pos0.c[2], -mean_pos0.c[0], -mean_pos0.c[1]);
            const auto p2 = make_coordinate(mean_pos1.c[2], -mean_pos1.c[0], -mean_pos1.c[1]);
            const auto LOR = record.event().get_LOR();
            bool swap = dp_pair.timing_pos() < 0; // get_LOR() will swap p1,p2 if timing_pos changes sign
            auto diff0 = (swap ? LOR.p1() : LOR.p2()) - p1;
            auto diff1 = (swap ? LOR.p2() : LOR.p1()) - p2;
            if (dp_pair.timing_pos() == 0 && norm(diff0) > 100) // if timing_pos == 0, it could be either case.
              {
                swap = !swap;
                diff0 = (swap ? LOR.p1() : LOR.p2()) - p1;
                diff1 = (swap ? LOR.p2() : LOR.p1()) - p2;
              }

            std::cout << dp_pair.timing_pos() << LOR.p1() << LOR.p2() << p1 << p2<< diff0 << diff1 << "\n";
          }
#endif
          if (record.event().is_prompt())
            {
              prompts_this_blk.push_back(e);
              ++num_events;
            }
          else
            {
              delayeds_this_blk.push_back(e);
            }
        } // end of spatial event processing
        if (num_events%100000 == 0)
        {
          std::cout << num_events << std::endl; 
        }
        // if (num_events == 250)
        //   break; 
    }     // end of while loop over all events
    writer.EndTimeBlocks();
    writer.Close();
    std::cout << "Done! Processed " << num_events << " events." << std::endl;
}


USING_NAMESPACE_STIR
int
main(int argc, char* argv[])
{
  std::cout << "Converting STIR to PETSIRD" << std::endl; 
  if (argc < 3 || argc > 4)
    {
      std::cout << "Converts list mode data from STIR to petsird format.\n"
                   "Usage: " << argv[0] << " <output_filename> <input_filename> [norm.par]\n";
      return EXIT_FAILURE;
    }

  STIRPETSIRDConvertor my_class(argv[1], argv[2]);
  if (argc == 4)
    {
      std::shared_ptr<BinNormalisation> normalisation_sptr;
      KeyParser parser;
      parser.add_start_key("Bin normalisation parameters");
      parser.add_parsing_key("type", &normalisation_sptr);
      parser.add_stop_key("END");
      parser.parse(argv[3]);
      my_class.set_normalisation_sptr(normalisation_sptr);
      
    }
  my_class.process_data();

  return EXIT_SUCCESS;
}
