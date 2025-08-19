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

#include "stir/listmode/ListModeData.h"
#include "stir/listmode/ListEvent.h"
#include "stir/listmode/ListRecord.h"
#include "stir/listmode/CListEventCylindricalScannerWithDiscreteDetectors.h"
#include "stir/ProjDataInfoCylindricalNoArcCorr.h"
#include "stir/ProjDataInfoBlocksOnCylindricalNoArcCorr.h"
#include "stir/recon_buildblock/BinNormalisation.h"
#include "stir/Bin.h"
#include "stir/DetectionPositionPair.h"
#include "stir/LORCoordinates.h"
#include "stir/Succeeded.h"
#include "stir/KeyParser.h"
#include "stir/IO/read_from_file.h"
#include "stir/error.h"
#include "binary/protocols.h"
#include "petsird_helpers.h"
#include "petsird_helpers/create.h"
#include "petsird_helpers/geometry.h"
#include "stir/ProjDataInfoGenericNoArcCorr.h"
#include "stir/ProjDataInfoBlocksOnCylindrical.h"
#include "stir/ProjDataInfoCylindrical.h"
#include "stir/stream.h"
#include "stir/DetectionPosition.h"

#include "STIR_PETSIRD_convertor.h"

constexpr float speed_of_light_mm_per_ps = 0.299792458F;

static petsird::Coordinate mean_position(const petsird::BoxShape& box_shape)
{
  petsird::Coordinate mean;
  mean.c = {0, 0, 0};
  for (auto& corner : box_shape.corners)
    {
      mean.c += corner.c;
    }
  mean.c /= box_shape.corners.size();
  return mean;
}

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

template <class T>
inline
std::array<T, 3> get_indices_from_id(T id, const std::array<T, 3>& sizes)
{
  const auto N0 = sizes[0];
  const auto N1 = sizes[1];
  const auto N2 = sizes[2];
  assert(id < N0 * N1 * N2);
  std::array<T, 3> inds;
  inds[2] = id % N2;
  id = id / N2;
  inds[1] = id %  N1;
  id = id / N1;
  inds[0] = id;
  return inds;
}

template <class T>
inline
T get_id_from_indices(const std::array<T, 3>& inds, const std::array<T, 3>& sizes)
{
  const auto N0 = sizes[0];
  const auto N1 = sizes[1];
  const auto N2 = sizes[2];
  assert(inds[0] < N0);
  assert(inds[1] < N1);
  assert(inds[2] < N2);
  const auto id = inds[2] + N2 * (inds[1] + N1 * inds[0]);
  // assert(inds == get_indices_from_id(id, sizes));
  return id;
}

inline
stir::DetectionPosition<> get_stir_det_pos_from_PETSIRD_id(const petsird_helpers::ExpandedDetectionBin& exp_det_bin, const stir::Scanner* const stir_scanner)
{
  // const auto num_det_els_in_module = stir_scanner->get_num_axial_crystals_per_block() * stir_scanner->get_num_transaxial_crystals_per_block();
  // const auto num_modules = stir_scanner->get_num_transaxial_blocks() * stir_scanner->get_num_axial_blocks();
  // const auto NUM_MODULES_ALONG_RING = stir_scanner->get_num_transaxial_blocks();
  const auto NUM_MODULES_ALONG_AXIS = stir_scanner->get_num_axial_blocks();
  const std::array<uint32_t, 3> NUM_CRYSTALS_PER_MODULE{ static_cast<unsigned>(stir_scanner->get_num_detector_layers()),
    static_cast<unsigned>(stir_scanner->get_num_transaxial_crystals_per_block()),
    static_cast<unsigned>(stir_scanner->get_num_axial_crystals_per_block())
      };
  const auto ax_mod = exp_det_bin.module_index % NUM_MODULES_ALONG_AXIS;
  const auto tang_mod = exp_det_bin.module_index / NUM_MODULES_ALONG_AXIS;

  const auto inds = get_indices_from_id(exp_det_bin.element_index, NUM_CRYSTALS_PER_MODULE);
  const stir::DetectionPosition<> pos(inds[1] + tang_mod * NUM_CRYSTALS_PER_MODULE[1],
                                      inds[2] + ax_mod * NUM_CRYSTALS_PER_MODULE[2],
                                      inds[0]);
  return pos;
}

inline
petsird::DetectionBin get_PETSIRD_id_from_stir_det_pos(const stir::DetectionPosition<>& det_pos, const stir::Scanner* const stir_scanner)
{
  const auto num_det_els_in_module = static_cast<unsigned>(stir_scanner->get_num_axial_crystals_per_block() * stir_scanner->get_num_transaxial_crystals_per_block());
  // const auto num_modules = stir_scanner->get_num_transaxial_blocks() * stir_scanner->get_num_axial_blocks();
  //const auto NUM_MODULES_ALONG_RING = stir_scanner->get_num_transaxial_blocks();
  const auto NUM_MODULES_ALONG_AXIS = static_cast<unsigned>(stir_scanner->get_num_axial_blocks());
  const std::array< std::size_t, 3> NUM_CRYSTALS_PER_MODULE{ static_cast<unsigned>(stir_scanner->get_num_detector_layers()),
    static_cast<unsigned>(stir_scanner->get_num_transaxial_crystals_per_block()),
    static_cast<unsigned>(stir_scanner->get_num_axial_crystals_per_block())
      };

  std::array<std::size_t, 3> inds;
  inds[1] = det_pos.tangential_coord() % NUM_CRYSTALS_PER_MODULE[1];
  const auto tang_mod = det_pos.tangential_coord() / NUM_CRYSTALS_PER_MODULE[1];
  inds[2] = det_pos.axial_coord() % NUM_CRYSTALS_PER_MODULE[2];
  const auto ax_mod = det_pos.axial_coord() / NUM_CRYSTALS_PER_MODULE[2];
  inds[0] = det_pos.radial_coord();
  const auto mod = ax_mod + tang_mod * NUM_MODULES_ALONG_AXIS;
  const auto det_el = get_id_from_indices(inds, NUM_CRYSTALS_PER_MODULE);
  // get_module_and_element uses { det / num_el_per_module, det % num_el_per_module }
  return det_el + mod * num_det_els_in_module;
}

void
check_id_conversion(const stir::Scanner* const stir_scanner)
{
  const auto num_det_els_in_module = stir_scanner->get_num_axial_crystals_per_block() * stir_scanner->get_num_transaxial_crystals_per_block();
  stir::DetectionPosition<> det_pos;
  for (det_pos.radial_coord() = 0; det_pos.radial_coord() < static_cast<unsigned>(stir_scanner->get_num_detector_layers()); ++det_pos.radial_coord())
    for (det_pos.axial_coord() = 0; det_pos.axial_coord() < static_cast<unsigned>(stir_scanner->get_num_rings()); ++det_pos.axial_coord())
      for (det_pos.tangential_coord() = 0; det_pos.tangential_coord() < static_cast<unsigned>(stir_scanner->get_num_detectors_per_ring()); ++det_pos.tangential_coord())
        {
          const auto id{get_PETSIRD_id_from_stir_det_pos(det_pos, stir_scanner)};
          const petsird_helpers::ExpandedDetectionBin exp_det_bin{id / num_det_els_in_module, id % num_det_els_in_module, 0};
          const auto new_det_pos{get_stir_det_pos_from_PETSIRD_id(exp_det_bin, stir_scanner)};
          if (det_pos != new_det_pos)
            stir::error("Error round-trip");
        }
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
                                                      { 0.0, 1.0, 0.0, -(rep1 - N1 / 2) * crystal_length[1] }, // somewhat surprising minus sign due to how we rotate (KT thinks)
                                                      { 0.0, 0.0, 1.0, (rep2 - N2 / 2) * crystal_length[2] } } };
            rep_volume.transforms.push_back(transform);
          }
  }

  petsird::DetectorModule detector_module;
  detector_module.detecting_elements = rep_volume;

  return detector_module;
}

// Convert from STIR scanner to petsird scanner info (for now, just cylindrical non-TOF scanners)
void
set_scanner_geometry(petsird::ScannerInformation& scanner_info,
                     const stir::ProjDataInfo& stir_proj_data_info, 
                     const stir::ExamInfo& stir_exam_info)
{
  const stir::Scanner* stir_scanner = stir_proj_data_info.get_scanner_ptr(); 
  const float radius = stir_scanner->get_inner_ring_radius();
  petsird::ReplicatedDetectorModule rep_module;

if (!stir::is_null_ptr(dynamic_cast<const stir::ProjDataInfoBlocksOnCylindrical * >(&stir_proj_data_info))) {

// NE: I would like to use some of the stuff in the norm branch, but first merge. 
      const std::array< int, 3> NUM_CRYSTALS_PER_BLOCK{ stir_scanner->get_num_detector_layers(),
        stir_scanner->get_num_transaxial_crystals_per_block(),
        stir_scanner->get_num_axial_crystals_per_block()
        };
    
    const std::array<float, 3> crystal_dims{
        stir_scanner->get_average_depth_of_interaction(),
        stir_scanner->get_transaxial_crystal_spacing(),
        stir_scanner->get_axial_crystal_spacing(),
        };

    {
      rep_module.object = get_detector_module_tmpl(crystal_dims, NUM_CRYSTALS_PER_BLOCK, radius);
      std::vector<float> angles;
      for (int i = 0; i < stir_scanner->get_num_transaxial_blocks(); ++i)
      {
        angles.push_back(static_cast<float>((2 * M_PI * i) / stir_scanner->get_num_transaxial_blocks()));
      }

      float MODULE_AXIS_SPACING = stir_scanner->get_num_rings() * stir_scanner->get_ring_spacing() / stir_scanner->get_num_axial_blocks(); 
      
      for (auto angle : angles)
        for (int ax_mod = 0; ax_mod < stir_scanner->get_num_axial_blocks(); ++ax_mod)
        {
          petsird::RigidTransformation transform{ { { std::cos(angle), std::sin(angle), 0.F, 0.F },
                                                    { -std::sin(angle), std::cos(angle), 0.F, 0.F },
                                                    { 0.F, 0.F, 1.F, MODULE_AXIS_SPACING * ax_mod } } };
          rep_module.transforms.push_back(transform);
        }
    }

} else if (!stir::is_null_ptr(dynamic_cast<const stir::ProjDataInfoCylindrical *>(&stir_proj_data_info))) {
     
    const std::array< int, 3> NUM_CRYSTALS_PER_MODULE{ stir_scanner->get_num_detector_layers(),
        stir_scanner->get_num_transaxial_crystals_per_block(),
        stir_scanner->get_num_axial_crystals_per_block()
        };
    
    const std::array<float, 3> crystal_dims{stir_scanner->get_average_depth_of_interaction(), 
      static_cast<float>(2*M_PI*radius / stir_scanner->get_num_detectors_per_ring()),
        stir_scanner->get_ring_spacing(),
        };

    {
      rep_module.object = get_detector_module_tmpl(crystal_dims, NUM_CRYSTALS_PER_MODULE, radius);
      std::vector<float> angles;
      for (int i = 0; i < stir_scanner->get_num_transaxial_blocks(); ++i)
      {
        angles.push_back(static_cast<float>((2 * M_PI * i) / stir_scanner->get_num_transaxial_blocks()));
      }

      float MODULE_AXIS_SPACING = stir_scanner->get_num_rings() * stir_scanner->get_ring_spacing() / stir_scanner->get_num_axial_blocks(); 
      
      for (auto angle : angles)
        for (int ax_mod = 0; ax_mod < stir_scanner->get_num_axial_blocks(); ++ax_mod)
        {
          petsird::RigidTransformation transform{ { { std::cos(angle), std::sin(angle), 0.F, 0.F },
                                                    { -std::sin(angle), std::cos(angle), 0.F, 0.F },
                                                    { 0.F, 0.F, 1.F, MODULE_AXIS_SPACING * ax_mod } } };
          rep_module.transforms.push_back(transform);
        }
    }

} else {
    std::cout << "This should never happen! Abort" << std::endl;
}

 auto& scanner_geometry = scanner_info.scanner_geometry;
 scanner_geometry.replicated_modules.push_back(rep_module);

  // TOF and energy information
  {
    auto& all_tof_bin_edges = scanner_info.tof_bin_edges;
    auto& all_tof_resolutions = scanner_info.tof_resolution;
    auto& all_event_energy_bin_edges = scanner_info.event_energy_bin_edges;
    auto& all_event_energy_resolutions = scanner_info.energy_resolution_at_511;

    // only 1 type of module in the current scanner
    const petsird::TypeOfModule type_of_module{ 0 };

    typedef yardl::NDArray<float, 1> FArray1D;
    // TOF info (in mm)
    // Variables in capitals are to do to get from scanner.
    long unsigned int NUMBER_OF_TOF_BINS = stir_proj_data_info.get_num_tof_poss();
    long unsigned int NUMBER_OF_EVENT_ENERGY_BINS = 1; //stir_scanner.get_num_energy_bins();
    const float TOF_RESOLUTION = stir_scanner->get_timing_resolution() * speed_of_light_mm_per_ps / 2;
    const float TOF_bin_width_mm = stir_scanner->get_size_of_timing_pos() * speed_of_light_mm_per_ps / 2;
    FArray1D tof_bin_edges_arr;
    yardl::resize(tof_bin_edges_arr, { NUMBER_OF_TOF_BINS + 1 });
    for (std::size_t i = 0; i < tof_bin_edges_arr.size(); ++i)
      tof_bin_edges_arr[i] = (i - NUMBER_OF_TOF_BINS / 2.F) * TOF_bin_width_mm;
    const petsird::BinEdges tof_bin_edges{ tof_bin_edges_arr };
    all_tof_bin_edges[type_of_module][type_of_module] = tof_bin_edges;

    all_tof_resolutions[type_of_module][type_of_module] = TOF_RESOLUTION;

    FArray1D event_energy_bin_edges_arr;
    yardl::resize(event_energy_bin_edges_arr, { NUMBER_OF_EVENT_ENERGY_BINS + 1 });
    const auto energy_LLD = stir_exam_info.get_low_energy_thres();
    const auto energy_ULD = stir_exam_info.get_high_energy_thres();
    for (std::size_t i = 0; i < event_energy_bin_edges_arr.size(); ++i)
      event_energy_bin_edges_arr[i] = energy_LLD + i * (energy_ULD - energy_LLD) / NUMBER_OF_EVENT_ENERGY_BINS;
    petsird::BinEdges event_energy_bin_edges{ event_energy_bin_edges_arr };
    all_event_energy_bin_edges[type_of_module] = event_energy_bin_edges;
    all_event_energy_resolutions[type_of_module] = stir_scanner->get_energy_resolution();    // as fraction of 511 (e.g. 0.11F)
  }
  // test
  {
    stir::CartesianCoordinate3D<float> coord_0, coord_1;
    auto& pdi = dynamic_cast<const stir::ProjDataInfoCylindricalNoArcCorr &>(stir_proj_data_info);
    for (unsigned int r=0; r< static_cast<unsigned>(stir_scanner->get_num_rings()); ++r)
      for (unsigned int d=0; d < static_cast<unsigned>(stir_scanner->get_num_detectors_per_ring()); ++d)
        {
          pdi.find_cartesian_coordinates_given_scanner_coordinates(coord_0, coord_1, r, 0, d, 0, 1);
          stir::DetectionPosition<> det_pos{d, r, 0};
          const auto det_bin = get_PETSIRD_id_from_stir_det_pos(det_pos, stir_scanner);
          const petsird::TypeOfModule type_of_module{0};
          const auto expanded_detection_bin
            = petsird_helpers::expand_detection_bin(scanner_info, type_of_module, det_bin);
          const auto box_shape = petsird_helpers::geometry::get_detecting_box(scanner_info, type_of_module, expanded_detection_bin);
          const auto mean_pos = mean_position(box_shape);
          const auto p0 = stir::make_coordinate(mean_pos.c[2], -mean_pos.c[0], -mean_pos.c[1]);
          const auto diff = coord_0 - p0;
          std::cout << det_pos << coord_0 << p0 << diff << "\n";
        }
  }
  // TODO scanner_info.coincidence_policy = petsird::CoincidencePolicy::kRejectMultiples;
  scanner_info.delayed_coincidences_are_stored = true;
  scanner_info.triple_events_are_stored = false;
}

template <class ProjDataInfoT>
void
get_detection_efficiencies_help(petsird::ScannerInformation& scanner,
                                const ProjDataInfoT& stir_proj_data_info,
                                const stir::ExamInfo& stir_exam_info,
                                const stir::BinNormalisation& norm)
{
  const stir::Scanner* stir_scanner = stir_proj_data_info.get_scanner_ptr(); 

  const auto num_modules = static_cast<unsigned>(stir_scanner->get_num_transaxial_blocks() * stir_scanner->get_num_axial_blocks());
  const auto NUM_MODULES_ALONG_RING = static_cast<unsigned>(stir_scanner->get_num_transaxial_blocks());
  const auto NUM_MODULES_ALONG_AXIS = static_cast<unsigned>(stir_scanner->get_num_axial_blocks());
  // TODO could do axial_fan_size based on max ring diff
  auto fan_size = std::ceil(static_cast<float>(stir_proj_data_info.get_num_tangential_poss()) / stir_scanner->get_num_transaxial_crystals_per_block());
  std::cerr << "Module fan_size along the ring : " << fan_size << std::endl;

  constexpr petsird::TypeOfModule type_of_module{0};
  const auto NZ = NUM_MODULES_ALONG_AXIS;
  auto& module_pair_SGID_LUT = (*scanner.detection_efficiencies.module_pair_sgidlut)[type_of_module][type_of_module];
  module_pair_SGID_LUT = yardl::NDArray<int, 2>({ num_modules, num_modules });

  int num_SGIDs = 0;
  {
    for (unsigned int mod1 = 0; mod1 < num_modules; ++mod1)
      {
        for (unsigned int mod2 = 0; mod2 < num_modules; ++mod2)
          {
            // const auto z1 = mod1 % NZ;
            const auto a1 = mod1 / NZ;
            // const auto z2 = mod2 % NZ;
            const auto a2 = mod2 / NZ;
            if (std::abs(2 * std::abs(int(a1) - int(a2)) - int(NUM_MODULES_ALONG_RING)) > fan_size)
              {
                module_pair_SGID_LUT(mod1, mod2) = -1;
              }
            else
              {
                module_pair_SGID_LUT(mod1, mod2) = num_SGIDs++;
              }
          }
      }
  }
  auto& module_pair_efficiencies_vector
      = (*scanner.detection_efficiencies.module_pair_efficiencies_vectors)[type_of_module][type_of_module];
  module_pair_efficiencies_vector.reserve(num_SGIDs);

  const auto num_det_els_in_module = static_cast<unsigned>(stir_scanner->get_num_axial_crystals_per_block() * stir_scanner->get_num_transaxial_crystals_per_block());
  const auto num_energy_bins = 1U; // stir_exam_info.get_num_energy_bins();
  for (int SGID = 0; SGID < num_SGIDs; ++SGID)
    {
      // extract first module_pair for this SGID.
      // TODO should check that it is found
      const auto LUT_iter = std::find(module_pair_SGID_LUT.begin(), module_pair_SGID_LUT.end(), SGID);
      const auto module_pair_idx =
        static_cast<unsigned>(LUT_iter - module_pair_SGID_LUT.begin());
      const auto mod0 = module_pair_idx / num_modules;
      const auto mod1 = module_pair_idx % num_modules;
      assert(module_pair_SGID_LUT(mod0, mod1) == SGID);
      
      petsird::ModulePairEfficiencies module_pair_efficiencies;
      module_pair_efficiencies.values = yardl::NDArray<float, 2>(
                                                                 { num_det_els_in_module * num_energy_bins, num_det_els_in_module * num_energy_bins });

      // find values from STIR
      {
        for (unsigned id0 = 0; id0 < num_det_els_in_module; ++id0)
          {
            const petsird_helpers::ExpandedDetectionBin expanded_det_bin0{mod0, id0, 0};
            const stir::DetectionPosition<> pos0{get_stir_det_pos_from_PETSIRD_id(expanded_det_bin0, stir_scanner)};
            for (unsigned id1 = 0; id1 < num_det_els_in_module; ++id1)
              {           
                const petsird_helpers::ExpandedDetectionBin expanded_det_bin1{mod1, id1, 0};
                const stir::DetectionPosition<> pos1{get_stir_det_pos_from_PETSIRD_id(expanded_det_bin1, stir_scanner)};

                const stir::DetectionPositionPair<> det_pos_pair(pos0, pos1);
                stir::Bin bin;
                float& eff = module_pair_efficiencies.values(id0, id1);
                if (stir_proj_data_info.get_bin_for_det_pos_pair(bin, det_pos_pair).succeeded())
                  {
                    // still need to check if within limits
                    if (bin.tangential_pos_num() >= stir_proj_data_info.get_min_tangential_pos_num() &&
                        bin.tangential_pos_num() <= stir_proj_data_info.get_max_tangential_pos_num() &&
                        bin.segment_num() >= stir_proj_data_info.get_min_segment_num() &&
                        bin.segment_num() <= stir_proj_data_info.get_max_segment_num() &&
                        bin.axial_pos_num() >= stir_proj_data_info.get_min_axial_pos_num(bin.segment_num()) &&
                        bin.axial_pos_num() <= stir_proj_data_info.get_max_axial_pos_num(bin.segment_num()))
                      eff = norm.get_bin_efficiency(bin);
                    else
                      eff = 0.F;
                  }
                else
                  {
                    eff = 0.F;
                  }
              }
          }
      }
      module_pair_efficiencies.sgid = SGID;
      module_pair_efficiencies_vector.emplace_back(module_pair_efficiencies);
    }

}

void
get_detection_efficiencies(petsird::ScannerInformation& scanner,
                           const stir::ProjDataInfo& stir_proj_data_info,
                           const stir::ExamInfo& stir_exam_info,
                           const stir::BinNormalisation& norm)
{
  if (auto stir_proj_data_info_ptr = dynamic_cast<const stir::ProjDataInfoBlocksOnCylindricalNoArcCorr * >(&stir_proj_data_info))
    {
      get_detection_efficiencies_help(scanner,
                                      *stir_proj_data_info_ptr,
                                      stir_exam_info,
                                      norm);
    }
  else if (auto stir_proj_data_info_ptr = dynamic_cast<const stir::ProjDataInfoCylindricalNoArcCorr * >(&stir_proj_data_info))
    {
        get_detection_efficiencies_help(scanner,
                                        *stir_proj_data_info_ptr,
                                        stir_exam_info,
                                        norm);
    }
  else
    {
      stir::error("Cannot only handle *NoArcCorr");
    }

}

petsird::Header
get_header()
{
  petsird::Subject subject;
  subject.id = "123456";
  petsird::Institution institution;
  institution.name = "ETSI Hackathon";
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

template <class ProjDataInfoT>
inline
stir::DetectionPositionPair<>
get_det_pos_pair_help(const stir::ListEvent& event, const ProjDataInfoT& proj_data_info)
{
  stir::DetectionPositionPair<> dp_pair;
  // first check if we can do a faster conversion
  if (auto event_discrete_scanner_ptr = dynamic_cast<stir::CListEventScannerWithDiscreteDetectors<ProjDataInfoT> const*>(&event))
    {
      event_discrete_scanner_ptr->get_detection_position(dp_pair);
    }
  else
    {
      // fall back to more general function
      stir::Bin curr_bin;
      event.get_bin(curr_bin, proj_data_info);
      proj_data_info.get_det_pos_pair_for_bin(dp_pair, curr_bin);
    }
  return dp_pair;
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
  event_time_blk.delayed_events = std::make_optional<std::vector<std::vector<petsird::ListOfCoincidenceEvents>>>();
  (*event_time_blk.delayed_events).resize(1);
  (*event_time_blk.delayed_events)[0].resize(1);


  petsird::ExternalSignalTimeBlock signal_time_blk;
  petsird::BedMovementTimeBlock bed_movement_time_blk;
  petsird::GantryMovementTimeBlock gantry_movement_time_blk;

  std::vector<petsird::CoincidenceEvent> prompts_this_blk;
  std::vector<petsird::CoincidenceEvent> delayeds_this_blk;

  double current_time = 0.0;
  unsigned long num_events = 0;

  // Setup the petsird header info
  petsird::Header header_info = get_header();
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
          (*event_time_blk.delayed_events)[0][0] = delayeds_this_blk;
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
      std::shared_ptr<stir::BinNormalisation> normalisation_sptr;
      stir::KeyParser parser;
      parser.add_start_key("Bin normalisation parameters");
      parser.add_parsing_key("type", &normalisation_sptr);
      parser.add_stop_key("END");
      parser.parse(argv[3]);
      my_class.set_normalisation_sptr(normalisation_sptr);
      
    }
  my_class.process_data();

  return EXIT_SUCCESS;
}
