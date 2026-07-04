/**
* @Function: Ambiguity Resolution
*
* @Author  : Cheng Chi
* @Email   : chichengcn@sjtu.edu.cn
*
* Copyright (C) 2023 by Cheng Chi, All rights reserved.
**/
#include "gici/gnss/ambiguity_resolution.h"

#include <sstream>

#include "gici/gnss/gnss_common.h"
#include "gici/gnss/phaserange_error_sd.h"
#include "gici/gnss/ambiguity_error.h"
#include "gici/gnss/gnss_parameter_blocks.h"

namespace gici {

// The default constructor
AmbiguityResolution::AmbiguityResolution(
    const AmbiguityResolutionOptions options,
    const std::shared_ptr<Graph>& graph) :
  options_(options), graph_(graph),
  cauchy_loss_function_(new ceres::CauchyLoss(1)),
  huber_loss_function_(new ceres::HuberLoss(1))
{
  ambiguities_.push_back(std::vector<Spec>());
  ambiguity_pairs_.push_back(std::vector<BsdPair>());
  ambiguity_lane_pairs_.push_back(std::vector<LanePair>());
}

// The default destructor
AmbiguityResolution::~AmbiguityResolution()
{}

// Solve ambiguity at given epoch
AmbiguityResolution::Result AmbiguityResolution::solvePpp(
             const BackendId& epoch_id, 
             const std::vector<BackendId>& ambiguity_ids,
             const Eigen::MatrixXd& ambiguity_covariance,
             const GnssMeasurement& measurements)
{
  // Prepare data 
  // Shift storage
  ambiguities_.push_back(std::vector<Spec>());
  ambiguity_pairs_.push_back(std::vector<BsdPair>());
  ambiguity_lane_pairs_.push_back(std::vector<LanePair>());
  if (ambiguities_.size() > 2) {
    ambiguities_.pop_front();
    ambiguity_pairs_.pop_front();
    ambiguity_lane_pairs_.pop_front();
  }
  full_parameters_store_.clear();
  lane_groups_.clear();
  current_epoch_id_for_stability_ = epoch_id.asInteger();
  current_ppp_ar_epoch_sequence_ = ++ppp_ar_epoch_sequence_;

  // Get parameters
  std::unordered_map<size_t, size_t> ambiguity_index_map;
  for (size_t i = 0; i < ambiguity_ids.size(); i++) {
    char system = ambiguity_ids[i].gSystem();
    std::string prn = ambiguity_ids[i].gPrn();
    int phase_id = ambiguity_ids[i].gPhaseId();
    if (!gnss_common::useSystem(options_.system_exclude, system)) continue;
    if (!gnss_common::useSatellite(options_.satellite_exclude, prn)) continue;
    if (!gnss_common::usePhase(options_.phase_exclude, system, phase_id)) continue;

    const Eigen::Vector3d& position = measurements.position;
    const Eigen::Vector3d& sat_position = 
      measurements.satellites.at(prn).sat_position;
    double elevation = gnss_common::satelliteElevation(sat_position, position);
    if (gnss_common::radToDegree(elevation) < options_.min_elevation) continue;

    uint64_t id = ambiguity_ids[i].asInteger();
    Spec ambiguity;
    ambiguity.id = ambiguity_ids[i];
    ambiguity.prn = prn;
    CHECK(graph_->parameterBlockExists(id));

    // parameter and residual block
    ambiguity.parameter_block = graph_->parameterBlockPtr(id);
    ambiguity.value = *ambiguity.parameter_block->parameters();
    Graph::ResidualBlockCollection residuals = graph_->residuals(id);
    CHECK(residuals.size() > 0);
    for (size_t r = 0; r < residuals.size(); ++r) {
      if (residuals[r].error_interface_ptr->typeInfo() == ErrorType::kPhaserangeError ||
          residuals[r].error_interface_ptr->typeInfo() == ErrorType::kPhaserangeErrorSD ||
          residuals[r].error_interface_ptr->typeInfo() == ErrorType::kPhaserangeErrorDD) {
        ambiguity.residual_block = residuals[r];
        break;
      }
    }

    // wavelength
    const Satellite& satellite = measurements.satellites.at(prn);
    for (auto obs : satellite.observations) {
      auto& observation = obs.second;
      if (gnss_common::getPhaseID(system, obs.first) == phase_id) {
        ambiguity.wavelength = observation.wavelength;
      }
    }

    // elevation angle for reference satellite selection
    ambiguity.elevation = elevation;

    curAmbs().push_back(ambiguity);
    ambiguity_index_map.insert(std::make_pair(curAmbs().size() - 1, i));
  }

  // Get covariance of ambiguity
  CHECK(ambiguity_covariance.cols() == ambiguity_ids.size());
  const size_t ambiguity_size = curAmbs().size();
  ambiguity_covariance_.resize(ambiguity_size, ambiguity_size);
  for (size_t i = 0; i < ambiguity_size; i++) {
    for (size_t j = 0; j < ambiguity_size; j++) {
      const size_t raw_i = ambiguity_index_map.at(i);
      const size_t raw_j = ambiguity_index_map.at(j);
      ambiguity_covariance_(i, j) = ambiguity_covariance(raw_i, raw_j);
    }
    curAmbs()[i].std = sqrt(ambiguity_covariance_(i, i));
  }

  // Collect all parameter blocks
  Graph::ParameterBlockCollection parameters = graph_->parameters();
  for (size_t p = 0; p < parameters.size(); p++) {
    Parameter parameter;
    parameter.id = BackendId(parameters[p].first);
    parameter.size = parameters[p].second->dimension();
    parameter.minimal_size = parameters[p].second->minimalDimension();
    parameter.value = Eigen::Map<Eigen::VectorXd>(
      parameters[p].second->parameters(), parameter.size);
    parameter.handle = parameters[p].second;
    if (full_parameters_store_.size() == 0) {
      parameter.covariance_start_index = 0;
    }
    else {
      Parameter& last_parameter = full_parameters_store_.back();
      parameter.covariance_start_index = 
        last_parameter.covariance_start_index + last_parameter.minimal_size;
    }
    full_parameters_store_.push_back(parameter);
  }

  // Sovle ambiguity
  // Apply Between-Satellite-Difference (BSD) and lane combination
  formSatellitePairPpp();

  // Sort lanes to groups
  groupingSatellitePair();
  if (lane_groups_.size() == 0) return Result::NoFix;

  // Try fix lanes
  int num_success_uwl = 0, num_candidate_uwl = 0; 
  int num_success_wl = 0, num_candidate_wl = 0; 
  int num_success_nl = 0, num_candidate_nl = 0;
  bool has_uwl = false, has_wl = false;
  for (int i = lane_groups_.size() - 1; i >= 0; i--) {
    LaneType type = curAmbLanePairs().at(lane_groups_.at(i).front()).laneType();
    double min_percentage_fixation;
    bool use_rounding = false;
    if (type == LaneType::UWL) {
      min_percentage_fixation = options_.min_percentage_fixation_uwl;
      use_rounding = true;
      has_uwl = true;
    }
    else if (type == LaneType::WL) {
      min_percentage_fixation = options_.min_percentage_fixation_wl;
      use_rounding = true;
      has_wl = true;
    }
    else {
      min_percentage_fixation = options_.min_percentage_fixation_nl;
      use_rounding = false;
    }

    int num_fixed = trySolveLanes(
      lane_groups_.at(i), min_percentage_fixation, use_rounding);

    if (type == LaneType::UWL) {
      num_success_uwl += num_fixed;
      num_candidate_uwl += lane_groups_.at(i).size();
    }
    else if (type == LaneType::WL) {
      num_success_wl += num_fixed;
      num_candidate_wl += lane_groups_.at(i).size();
    }
    else {
      num_success_nl += num_fixed;
      num_candidate_nl += lane_groups_.at(i).size();
    }
  }

  // WL/UWL are only auxiliary for PPP AR. Do not keep them if NL is not fixed.
  if (num_success_nl == 0) {
    setGraphParameters(full_parameters_store_);
    eraseAmbiguityResidualBlocks(curAmbLanePairs());
    refreshCurrentAmbiguityValuesFromGraph();
    if (tryApplyHeldPppNlFixes(epoch_id, measurements)) {
      return Result::NlFix;
    }
    return Result::NoFix;
  }

  int num_removed_auxiliary_constraints = 0;
  for (auto& lane_pair : curAmbLanePairs()) {
    if (lane_pair.laneType() == LaneType::NL) continue;
    if (lane_pair.residual_id) {
      graph_->removeResidualBlock(lane_pair.residual_id);
      lane_pair.residual_id = nullptr;
      lane_pair.is_fixed = false;
      num_removed_auxiliary_constraints++;
    }
  }
  LOG(INFO) << "PPP AR final graph keeps NL only: removed_auxiliary_constraints="
            << num_removed_auxiliary_constraints;

  if (!validatePppFixedSolution(epoch_id, "fresh")) {
    if (tryApplyHeldPppNlFixes(epoch_id, measurements)) {
      return Result::NlFix;
    }
    return Result::NoFix;
  }
  updateHeldPppNlFixesFromAccepted();

  // Check status
  return Result::NlFix;
}

// Apply Between-Satellite-Difference (BSD) and lane combination for PPP
void AmbiguityResolution::formSatellitePairPpp()
{
  // Prepare data
  std::map<char, int> system_to_num_phases;
  std::multimap<char, double> system_to_wavelengths;
  std::map<std::string, int> prn_to_number_phases; 
  std::multimap<std::string, double> prn_to_wavelengths;
  std::multimap<std::string, Spec> prn_to_specs;
  std::map<uint64_t, size_t> id_to_spec_id;
  bool has_bds_meo = false;  // we prefer to select BDS MEO satellite as reference
  for (size_t i = 0; i < curAmbs().size(); i++) {
    Spec& ambiguity = curAmbs()[i];
    std::string prn = ambiguity.id.gPrn();

    auto it = prn_to_number_phases.find(prn);
    if (it == prn_to_number_phases.end()) {
      prn_to_number_phases.insert(std::make_pair(prn, 1));
    }
    else it->second++;
    prn_to_specs.insert(std::make_pair(prn, ambiguity));
    prn_to_wavelengths.insert(std::make_pair(prn, ambiguity.wavelength));
    id_to_spec_id.insert(std::make_pair(ambiguity.id.asInteger(), i));

    // check if we have a BDS MEO satellite
    if (prn[0] == 'C' && !gnss_common::isBdsMeo(prn)) {
      has_bds_meo = true;
    }
  }
  for (size_t i = 0; i < getGnssSystemList().size(); i++) {
    char system = getGnssSystemList()[i];

    system_to_num_phases.insert(std::make_pair(system, 0));
    for (auto it : prn_to_number_phases) {
      if (it.first[0] != system) continue;
      if (system == 'C' && has_bds_meo && 
          !gnss_common::isBdsMeo(it.first)) continue;
      if (system_to_num_phases.at(system) < it.second) {
        system_to_num_phases.at(system) = it.second;
      }
    }

    for (auto it : prn_to_wavelengths) {
      if (it.first[0] != system) continue;
      if (system_to_wavelengths.find(system) == system_to_wavelengths.end()) {
        system_to_wavelengths.insert(std::make_pair(system, it.second));
      }
      bool found = false;
      for (auto it_wave = system_to_wavelengths.lower_bound(system); 
          it_wave != system_to_wavelengths.upper_bound(system); it_wave++) {
        if (it_wave->second == it.second) {
          found = true; break;
        }
      }
      if (!found) system_to_wavelengths.insert(std::make_pair(system, it.second));
    }
  }

  // Find reference satellites for each system and frequencies
  std::map<char, std::string> system_to_reference_prn;
  for (size_t i = 0; i < getGnssSystemList().size(); i++) {
    char system = getGnssSystemList()[i];

    // find satellite with maximum elevation angle
    double max_elevation = 0.0;
    for (size_t j = 0; j < curAmbs().size(); j++) {
      Spec& ambiguity = curAmbs()[j];
      if (ambiguity.id.gSystem() != system) continue;

      // we prefer to use BDS MEO satellites
      if (system == 'C' && has_bds_meo && 
          !gnss_common::isBdsMeo(ambiguity.prn)) continue;

      // we only select satellites with max phase number
      if (prn_to_number_phases.at(ambiguity.id.gPrn()) != 
          system_to_num_phases.at(system)) continue;

      if (max_elevation < ambiguity.elevation) {
        system_to_reference_prn[system] = ambiguity.id.gPrn();
        max_elevation = ambiguity.elevation;
      }
    }
  }

  // Form BSD pair
  for (size_t i = 0; i < curAmbs().size(); i++) {
    Spec& ambiguity = curAmbs()[i];
    char system = ambiguity.id.gSystem();
    std::string prn = ambiguity.id.gPrn();
    if (system_to_reference_prn.find(system) == 
        system_to_reference_prn.end()) continue;
    std::string prn_ref = system_to_reference_prn.at(system);

    if (prn == prn_ref) {
      ambiguity.is_reference = true;
      continue;
    }

    for (auto it = prn_to_specs.lower_bound(prn_ref); 
         it != prn_to_specs.upper_bound(prn_ref); it++) {
      Spec& ambiguity_ref = it->second;
      if (ambiguity_ref.id.gPhaseId() == ambiguity.id.gPhaseId()) {
        BsdPair ambiguity_pair(curAmbs(), i, id_to_spec_id.at(ambiguity_ref.id.asInteger()));
        curAmbPairs().push_back(ambiguity_pair);
        break;
      }
    }
  }

  // Widelane combination 
  // We set the phase with highest frequency as base, and sequentially apply dual-freqency 
  // widelane combination between two adjacent frequecies from the lowest frequecy to base.
  // Sort frequencies for each system
  std::map<char, double> system_to_base_wavelength;
  for (size_t i = 0; i < getGnssSystemList().size(); i++) {
    char system = getGnssSystemList()[i];
    std::vector<double> wavelengths;
    auto it = system_to_wavelengths.lower_bound(system);
    if (it == system_to_wavelengths.end()) continue;
    for (; it != system_to_wavelengths.upper_bound(system); it++) {
      wavelengths.push_back(it->second);
    }
    if (wavelengths.size() == 0) continue;
    system_to_wavelengths.erase(system);
    std::sort(wavelengths.begin(), wavelengths.end());
    for (std::vector<double>::reverse_iterator wavelength = wavelengths.rbegin(); 
         wavelength != wavelengths.rend(); wavelength++) {
      system_to_wavelengths.insert(std::make_pair(system, *wavelength));
    }
    system_to_base_wavelength.insert(std::make_pair(system, wavelengths.front()));
  }

  // Form combinations
  for (size_t i = 0; i < curAmbPairs().size(); i++) {
    char system = curAmbs()[curAmbPairs()[i].spec_id].id.gSystem();
    std::string prn = curAmbs()[curAmbPairs()[i].spec_id].id.gPrn();

    // base frequency observation for NL AR
    if (curAmbPairs()[i].wavelength == system_to_base_wavelength.at(system)) {
      LanePair nl_pair(curAmbPairs(), i);
      curAmbLanePairs().push_back(nl_pair);
      curAmbPairs()[i].is_base_frequency = true;
    }

    // Widelanes
    double wavelength = curAmbPairs()[i].wavelength;
    // select the nearest frequency to combine
    double desired_wavelength_to_combine = 1.0e6;
    for (auto it = system_to_wavelengths.lower_bound(system); 
         it != system_to_wavelengths.upper_bound(system); it++) {
      if (it->second <= wavelength) continue;
      if (it->second < desired_wavelength_to_combine) {
        desired_wavelength_to_combine = it->second;
      }
    } 
    if (desired_wavelength_to_combine == 1.0e6) continue;
    for (size_t j = 0; j < curAmbPairs().size(); j++) {
      char system_j = curAmbs()[curAmbPairs()[j].spec_id].id.gSystem();
      if (system_j != system) continue;

      std::string prn_j = curAmbs()[curAmbPairs()[j].spec_id].id.gPrn();
      if (prn_j != prn) continue;
      
      double wavelength_j = curAmbPairs()[j].wavelength;
      if (wavelength_j != desired_wavelength_to_combine) continue;

      LanePair wl_pair(curAmbPairs(), i, j);
      curAmbLanePairs().push_back(wl_pair);
    }
  }
}

// Sort lanes to groups
void AmbiguityResolution::groupingSatellitePair()
{
  // Grouping lanes
  const size_t min_num = options_.min_num_satellite_pairs_fixation;
  pushBatchBack(lane_groups_, groupingSatellitePairForBand(LaneType::NL, min_num));
  pushBatchBack(lane_groups_, groupingSatellitePairForBand(LaneType::WL, min_num));
  pushBatchBack(lane_groups_, groupingSatellitePairForBand(LaneType::UWL, min_num));

  // Check if we need clear invalid groups
  while (lane_groups_.size() > 0 && lane_groups_.back().size() < min_num) {
    lane_groups_.pop_back();
  }
} 

// Sort lanes to groups for given band
std::vector<std::vector<size_t>> 
AmbiguityResolution::groupingSatellitePairForBand(
  const LaneType& type, const size_t min_num)
{
  // Get all wavelengths
  std::multimap<double, size_t> wavelength_to_indexes;
  std::unordered_map<size_t, char> index_to_systems;
  for (size_t i = 0; i < curAmbLanePairs().size(); i++) {
    if (curAmbLanePairs()[i].laneType() != type) continue;
    double wavelength = curAmbLanePairs()[i].wavelength;
    size_t id = curAmbPairs()[curAmbLanePairs()[i].bsd_pair_id_higher].spec_id_ref;
    char system = curAmbs()[id].prn[0];
    wavelength_to_indexes.insert(std::make_pair(wavelength, i));
    index_to_systems.insert(std::make_pair(i, system));
  }

  // The lanes with the same wavelength will be putted into one group. 
  // The lane number of every groups should be no less than min_num.
  double current_wavelength = 0.0;
  size_t current_group_index = 0;
  char current_system = 0;
  std::vector<std::vector<size_t>> lane_groups;
  std::vector<std::vector<double>> lane_group_wavelengths;
  for (auto it = wavelength_to_indexes.begin(); 
       it != wavelength_to_indexes.end();) {
    double wavelength = it->first;
    size_t index = it->second;
    char system = index_to_systems.at(index);
    if (current_wavelength == 0.0) current_wavelength = wavelength;
    if (current_system == 0) current_system = system;
    // grouping current wavelength
    if (wavelength == current_wavelength && current_system == system) {
      if (lane_groups.size() <= current_group_index) {
        lane_groups.push_back(std::vector<size_t>());
        lane_group_wavelengths.push_back(std::vector<double>());
        lane_group_wavelengths.back().push_back(wavelength);
      }
      std::vector<size_t>& group = lane_groups.at(current_group_index);
      group.push_back(index);
      it = wavelength_to_indexes.erase(it);
      continue;
    }
    else {
      // check if we have sufficient lanes in current group
      std::vector<size_t>& group = lane_groups.at(current_group_index);
      if (group.size() < min_num) {
        // merge into the last group, or continue adding lanes behind
        bool do_merge = false;
        if (current_group_index == 0) {
          do_merge = false;
        }
        else {
          double last_group_wavelength = 
            lane_group_wavelengths.at(current_group_index - 1).back();
          double current_group_wavelength = 
            lane_group_wavelengths.at(current_group_index).front();
          double dw_last = current_group_wavelength - last_group_wavelength;
          double dw_next = wavelength - current_group_wavelength;
          if (dw_last < dw_next) do_merge = true;
          else do_merge = false;
        }
        if (do_merge) {
          pushBatchBack(lane_groups.at(current_group_index - 1), group);
          pushBatchBack(lane_group_wavelengths.at(current_group_index - 1), 
                        lane_group_wavelengths.at(current_group_index));
          lane_group_wavelengths.pop_back();
          lane_groups.pop_back();
        }
      }
      else {
        current_group_index++;
      }
      current_wavelength = wavelength;
      current_system = system;
      continue;
    }
    it++;
  }
  // check if the last group has insufficient number of lanes
  if (lane_groups.size() > 1 && lane_groups.back().size() < min_num) {
    pushBatchBack(lane_groups.at(current_group_index - 1), lane_groups.back());
    pushBatchBack(lane_group_wavelengths.at(current_group_index - 1), 
                  lane_group_wavelengths.at(current_group_index));
    lane_group_wavelengths.pop_back();
    lane_groups.pop_back();
  }

  return lane_groups;
}

// Solve ambiguities
int AmbiguityResolution::solveLanes(
  const std::vector<size_t>& indexes, 
  const double min_percentage_fixation,
  const DelMethod method,
  const bool use_rounding,
  const bool skip_full)
{
  // Get combinations
  std::vector<LanePair> lane_pairs;
  for (size_t i = 0; i < indexes.size(); i++) {
    lane_pairs.push_back(curAmbLanePairs().at(indexes[i]));
  }
  LaneType lane_type = lane_pairs.front().laneType();

  // Sort them for partial AR
  // by variance
  if (method == DelMethod::Variance) {
    std::sort(lane_pairs.begin(), lane_pairs.end(), [](LanePair lhs, LanePair rhs) 
      -> bool { return lhs.std < rhs.std; });
  }
  // by elevation
  else if (method == DelMethod::Elevation) {
    std::sort(lane_pairs.begin(), lane_pairs.end(), [](LanePair lhs, LanePair rhs) 
      -> bool { return lhs.elevation > rhs.elevation; });
  }
  // by frational part
  else if (method == DelMethod::Fractional) {
    std::sort(lane_pairs.begin(), lane_pairs.end(), [](LanePair lhs, LanePair rhs) 
      -> bool { 
        const double cycle_lhs = lhs.value / lhs.wavelength;
        const double cycle_rhs = rhs.value / rhs.wavelength;
        const double fraction_lhs = fabs(cycle_lhs - round(cycle_lhs));
        const double fraction_rhs = fabs(cycle_rhs - round(cycle_rhs));
        return fraction_lhs < fraction_rhs; });
  }

  // Get float ambiguities and its covariance
  Eigen::VectorXd float_ambiguities = Eigen::VectorXd::Zero(lane_pairs.size());
  Eigen::MatrixXd float_covariance = 
    Eigen::MatrixXd::Zero(lane_pairs.size(), lane_pairs.size());
  Eigen::MatrixXd differential_jacobian = 
    Eigen::MatrixXd::Zero(lane_pairs.size(), curAmbs().size());
  for (size_t i = 0; i < lane_pairs.size(); i++) {
    float_ambiguities(i) = lane_pairs[i].value / lane_pairs[i].wavelength;
    size_t id_higher_raw = curAmbPairs()[lane_pairs[i].bsd_pair_id_higher].spec_id;
    size_t id_higher_ref = curAmbPairs()[lane_pairs[i].bsd_pair_id_higher].spec_id_ref;
    const LaneType lane_type_i = lane_pairs[i].laneType();
    const double fractional = fabs(float_ambiguities(i) - round(float_ambiguities(i)));
    const char* lane_type_name = lane_type_i == LaneType::NL ? "NL" :
      (lane_type_i == LaneType::WL ? "WL" : "UWL");

    if (lane_type_i != LaneType::NL) {
      size_t id_lower_raw = curAmbPairs()[lane_pairs[i].bsd_pair_id_lower].spec_id;
      size_t id_lower_ref = curAmbPairs()[lane_pairs[i].bsd_pair_id_lower].spec_id_ref;
      double wave_higher = curAmbPairs()[lane_pairs[i].bsd_pair_id_higher].wavelength;
      double wave_lower = curAmbPairs()[lane_pairs[i].bsd_pair_id_lower].wavelength;
      differential_jacobian(i, id_higher_raw) = 1.0 / wave_higher;
      differential_jacobian(i, id_higher_ref) = -1.0 / wave_higher;
      differential_jacobian(i, id_lower_raw) = -1.0 / wave_lower;
      differential_jacobian(i, id_lower_ref) = 1.0 / wave_lower;
      LOG(INFO) << "AR float ambiguity: lane=" << lane_type_name
                << ", high=" << curAmbs()[id_higher_raw].id.gPrn()
                << "-" << curAmbs()[id_higher_ref].id.gPrn()
                << "/L" << curAmbs()[id_higher_raw].id.gPhaseId()
                << ", low=" << curAmbs()[id_lower_raw].id.gPrn()
                << "-" << curAmbs()[id_lower_ref].id.gPrn()
                << "/L" << curAmbs()[id_lower_raw].id.gPhaseId()
                << ", cycle=" << float_ambiguities(i)
                << ", frac=" << fractional
                << ", wavelength=" << lane_pairs[i].wavelength
                << ", std=" << lane_pairs[i].std / lane_pairs[i].wavelength;
    }
    else {
      double wave = lane_pairs[i].wavelength;
      differential_jacobian(i, id_higher_raw) = 1.0 / wave;
      differential_jacobian(i, id_higher_ref) = -1.0 / wave;
      LOG(INFO) << "AR float ambiguity: lane=" << lane_type_name
                << ", pair=" << curAmbs()[id_higher_raw].id.gPrn()
                << "-" << curAmbs()[id_higher_ref].id.gPrn()
                << "/L" << curAmbs()[id_higher_raw].id.gPhaseId()
                << ", cycle=" << float_ambiguities(i)
                << ", frac=" << fractional
                << ", wavelength=" << lane_pairs[i].wavelength
                << ", std=" << lane_pairs[i].std / lane_pairs[i].wavelength;
    }
  }
  float_covariance = differential_jacobian * 
      ambiguity_covariance_ * differential_jacobian.transpose();

  // Sovle ambiguity
  Eigen::VectorXd fixed_ambiguities;
  int num_active = float_ambiguities.size();
  int min_num_fixation = static_cast<int>(
    min_percentage_fixation * static_cast<double>(num_active));
  if (min_num_fixation < options_.min_num_satellite_pairs_fixation) {
    min_num_fixation = options_.min_num_satellite_pairs_fixation;
  }
  Eigen::VectorXd active_float_ambiguities = float_ambiguities;
  Eigen::MatrixXd active_float_covariance = float_covariance;
  // solve by LAMBDA
  if (!use_rounding) {
    // Try full AR and then partial AR, until it successed or active number 
    // of ambiguities reaches the minimum number of fixation judgement.
    double ratio = 0.0;
    while (num_active >= min_num_fixation) {
      if (!(skip_full && num_active == float_ambiguities.size())) {
        if (sqrt(active_float_covariance(num_active - 1, num_active - 1)) < 0.25)
        if (solveAmbiguityLambda(active_float_ambiguities, 
            active_float_covariance, options_.ratio, fixed_ambiguities, ratio)) break;
      }
      // reduce subsets
      --num_active;
      active_float_ambiguities = float_ambiguities.topRows(num_active);
      active_float_covariance = float_covariance.topLeftCorner(num_active, num_active);
    }
    // if (ratio != 0.0) LOG(INFO) << (int)lane_type << ": " << ratio << " | " << active_float_ambiguities.transpose();
  }
  // solve by rounding
  else {
    num_active = 0;
    for (size_t i = 0; i < float_ambiguities.size(); i++) {
      double integer = round(active_float_ambiguities(i));
      double fractional = active_float_ambiguities(i) - integer;
      if (sqrt(float_covariance(i, i)) > 0.25) break;
      if (fabs(fractional) > 0.25) break;
      fixed_ambiguities.conservativeResize(num_active + 1);
      fixed_ambiguities(num_active) = integer;
      num_active++;
    }
  }

  if (num_active < min_num_fixation) return 0;

  if (lane_type == LaneType::NL && options_.min_num_continuous_nl_fix > 1) {
    int num_stable = 0;
    for (int i = 0; i < num_active; i++) {
      const LanePair& lane_pair = lane_pairs[i];
      const BsdPair& bsd_pair = curAmbPairs().at(lane_pair.bsd_pair_id_higher);
      const Spec& raw = curAmbs().at(bsd_pair.spec_id);
      const Spec& ref = curAmbs().at(bsd_pair.spec_id_ref);
      const long long integer = static_cast<long long>(llround(fixed_ambiguities(i)));
      std::ostringstream oss;
      oss << raw.id.gSystem() << ":" << raw.id.gPrn() << "-"
          << ref.id.gPrn() << ":L" << raw.id.gPhaseId();
      const std::string key = oss.str();

      IntegerStability& stability = nl_integer_stability_[key];
      if (stability.last_epoch_id == current_epoch_id_for_stability_) {
        if (stability.integer != integer) {
          stability.integer = integer;
          stability.count = 1;
        }
      }
      else {
        if (stability.count > 0 && stability.integer == integer) {
          stability.count++;
        }
        else {
          stability.integer = integer;
          stability.count = 1;
        }
        stability.last_epoch_id = current_epoch_id_for_stability_;
      }

      if (stability.integer == integer &&
          stability.count >= options_.min_num_continuous_nl_fix) {
        num_stable++;
      }
      else {
        break;
      }
    }

    if (num_stable < num_active) {
      LOG(INFO) << "PPP AR NL integer stability: proposed=" << num_active
                << ", stable_prefix=" << num_stable
                << ", required_continuous="
                << options_.min_num_continuous_nl_fix;
    }
    if (num_stable < min_num_fixation) return 0;
    if (num_stable < num_active) {
      num_active = num_stable;
      fixed_ambiguities = fixed_ambiguities.head(num_active).eval();
    }
  }

  // Contraint the fixed ambiguities to temporary parameters and check phaserange residual. 
  // ambiguity parameters
  Eigen::VectorXd ambiguity_parameters;
  ambiguity_parameters.resize(curAmbs().size());
  for (size_t i = 0; i < curAmbs().size(); i++) {
    ambiguity_parameters(i) = curAmbs()[i].value;
  }

  // ambiguity measurements
  CHECK(fixed_ambiguities.size() == num_active);
  size_t ambiguity_size = ambiguity_covariance_.cols();
  Eigen::MatrixXd jacobian; jacobian.resize(num_active, ambiguity_size);
  jacobian = differential_jacobian.topRows(num_active);
  Eigen::MatrixXd fix_ambiguity_covariance = 
    Eigen::MatrixXd::Identity(num_active, num_active) * 1e-6; // 0.001 cycles

  // apply Kalman upadte to constraint the parameters
  Eigen::MatrixXd kalman_gain = ambiguity_covariance_ * jacobian.transpose() * (jacobian * 
    ambiguity_covariance_ * jacobian.transpose() + fix_ambiguity_covariance).inverse();
  Eigen::VectorXd dx = kalman_gain * 
    (fixed_ambiguities - jacobian * ambiguity_parameters);
  ambiguity_covariance_ = (Eigen::MatrixXd::Identity(ambiguity_size, ambiguity_size) - 
    kalman_gain * jacobian) * ambiguity_covariance_;
  ambiguity_covariance_ = (ambiguity_covariance_ + ambiguity_covariance_.transpose()) / 2.0;
  // update dx
  ambiguity_parameters += dx;

  // Put ambiguities on global parameters
  // set to ambiguity handles
  auto& ambiguities = curAmbs();
  auto& ambiguity_pairs = curAmbPairs();
  auto& ambiguity_lane_pairs = curAmbLanePairs();
  for (size_t i = 0; i < ambiguities.size(); i++) {
    ambiguities[i].value = ambiguity_parameters(i);
  }
  for (size_t i = 0; i < ambiguity_pairs.size(); i++) {
    ambiguity_pairs[i].update(ambiguities);
  }
  for (size_t i = 0; i < ambiguity_lane_pairs.size(); i++) {
    ambiguity_lane_pairs[i].update(ambiguity_pairs);
  }
  // set fix information
  for (int i = 0; i < num_active; i++) {
    auto& fixed_lane_pair = lane_pairs[i];
    for (size_t j = 0; j < ambiguity_lane_pairs.size(); j++) {
      if (ambiguity_lane_pairs[j].bsd_pair_id_higher == 
          fixed_lane_pair.bsd_pair_id_higher && 
          ambiguity_lane_pairs[j].bsd_pair_id_lower == 
          fixed_lane_pair.bsd_pair_id_lower &&
          ambiguity_lane_pairs[j].laneType() == fixed_lane_pair.laneType()) {
        ambiguity_lane_pairs[j].is_fixed = true;
        ambiguity_lane_pairs[j].fixed_integer =
          static_cast<long long>(llround(fixed_ambiguities(i)));
        if (ambiguity_lane_pairs[j].laneType() == LaneType::NL) {
          size_t id = ambiguity_lane_pairs[j].bsd_pair_id_higher;
          ambiguity_pairs[id].is_fixed = true;
        }
        break;
      }
    }
  }

  // Add to graph
  const double information = 1.0e6; 
  for (int i = 0; i < num_active; i++) {
    auto& lane_pair = lane_pairs[i];
    std::vector<double> coefficients;
    size_t id_higher_raw = 
      ambiguity_pairs[lane_pair.bsd_pair_id_higher].spec_id;
    size_t id_higher_ref = 
      ambiguity_pairs[lane_pair.bsd_pair_id_higher].spec_id_ref;

    if (lane_pair.laneType() != LaneType::NL) {
      size_t id_lower_raw = 
        ambiguity_pairs[lane_pair.bsd_pair_id_lower].spec_id;
      size_t id_lower_ref = 
        ambiguity_pairs[lane_pair.bsd_pair_id_lower].spec_id_ref;
      double wave_higher = 
        ambiguity_pairs[lane_pair.bsd_pair_id_higher].wavelength;
      double wave_lower = 
        ambiguity_pairs[lane_pair.bsd_pair_id_lower].wavelength;
      coefficients.push_back(1.0 / wave_higher);
      coefficients.push_back(-1.0 / wave_higher);
      coefficients.push_back(-1.0 / wave_lower);
      coefficients.push_back(1.0 / wave_lower);
      std::shared_ptr<AmbiguityError4Coef> ambiguity_error = 
        std::make_shared<AmbiguityError4Coef>(
        fixed_ambiguities(i), information, coefficients);
        CHECK(lane_pair.residual_id == nullptr);
      lane_pair.residual_id = graph_->addResidualBlock(
        ambiguity_error,
        cauchy_loss_function_ ? cauchy_loss_function_.get() : nullptr,
        graph_->parameterBlockPtr(curAmbs()[id_higher_raw].id.asInteger()),
        graph_->parameterBlockPtr(curAmbs()[id_higher_ref].id.asInteger()),
        graph_->parameterBlockPtr(curAmbs()[id_lower_raw].id.asInteger()),
        graph_->parameterBlockPtr(curAmbs()[id_lower_ref].id.asInteger()));
    }
    else {
      double wave = lane_pair.wavelength;
      coefficients.push_back(1.0 / wave);
      coefficients.push_back(-1.0 / wave);

      std::shared_ptr<AmbiguityError2Coef> ambiguity_error = 
        std::make_shared<AmbiguityError2Coef>(
        fixed_ambiguities(i), information, coefficients);
      CHECK(lane_pair.residual_id == nullptr);
      lane_pair.residual_id = graph_->addResidualBlock(
        ambiguity_error,
        cauchy_loss_function_ ? cauchy_loss_function_.get() : nullptr,
        graph_->parameterBlockPtr(curAmbs()[id_higher_raw].id.asInteger()),
        graph_->parameterBlockPtr(curAmbs()[id_higher_ref].id.asInteger()));
    }

    bool residual_id_written = false;
    for (size_t j = 0; j < ambiguity_lane_pairs.size(); j++) {
      if (ambiguity_lane_pairs[j].bsd_pair_id_higher ==
          lane_pair.bsd_pair_id_higher &&
          ambiguity_lane_pairs[j].bsd_pair_id_lower ==
          lane_pair.bsd_pair_id_lower &&
          ambiguity_lane_pairs[j].laneType() == lane_pair.laneType()) {
        ambiguity_lane_pairs[j].residual_id = lane_pair.residual_id;
        ambiguity_lane_pairs[j].fixed_integer =
          static_cast<long long>(llround(fixed_ambiguities(i)));
        residual_id_written = true;
        break;
      }
    }
    CHECK(residual_id_written);
  }

  return num_active;
}

// Try to solve ambiguity by all partial AR methods
int AmbiguityResolution::trySolveLanes(
  const std::vector<size_t>& indexes, 
  const double min_percentage_fixation,
  const bool use_rounding)
{
  static const std::vector<DelMethod> sequence = {
    DelMethod::Elevation, DelMethod::Variance, DelMethod::Fractional};
  // static const std::vector<DelMethod> sequence = {DelMethod::Elevation};
  int num_fixed = 0;
  for (size_t i = 0; i < sequence.size(); i++) {
    if (num_fixed) break;
    DelMethod method = sequence[i];
    bool skip_full = i == 0 ? false : true;
    num_fixed = solveLanes(indexes, min_percentage_fixation, 
      method, use_rounding, skip_full);
  }

  return num_fixed;
}

// Search match on the last pairs
bool AmbiguityResolution::findMatch(
    LanePair& lane_pair,
    std::vector<LanePair>& matches,
    std::vector<double>& coefficients)
{
  if (!lane_pair.is_fixed) return false;

  matches.clear();
  coefficients.clear();

  auto& pair_higher = curAmbPairs()[lane_pair.bsd_pair_id_higher];
  auto& spec_higher_raw = curAmbs()[pair_higher.spec_id];
  auto& spec_higher_ref = curAmbs()[pair_higher.spec_id_ref];
  char system = spec_higher_raw.id.gSystem();
  std::string prn_raw = spec_higher_raw.id.gPrn();
  std::string prn_ref = spec_higher_ref.id.gPrn();
  double wavelength = lane_pair.wavelength;
  double phase_id_higher = spec_higher_raw.id.gPhaseId();
  double phase_id_lower = 0.0;
  if (lane_pair.laneType() != LaneType::NL) {
    auto& pair_lower = curAmbPairs()[lane_pair.bsd_pair_id_lower];
    auto& spec_lower_raw = curAmbs()[pair_lower.spec_id];
    auto& spec_lower_ref = curAmbs()[pair_lower.spec_id_ref];
    phase_id_lower = spec_lower_raw.id.gPhaseId();
  }
  auto& last_lane_pairs = lastAmbLanePairs();

  for (size_t i = 0; i < last_lane_pairs.size(); i++) {
    if (last_lane_pairs[i].wavelength != wavelength) continue;
    if (!last_lane_pairs[i].is_fixed) continue;

    auto& last_pair_higher = lastAmbPairs()[last_lane_pairs[i].bsd_pair_id_higher];
    auto& last_spec_higher_raw = lastAmbs()[last_pair_higher.spec_id];
    auto& last_spec_higher_ref = lastAmbs()[last_pair_higher.spec_id_ref];
    char last_system = last_spec_higher_raw.id.gSystem();
    std::string last_prn_raw = last_spec_higher_raw.id.gPrn();
    std::string last_prn_ref = last_spec_higher_ref.id.gPrn();
    double last_phase_id_higher = last_spec_higher_raw.id.gPhaseId();
    double last_phase_id_lower = 0.0;
    if (last_lane_pairs[i].laneType() != LaneType::NL) {
      auto& last_pair_lower = lastAmbPairs()[last_lane_pairs[i].bsd_pair_id_lower];
      auto& last_spec_lower_raw = lastAmbs()[last_pair_lower.spec_id];
      auto& last_spec_lower_ref = lastAmbs()[last_pair_lower.spec_id_ref];
      last_phase_id_lower = last_spec_lower_raw.id.gPhaseId();
    }

    if (last_system != system) continue;

    if (lane_pair.laneType() != LaneType::NL) {
      if (last_phase_id_higher != phase_id_higher || 
          last_phase_id_lower != phase_id_lower) {
        continue;
      }
    }
    else {
      if (last_phase_id_higher != phase_id_higher) {
        continue;
      }
    }

    if (last_prn_raw == prn_raw && last_prn_ref == prn_ref) {
      matches.push_back(last_lane_pairs[i]);
      coefficients.push_back(1.0);
      return true;
    }

    // reference satellite changed
    if (prn_ref != last_prn_ref) {
      // last reference satellite is current raw satellite,
      // we find another satellite
      if (prn_raw == last_prn_ref && prn_ref == last_prn_raw) {
        matches.push_back(last_lane_pairs[i]);
        coefficients.push_back(-1.0);
        return true;
      }
      // not that case
      // we find all candidates with the satellite PRNs
      else {
        if (prn_raw == last_prn_raw) {
          matches.push_back(last_lane_pairs[i]);
          coefficients.push_back(1.0);
        }
        if (prn_ref == last_prn_raw) {
          matches.push_back(last_lane_pairs[i]);
          coefficients.push_back(-1.0);
        }
      }
    }
  }

  // Check candidates
  CHECK(matches.size() == coefficients.size());
  CHECK(matches.size() < 3);
  // we did not find the pair
  if (matches.size() < 2) return false;
  else return true;
}

// Get a stable key for a PPP NL BSD lane
std::string AmbiguityResolution::pppNlLaneKey(const LanePair& lane_pair)
{
  CHECK(lane_pair.laneType() == LaneType::NL);
  const BsdPair& bsd_pair = curAmbPairs().at(lane_pair.bsd_pair_id_higher);
  const Spec& raw = curAmbs().at(bsd_pair.spec_id);
  const Spec& ref = curAmbs().at(bsd_pair.spec_id_ref);

  std::ostringstream oss;
  oss << raw.id.gSystem() << ":" << raw.id.gPrn() << "-"
      << ref.id.gPrn() << ":L" << raw.id.gPhaseId();
  return oss.str();
}

// Check whether either satellite of a PPP NL BSD lane has a phase slip
bool AmbiguityResolution::pppNlLaneHasSlip(
  const LanePair& lane_pair, const GnssMeasurement& measurements)
{
  CHECK(lane_pair.laneType() == LaneType::NL);
  const BsdPair& bsd_pair = curAmbPairs().at(lane_pair.bsd_pair_id_higher);
  const Spec& raw = curAmbs().at(bsd_pair.spec_id);
  const Spec& ref = curAmbs().at(bsd_pair.spec_id_ref);

  auto has_slip = [&](const Spec& spec) {
    auto sat_it = measurements.satellites.find(spec.id.gPrn());
    if (sat_it == measurements.satellites.end()) return true;
    const Satellite& satellite = sat_it->second;
    for (const auto& obs : satellite.observations) {
      if (gnss_common::getPhaseID(satellite.getSystem(), obs.first) ==
          spec.id.gPhaseId()) {
        return obs.second.slip;
      }
    }
    return true;
  };

  return has_slip(raw) || has_slip(ref);
}

// Validate PPP fixed constraints before accepting them
bool AmbiguityResolution::validatePppFixedSolution(
  const BackendId& epoch_id, const std::string& source)
{
  double range_cost_before = computeRangeCost(epoch_id);

  // Get position before solve - handle both PositionParameterBlock (3D, ECEF)
  // and PoseParameterBlock (7D, ENU position + quaternion)
  auto param_block = graph_->parameterBlockPtr(epoch_id.asInteger());
  CHECK(param_block != nullptr);
  const bool is_pose_block = (param_block->dimension() == 7);

  Eigen::Vector3d position_before;
  if (!is_pose_block) {
    std::shared_ptr<PositionParameterBlock> position_block =
      std::static_pointer_cast<PositionParameterBlock>(param_block);
    position_before = position_block->estimate();
  } else {
    // PoseParameterBlock: first 3 elements are position in ENU
    Eigen::Vector3d pos_enu(param_block->parameters()[0],
                            param_block->parameters()[1],
                            param_block->parameters()[2]);
    if (coordinate_) {
      position_before = coordinate_->convert(pos_enu, GeoType::ENU, GeoType::ECEF);
    } else {
      position_before = pos_enu;
    }
  }

  graph_->options.max_num_iterations = 1;
  graph_->options.logging_type = ceres::LoggingType::SILENT;
  graph_->options.minimizer_progress_to_stdout = false;
  graph_->solve();
  double range_cost = computeRangeCost(epoch_id);

  Eigen::Vector3d position_after;
  if (!is_pose_block) {
    std::shared_ptr<PositionParameterBlock> position_block =
      std::static_pointer_cast<PositionParameterBlock>(param_block);
    position_after = position_block->estimate();
  } else {
    Eigen::Vector3d pos_enu(param_block->parameters()[0],
                            param_block->parameters()[1],
                            param_block->parameters()[2]);
    if (coordinate_) {
      position_after = coordinate_->convert(pos_enu, GeoType::ENU, GeoType::ECEF);
    } else {
      position_after = pos_enu;
    }
  }

  double position_jump_horizontal = 0.0;
  double position_jump_vertical = 0.0;
  if (!checkZero(position_before) && !checkZero(position_after)) {
    double pos[3], dxyz[3], denu[3];
    ecef2pos(position_before.data(), pos);
    Eigen::Vector3d dposition = position_after - position_before;
    dxyz[0] = dposition.x();
    dxyz[1] = dposition.y();
    dxyz[2] = dposition.z();
    ecef2enu(pos, dxyz, denu);
    position_jump_horizontal = sqrt(denu[0] * denu[0] + denu[1] * denu[1]);
    position_jump_vertical = fabs(denu[2]);
  }
  if ((options_.max_position_jump_ar_horizontal > 0.0 &&
       position_jump_horizontal > options_.max_position_jump_ar_horizontal) ||
      (options_.max_position_jump_ar_vertical > 0.0 &&
       position_jump_vertical > options_.max_position_jump_ar_vertical)) {
    LOG(INFO) << "Invalid PPP AR " << source
              << ": Position jump is too large. horizontal="
              << std::fixed << std::setprecision(3)
              << position_jump_horizontal << " m, vertical="
              << position_jump_vertical << " m, thresholds="
              << options_.max_position_jump_ar_horizontal << "/"
              << options_.max_position_jump_ar_vertical << " m.";
    setGraphParameters(full_parameters_store_);
    eraseAmbiguityResidualBlocks(curAmbLanePairs());
    refreshCurrentAmbiguityValuesFromGraph();
    return false;
  }

  const double tolerance = 0.01;
  if (range_cost > range_cost_before + tolerance) {
    LOG(INFO) << "Invalid PPP AR " << source
              << ": Total cost increases from "
              << std::scientific << std::setprecision(3)
              << range_cost_before << " to " << range_cost << ".";
    setGraphParameters(full_parameters_store_);
    eraseAmbiguityResidualBlocks(curAmbLanePairs());
    refreshCurrentAmbiguityValuesFromGraph();
    return false;
  }

  refreshCurrentAmbiguityValuesFromGraph();
  return true;
}

// Store accepted fresh PPP NL fixed integers for short hold
void AmbiguityResolution::updateHeldPppNlFixesFromAccepted()
{
  if (!options_.use_ppp_nl_fix_hold) return;

  int num_held = 0;
  for (const auto& lane_pair : curAmbLanePairs()) {
    if (lane_pair.laneType() != LaneType::NL) continue;
    if (!lane_pair.is_fixed || lane_pair.residual_id == nullptr) continue;

    const std::string key = pppNlLaneKey(lane_pair);
    HeldNlFix& held = held_ppp_nl_fixes_[key];
    held.integer = lane_pair.fixed_integer;
    held.num_accepted++;
    held.num_failed = 0;
    held.last_epoch_sequence = current_ppp_ar_epoch_sequence_;
    num_held++;
  }
  LOG(INFO) << "PPP AR NL hold update: held_pairs=" << num_held;
}

// Apply previously accepted PPP NL integers when fresh AR drops out
bool AmbiguityResolution::tryApplyHeldPppNlFixes(
  const BackendId& epoch_id, const GnssMeasurement& measurements)
{
  if (!options_.use_ppp_nl_fix_hold) return false;
  if (held_ppp_nl_fixes_.empty()) return false;

  refreshCurrentAmbiguityValuesFromGraph();

  std::vector<size_t> held_indexes;
  std::vector<long long> held_integers;
  std::vector<std::string> held_keys;
  for (size_t i = 0; i < curAmbLanePairs().size(); i++) {
    LanePair& lane_pair = curAmbLanePairs()[i];
    if (lane_pair.laneType() != LaneType::NL) continue;

    const std::string key = pppNlLaneKey(lane_pair);
    auto held_it = held_ppp_nl_fixes_.find(key);
    if (held_it == held_ppp_nl_fixes_.end()) continue;

    HeldNlFix& held = held_it->second;
    if (options_.max_num_ppp_nl_fix_hold_epoch_gap >= 0 &&
        current_ppp_ar_epoch_sequence_ >
        held.last_epoch_sequence +
        static_cast<uint64_t>(options_.max_num_ppp_nl_fix_hold_epoch_gap)) {
      held_ppp_nl_fixes_.erase(held_it);
      continue;
    }
    if (pppNlLaneHasSlip(lane_pair, measurements)) {
      held_ppp_nl_fixes_.erase(held_it);
      continue;
    }

    const double float_cycles = lane_pair.value / lane_pair.wavelength;
    const double residual_cycles =
      fabs(float_cycles - static_cast<double>(held.integer));
    if (options_.max_ppp_nl_fix_hold_residual_cycles > 0.0 &&
        residual_cycles > options_.max_ppp_nl_fix_hold_residual_cycles) {
      LOG(INFO) << "PPP AR NL hold skip: key=" << key
                << ", float=" << float_cycles
                << ", held=" << held.integer
                << ", residual_cycles=" << residual_cycles
                << ", threshold="
                << options_.max_ppp_nl_fix_hold_residual_cycles;
      continue;
    }

    held_indexes.push_back(i);
    held_integers.push_back(held.integer);
    held_keys.push_back(key);
  }

  if (static_cast<int>(held_indexes.size()) <
      options_.min_num_satellite_pairs_fixation) {
    return false;
  }

  const double information = 1.0e6;
  for (size_t i = 0; i < held_indexes.size(); i++) {
    LanePair& lane_pair = curAmbLanePairs()[held_indexes[i]];
    CHECK(lane_pair.laneType() == LaneType::NL);
    CHECK(lane_pair.residual_id == nullptr);

    std::vector<double> coefficients;
    const size_t id_higher_raw =
      curAmbPairs()[lane_pair.bsd_pair_id_higher].spec_id;
    const size_t id_higher_ref =
      curAmbPairs()[lane_pair.bsd_pair_id_higher].spec_id_ref;
    const double wave = lane_pair.wavelength;
    coefficients.push_back(1.0 / wave);
    coefficients.push_back(-1.0 / wave);

    std::shared_ptr<AmbiguityError2Coef> ambiguity_error =
      std::make_shared<AmbiguityError2Coef>(
        static_cast<double>(held_integers[i]), information, coefficients);
    lane_pair.residual_id = graph_->addResidualBlock(
      ambiguity_error,
      cauchy_loss_function_ ? cauchy_loss_function_.get() : nullptr,
      graph_->parameterBlockPtr(curAmbs()[id_higher_raw].id.asInteger()),
      graph_->parameterBlockPtr(curAmbs()[id_higher_ref].id.asInteger()));
    lane_pair.is_fixed = true;
    lane_pair.fixed_integer = held_integers[i];
    curAmbPairs()[lane_pair.bsd_pair_id_higher].is_fixed = true;
  }

  if (!validatePppFixedSolution(epoch_id, "hold")) {
    for (const auto& key : held_keys) held_ppp_nl_fixes_.erase(key);
    LOG(INFO) << "PPP AR NL hold rejected: removed_held_pairs="
              << held_keys.size();
    return false;
  }

  for (const auto& key : held_keys) {
    HeldNlFix& held = held_ppp_nl_fixes_[key];
    held.num_accepted++;
    held.num_failed = 0;
    held.last_epoch_sequence = current_ppp_ar_epoch_sequence_;
  }
  LOG(INFO) << "PPP AR NL hold accepted: fixed_pairs=" << held_keys.size();
  return true;
}

// Erase all added ambiguity residual blocks 
void AmbiguityResolution::eraseAmbiguityResidualBlocks(
  std::vector<LanePair>& lane_pairs)
{
  for (auto& lane_pair : lane_pairs) {
    if (lane_pair.residual_id) {
      graph_->removeResidualBlock(lane_pair.residual_id);
      lane_pair.residual_id = nullptr;
    }
    lane_pair.is_fixed = false;
    lane_pair.fixed_integer = 0;
  }
}

// Set graph parameter values
void AmbiguityResolution::setGraphParameters(std::vector<Parameter>& parameters)
{
  for (auto& parameter : parameters) {
    Eigen::Map<Eigen::VectorXd> value(
      parameter.handle->parameters(), parameter.size);
    value = parameter.value;
  }
}

// Refresh current ambiguity, BSD, and lane values from graph parameters
void AmbiguityResolution::refreshCurrentAmbiguityValuesFromGraph()
{
  for (auto& ambiguity : curAmbs()) {
    CHECK(ambiguity.parameter_block != nullptr);
    ambiguity.value = *ambiguity.parameter_block->parameters();
  }
  for (auto& pair : curAmbPairs()) pair.update(curAmbs());
  for (auto& lane_pair : curAmbLanePairs()) lane_pair.update(curAmbPairs());
}

// Compute pseudorange and phasernage total cost in current epoch
double AmbiguityResolution::computeRangeCost(const BackendId& epoch_id)
{
  Graph::ResidualBlockCollection residual_blocks = 
    graph_->residuals(epoch_id.asInteger());
  double total_cost_square = 0.0;
  std::vector<double> residuals;
  for (size_t i = 0; i < residual_blocks.size(); i++) {
    auto& residual_block = residual_blocks[i];
    ErrorType type = residual_block.error_interface_ptr->typeInfo();
    if (!(type == ErrorType::kPseudorangeError || 
          type == ErrorType::kPseudorangeErrorSD || 
          type == ErrorType::kPseudorangeErrorDD ||
          type == ErrorType::kPhaserangeError || 
          type == ErrorType::kPhaserangeErrorSD || 
          type == ErrorType::kPhaserangeErrorDD)) continue;
    double residual[1];
    graph_->problem()->EvaluateResidualBlock(residual_block.residual_block_id, 
      false, nullptr, residual, nullptr);
    residuals.push_back(*residual);
    total_cost_square += square(*residual);
  }
  return sqrt(total_cost_square);
}

}
