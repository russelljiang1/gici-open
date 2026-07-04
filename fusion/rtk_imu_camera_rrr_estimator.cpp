/**
* @Function: RTK/IMU/Camera tightly couple estimator (GNSS raw (RTK formula) + IMU raw + camera raw)
*
* @Author  : Cheng Chi
* @Email   : chichengcn@sjtu.edu.cn
*
* Copyright (C) 2023 by Cheng Chi, All rights reserved.
**/
#include "gici/fusion/rtk_imu_camera_rrr_estimator.h"

#include "gici/gnss/position_error.h"

namespace gici {

// ===== 全局计时统计变量 =====
static double cumulative_constraint_time_ms = 0.0;
static double cumulative_optimize_time_ms = 0.0; 
static double cumulative_outlier_time_ms = 0.0;
static double cumulative_margin_time_ms = 0.0;
static double cumulative_visual_total_ms = 0.0;
static double cumulative_total_time_ms = 0.0;
static int constraint_count = 0;
static int optimize_count = 0;
static int outlier_count = 0;
static int margin_count = 0;

// The default constructor
RtkImuCameraRrrEstimator::RtkImuCameraRrrEstimator(
               const RtkImuCameraRrrEstimatorOptions& options, 
               const GnssImuInitializerOptions& init_options, 
               const RtkEstimatorOptions rtk_options,
               const GnssEstimatorBaseOptions& gnss_base_options, 
               const GnssLooseEstimatorBaseOptions& gnss_loose_base_options, 
               const VisualEstimatorBaseOptions& visual_base_options,
               const ImuEstimatorBaseOptions& imu_base_options,
               const EstimatorBaseOptions& base_options,
               const AmbiguityResolutionOptions& ambiguity_options) :
  rrr_options_(options), rtk_options_(rtk_options),
  GnssEstimatorBase(gnss_base_options, base_options),
  VisualEstimatorBase(visual_base_options, base_options),
  ImuEstimatorBase(imu_base_options, base_options),
  EstimatorBase(base_options)
{
  
  type_ = EstimatorType::RtkImuCameraRrr;
  is_use_phase_ = true;
  states_.push_back(State());
  gnss_measurement_pairs_.push_back(
    std::make_pair(GnssMeasurement(), GnssMeasurement()));
  frame_bundles_.push_back(nullptr);
  num_satellites_ = 0;

  // Initialization control
  initializer_sub_estimator_.reset(new RtkEstimator(
    rtk_options, gnss_base_options, base_options, ambiguity_options));
  gnss_imu_initializer_.reset(new GnssImuInitializer(
    init_options, gnss_loose_base_options, imu_base_options, 
    base_options, graph_, initializer_sub_estimator_));

  // Ambiguity resolution
  ambiguity_resolution_.reset(new AmbiguityResolution(ambiguity_options, graph_));
  
  // RTK estimator used for ambiguity covariance estimation
  RtkEstimatorOptions sub_rtk_options = rtk_options;
  EstimatorBaseOptions sub_base_options = base_options;
  GnssEstimatorBaseOptions sub_gnss_base_options = gnss_base_options;
  sub_rtk_options.use_ambiguity_resolution = false;
  sub_rtk_options.max_window_length = 2;
  sub_base_options.verbose_output = false;
  sub_gnss_base_options.use_outlier_rejection = false;
  ambiguity_covariance_estimator_.reset(new RtkEstimator(sub_rtk_options, 
    sub_gnss_base_options, sub_base_options, ambiguity_options));

}

// The default destructor
RtkImuCameraRrrEstimator::~RtkImuCameraRrrEstimator()
{}

// Add measurement
bool RtkImuCameraRrrEstimator::addMeasurement(const EstimatorDataCluster& measurement)
{
  // GNSS/IMU initialization
  if (coordinate_ == nullptr || !gravity_setted_) return false;
  if (!gnss_imu_initializer_->finished()) {
    if (gnss_imu_initializer_->getCoordinate() == nullptr) {
      gnss_imu_initializer_->setCoordinate(coordinate_);
      initializer_sub_estimator_->setCoordinate(coordinate_);
      gnss_imu_initializer_->setGravity(imu_base_options_.imu_parameters.g);
    }
    if (gnss_imu_initializer_->addMeasurement(measurement)) {
      gnss_imu_initializer_->estimate();
      // set result to estimator
      setInitializationResult(gnss_imu_initializer_);
    }
    return false;
  }

  // Add IMU
  if (measurement.imu && measurement.imu_role == ImuRole::Major) {
    addImuMeasurement(*measurement.imu);
  }

  // Add GNSS
  if (measurement.gnss) {
    // feed to covariance estimator
    if (!ambiguity_covariance_coordinate_setted_ && coordinate_) {
      ambiguity_covariance_estimator_->setCoordinate(coordinate_);
      ambiguity_covariance_coordinate_setted_ = true;
    }
    if (coordinate_ && ambiguity_covariance_estimator_->addMeasurement(measurement)) {
      ambiguity_covariance_estimator_->estimate();
    }
    // feed to local
    GnssMeasurement rov, ref;
    meausrement_align_.add(measurement);
    if (meausrement_align_.get(rtk_options_.max_age, rov, ref)) {
      return addGnssMeasurementAndState(rov, ref);
    }
  }

  // Add images
  if (measurement.frame_bundle) {
    if (!visual_initialized_) return visualInitialization(measurement.frame_bundle);
    return addImageMeasurementAndState(measurement.frame_bundle);
  }

  return false;
}

// Add GNSS measurements and state
bool RtkImuCameraRrrEstimator::addGnssMeasurementAndState(
    const GnssMeasurement& measurement_rov, 
    const GnssMeasurement& measurement_ref)
{
  // Get prior states
  Eigen::Vector3d position_prior = coordinate_->convert(
    getPoseEstimate().getPosition(), GeoType::ENU, GeoType::ECEF);

  // Set to local measurement handle
  curGnssRov() = measurement_rov;
  curGnssRov().position = position_prior;
  curGnssRef() = measurement_ref;

  // Erase duplicated phases, arrange to one observation per phase
  gnss_common::rearrangePhasesAndCodes(curGnssRov());
  gnss_common::rearrangePhasesAndCodes(curGnssRef());

  // Form double difference pair
  std::map<char, std::string> system_to_base_prn;
  GnssMeasurementDDIndexPairs phase_index_pairs = gnss_common::formPhaserangeDDPair(
    curGnssRov(), curGnssRef(), system_to_base_prn, gnss_base_options_.common);
  GnssMeasurementDDIndexPairs code_index_pairs = gnss_common::formPseudorangeDDPair(
    curGnssRov(), curGnssRef(), system_to_base_prn, gnss_base_options_.common);

  // Cycle-slip detection
  if (!isFirstEpoch()) {
    cycleSlipDetectionSD(lastGnssRov(), lastGnssRef(), 
      curGnssRov(), curGnssRef(), gnss_base_options_.common);
  }

  // Add parameter blocks
  double timestamp = curGnssRov().timestamp;
  // pose and speed and bias block
  const int32_t bundle_id = curGnssRov().id;
  BackendId pose_id = createGnssPoseId(bundle_id);
  size_t index = insertImuState(timestamp, pose_id);
  states_[index].status = GnssSolutionStatus::Single;
  latest_state_index_ = index;
  // GNSS extrinsics, it should be added at initialization step
  CHECK(gnss_extrinsics_id_.valid());
  // ambiguity blocks
  addSdAmbiguityParameterBlocks(curGnssRov(), 
    curGnssRef(), phase_index_pairs, curGnssRov().id, curAmbiguityState());
  // frequency block
  int num_valid_doppler_system = 0;
  addFrequencyParameterBlocks(curGnssRov(), curGnssRov().id, num_valid_doppler_system);

  // Add pseudorange residual blocks
  int num_valid_satellite = 0;
  addDdPseudorangeResidualBlocks(curGnssRov(), 
    curGnssRef(), code_index_pairs, states_[index], num_valid_satellite);
  
  // We do not need to check if the number of satellites is sufficient in tightly fusion.
  if (!checkSufficientSatellite(num_valid_satellite, 0)) {
    // do nothing
  }
  num_satellites_ = num_valid_satellite;

  // No satellite
  if (num_satellites_ == 0) {
    // erase parameters in current state
    eraseFrequencyParameterBlocks(states_[index]);
    eraseAmbiguityParameterBlocks(curAmbiguityState());
    eraseImuState(states_[index]);
    return false;
  }

  // Add phaserange residual blocks
  addDdPhaserangeResidualBlocks(
    curGnssRov(), curGnssRef(), phase_index_pairs, states_[index]);

  // Add doppler residual blocks
  addDopplerResidualBlocks(curGnssRov(), states_[index], num_valid_satellite, 
    false, getImuMeasurementNear(timestamp).angular_velocity);

  // Add relative errors
  if (lastGnssState().valid()) {  // maybe invalid here because of long term GNSS absent
    // frequency
    addRelativeFrequencyResidualBlock(lastGnssState(), states_[index]);
    // ambiguity
    addRelativeAmbiguityResidualBlock(
      lastGnssRov(), curGnssRov(), lastAmbiguityState(), curAmbiguityState());
  }

  // ZUPT
  addZUPTResidualBlock(states_[index]);

  // Car motion
  if (imu_base_options_.car_motion) {
    // heading measurement constraint
    addHMCResidualBlock(states_[index]);
    // non-holonomic constraint
    addNHCResidualBlock(states_[index]);
  }

  // Compute DOP
  updateGdop(curGnssRov(), code_index_pairs);

  return true;
}

// Add image measurements and state
bool RtkImuCameraRrrEstimator::addImageMeasurementAndState(
  const FrameBundlePtr& frame_bundle, const SpeedAndBias& speed_and_bias)
{
  // If initialized, we are supposed not at the first epoch
  CHECK(!isFirstEpoch());

  // Check frequency
  if (!frame_bundle->isKeyframe() && frame_bundles_.size() > 1) {
    const double t_last = lastFrameBundle()->getMinTimestampSeconds();
    const double t_cur = frame_bundle->getMinTimestampSeconds();
    if (t_cur - t_last < 1.0 / visual_base_options_.max_frequency) return false;
  }

  // Set to local measurement handle
  curFrameBundle() = frame_bundle;

  // Add parameter blocks
  double timestamp = curFrame()->getTimestampSec();
  // pose and speed and bias block
  const int32_t bundle_id = curFrame()->bundleId();
  BackendId pose_id = createNFrameId(bundle_id);
  size_t index;
  if (speed_and_bias != SpeedAndBias::Zero()) {
    index = insertImuState(timestamp, pose_id, 
      curFrame()->T_world_imu(), speed_and_bias, true);
  }
  else {
    index = insertImuState(timestamp, pose_id);
  }
  states_[index].status = latestGnssState().status;
  states_[index].is_keyframe = curFrame()->isKeyframe();
  latest_state_index_ = index;
  // camera extrinsics
  if (!camera_extrinsics_id_.valid()) {
    camera_extrinsics_id_ = 
      addCameraExtrinsicsParameterBlock(bundle_id, curFrame()->T_imu_cam());
    addCameraExtrinsicsResidualBlock(camera_extrinsics_id_, curFrame()->T_imu_cam(), 
      visual_base_options_.camera_extrinsics_initial_std.head<3>(), 
      visual_base_options_.camera_extrinsics_initial_std.tail<3>() * D2R);
  }

  // ===== 约束构建计时 =====
  auto constraint_start = std::chrono::high_resolution_clock::now();

  // 所有帧：更新特征轨迹
  updateFeatureTracksWithNewFrame(states_[index].id, curFrameBundle());

  // 关键帧：同时构建两种约束
  if (visual_initialized_ && curFrame()->isKeyframe()) {
    curFrame()->set_T_w_imu(getPoseEstimate(states_[index]));
    
    // 构建DPO约束
    buildDPOConstraints(states_[index], curFrameBundle());    
    
    /*
    // 构建landmark约束
    feature_handler_->initializeLandmarks(curFrame());
    addLandmarkParameterBlocksWithResiduals(curFrame());
    addReprojectionErrorResidualBlocks(states_[index], curFrame());
    */
  }

  auto constraint_end = std::chrono::high_resolution_clock::now();
  auto constraint_duration = std::chrono::duration_cast<std::chrono::milliseconds>(constraint_end - constraint_start);

  // ===== 约束构建加合统计 =====
  cumulative_constraint_time_ms += constraint_duration.count();
  cumulative_visual_total_ms += constraint_duration.count();
  cumulative_total_time_ms += constraint_duration.count();
  constraint_count++;
  
  LOG(INFO) << "=== Constraint Building: took " << constraint_duration.count() << "ms ===";
  LOG(INFO) << "=== Constraint Cumulative: Total " << cumulative_constraint_time_ms << "ms"
            << " (" << cumulative_constraint_time_ms/1000.0 << "s)";
  
  // 写入CSV
  static std::ofstream constraint_csv("constraint_timing.csv", std::ios::app);
  if (constraint_csv.is_open()) {
      constraint_csv << constraint_count << "," << constraint_duration.count() 
                    << "," << cumulative_constraint_time_ms/1000.0 
                    << "," << "DPO" << std::endl;
  }

  return true;
}

// Visual initialization
bool RtkImuCameraRrrEstimator::visualInitialization(const FrameBundlePtr& frame_bundle)
{
  // store poses
  do_not_remove_imu_measurements_ = true;
  Solution solution;
  solution.timestamp = getTimestamp();
  solution.pose = getPoseEstimate();
  solution.speed_and_bias = getSpeedAndBiasEstimate();
  init_solution_store_.push_back(solution);
  
  // store keyframes
  if (frame_bundle->isKeyframe()) {
    init_keyframes_.push_back(frame_bundle);
    while (init_keyframes_.size() > 2) init_keyframes_.pop_front();
  } 
  if (init_keyframes_.size() < 2) return false;

  // check if GNSS/IMU estimator fully converged
  double std_yaw = sqrt(computeAndGetCovariance(lastState())(5, 5));
  if (std_yaw > rrr_options_.min_yaw_std_init_visual * D2R) return false;

  // set poses
  std::vector<SpeedAndBias> speed_and_biases;
  for (auto& frame_bundle : init_keyframes_) {
    bool found = false;
    double timestamp = frame_bundle->getMinTimestampSeconds();
    for (size_t i = 1; i < init_solution_store_.size(); i++) {
      double dt1 = init_solution_store_[i].timestamp - timestamp;
      double dt2 = init_solution_store_[i - 1].timestamp - timestamp;
      if ((dt1 >= 0 && dt2 <= 0) || 
          (i == init_solution_store_.size() - 1) && dt1 < 0 && fabs(dt1) < 2.0) {
        size_t idx = (dt1 >= 0 && dt2 <= 0) ? (i - 1) : i;
        Transformation T_WS = init_solution_store_[idx].pose;
        SpeedAndBias speed_and_bias = init_solution_store_[idx].speed_and_bias;
        imuIntegration(init_solution_store_[idx].timestamp, 
          timestamp, T_WS, speed_and_bias);
        speed_and_biases.push_back(speed_and_bias);
        frame_bundle->set_T_W_B(T_WS);
        found = true;
        break;
      }
    }
    CHECK(found);
  }
  // initialize landmarks
  feature_handler_->initializeLandmarks(init_keyframes_.back()->at(0));
  feature_handler_->setGlobalScaleInitialized();
  // add two keyframes
  CHECK(addImageMeasurementAndState(init_keyframes_.front(), speed_and_biases.front()));
  CHECK(addImageMeasurementAndState(init_keyframes_.back(), speed_and_biases.back()));

  // set flag
  visual_initialized_ = true;
  do_not_remove_imu_measurements_ = false;
  can_compute_covariance_ = true;
  return true;
}

// Solve current graph
bool RtkImuCameraRrrEstimator::estimate()
{
  // ===== 开始计时 =====
  auto start_time = std::chrono::high_resolution_clock::now();

  status_ = EstimatorStatus::Converged;

  // ===== optimize函数单独计时 =====
  auto optimize_start = std::chrono::high_resolution_clock::now();

  // Optimize
  optimize();

  auto optimize_end = std::chrono::high_resolution_clock::now();
  auto optimize_duration = std::chrono::duration_cast<std::chrono::milliseconds>(optimize_end - optimize_start);
  
  cumulative_optimize_time_ms += optimize_duration.count();
  optimize_count++;
  
  LOG(INFO) << "=== Optimize: took " << optimize_duration.count() << "ms ===";
  LOG(INFO) << "=== Optimize Cumulative: Total " << cumulative_optimize_time_ms << "ms"
            << " (" << cumulative_optimize_time_ms/1000.0 << "s)";

  // Get current sensor type
  IdType new_state_type = states_[latest_state_index_].id.type();

  // ===== 声明outlier_duration变量 =====
  std::chrono::milliseconds outlier_duration(0);

  // GNSS processes
  double gdop = 0.0;
  if (new_state_type == IdType::gPose)
  {
    // Reject GNSS outliers
    size_t n_pseudorange = numPseudorangeError(states_[latest_state_index_]);
    size_t n_phaserange = numPhaserangeError(states_[latest_state_index_]);
    size_t n_doppler = numDopplerError(states_[latest_state_index_]);
    if (gnss_base_options_.use_outlier_rejection) {
      rejectPseudorangeOutlier(states_[latest_state_index_], curAmbiguityState());
      rejectDopplerOutlier(states_[latest_state_index_]);
      rejectPhaserangeOutlier(states_[latest_state_index_], curAmbiguityState());
    }

    // Check if we rejected too many GNSS residuals
    double ratio_pseudorange = n_pseudorange == 0.0 ? 0.0 : 1.0 - 
      getDivide(numPseudorangeError(states_[latest_state_index_]), n_pseudorange);
    double ratio_phaserange = n_phaserange == 0.0 ? 0.0 : 1.0 - 
      getDivide(numPhaserangeError(states_[latest_state_index_]), n_phaserange);
    double ratio_doppler = n_doppler == 0.0 ? 0.0 : 1.0 - 
      getDivide(numDopplerError(states_[latest_state_index_]), n_doppler);
    const double thr = gnss_base_options_.diverge_max_reject_ratio;
    if (isGnssGoodObservation() && 
        (ratio_pseudorange > thr || ratio_phaserange > thr || ratio_doppler > thr)) {
      num_cotinuous_reject_gnss_++;
    }
    else num_cotinuous_reject_gnss_ = 0;
    if (num_cotinuous_reject_gnss_ > 
        gnss_base_options_.diverge_min_num_continuous_reject) {
      LOG(WARNING) << "Estimator diverge: Too many GNSS outliers rejected!";
      status_ = EstimatorStatus::Diverged;
      num_cotinuous_reject_gnss_ = 0;
    }

    // Ambiguity resolution
    for (size_t i = latest_state_index_; i < states_.size(); i++) {
      states_[i].status = GnssSolutionStatus::Float;
    }
    if (rtk_options_.use_ambiguity_resolution) {
      // get covariance of ambiguities
      Eigen::MatrixXd ambiguity_covariance;
      if (estimateAmbiguityCovariance(
        states_[latest_state_index_], ambiguity_covariance))
      {
        // solve
        AmbiguityResolution::Result ret = ambiguity_resolution_->solveRtk(
          states_[latest_state_index_].id, curAmbiguityState().ids, 
          ambiguity_covariance, gnss_measurement_pairs_.back());
        if (ret == AmbiguityResolution::Result::NlFix) {
          for (size_t i = latest_state_index_; i < states_.size(); i++) {
            states_[i].status = GnssSolutionStatus::Fixed;
          }
        }
      }
    }

    // Check if we continuously cannot fix ambiguity, while we have good observations
    if (rtk_options_.use_ambiguity_resolution) {
      const double thr = gnss_base_options_.good_observation_max_reject_ratio;
      if (isGnssGoodObservation() && ratio_pseudorange < thr && 
          ratio_phaserange < thr && ratio_doppler < thr) {
        if (curState().status != GnssSolutionStatus::Fixed) num_continuous_unfix_++;
        else num_continuous_unfix_ = 0; 
      }
      else num_continuous_unfix_ = 0;
      if (num_continuous_unfix_ > 
          gnss_base_options_.reset_ambiguity_min_num_continuous_unfix) {
        LOG(INFO) << "Continuously unfix under good observations. Clear current ambiguities.";
        resetAmbiguityEstimation();
        ambiguity_covariance_estimator_->resetAmbiguityEstimation();
        num_continuous_unfix_ = 0;
      }
    }
  }


  // === 视觉处理部分替换 ===
  // Image processes
  if (new_state_type == IdType::cPose) {
    // update landmarks to frontend
    // updateLandmarks();
    // update states to frontend
    // updateFrameStateToFrontend(states_[latest_state_index_], curFrame());
    auto outlier_start = std::chrono::high_resolution_clock::now();
    // DPO outlier rejection
    size_t n_dpo_constraints = numDPOConstraints(states_[latest_state_index_]);
    rejectDPOConstraintOutlier(states_[latest_state_index_]);

    // Landmark outlier rejection
    // size_t n_reprojection = numReprojectionError(curFrame());
    // rejectReprojectionErrorOutlier(curFrame());
    
    // 检查rejection ratio
    double ratio_dpo = n_dpo_constraints == 0.0 ? 0.0 : 1.0 - 
      getDivide(numDPOConstraints(states_[latest_state_index_]), n_dpo_constraints);

    // double ratio_reprojection = n_reprojection == 0.0 ? 0.0 : 1.0 - 
    //   getDivide(numReprojectionError(curFrame()), n_reprojection);
    auto outlier_end = std::chrono::high_resolution_clock::now();
    auto outlier_duration = std::chrono::duration_cast<std::chrono::milliseconds>(outlier_end - outlier_start);
    
    // ===== outlier rejection加合统计 =====
    cumulative_outlier_time_ms += outlier_duration.count();
    cumulative_visual_total_ms += outlier_duration.count();
    outlier_count++;
    
    LOG(INFO) << "=== Visual Outlier Rejection: took " << outlier_duration.count() << "ms ===";
    LOG(INFO) << "=== Outlier Cumulative: Total " << cumulative_outlier_time_ms << "ms"
              << " (" << cumulative_outlier_time_ms/1000.0 << "s)";
              
    if (ratio_dpo > visual_base_options_.diverge_max_reject_ratio) {
      num_cotinuous_reject_visual_++;
    } else {
      num_cotinuous_reject_visual_ = 0;
    }

    if (num_cotinuous_reject_visual_ > 
        visual_base_options_.diverge_min_num_continuous_reject) {
      LOG(WARNING) << "Estimator diverge: Too many DPO outliers rejected!";
      status_ = EstimatorStatus::Diverged;
      num_cotinuous_reject_visual_ = 0;
    }
  }


  // Log information
  State& new_state = states_[latest_state_index_];
  if (base_options_.verbose_output) {
    int dif_distance = static_cast<int>(round(
      (curGnssRov().position - curGnssRef().position).norm() / 1.0e3));
    LOG(INFO) << estimatorTypeToString(type_) << ": " 
      << "Iterations: " << graph_->summary.iterations.size() << ", "
      << std::scientific << std::setprecision(3) 
      << "Initial cost: " << graph_->summary.initial_cost << ", "
      << "Final cost: " << graph_->summary.final_cost
      << ", Sensor type: " << std::setw(1) << 
        static_cast<int>(BackendId::sensorType(new_state.id.type()))
      << ", Sat number: " << std::setw(2) << num_satellites_
      << ", GDOP: " << std::setprecision(1) << std::fixed << gdop_
      << ", Fix status: " << std::setw(1) << static_cast<int>(new_state.status);
  }

  // Apply marginalization
  marginalization(new_state_type);

  // Shift memory for states and measurements
  if (new_state_type == IdType::gPose) {
    gnss_measurement_pairs_.push_back(
      std::make_pair(GnssMeasurement(), GnssMeasurement()));
    ambiguity_states_.push_back(AmbiguityState());
  }
  if (new_state_type == IdType::cPose) frame_bundles_.push_back(nullptr);
  states_.push_back(State());
  while (!visual_initialized_ && 
         states_.size() > rrr_options_.max_gnss_window_length_minor) {
    states_.pop_front();
    ambiguity_states_.pop_front();
    gnss_measurement_pairs_.pop_front();
  }
  // only keep frame measurement data for two epochs
  while (frame_bundles_.size() > 2) frame_bundles_.pop_front();

  // ===== 计时结束 =====
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
  
  // ===== 添加累积时间记录 =====
  static int total_frame_count = 0;
  
  // ===== 添加累积时间记录 =====
  cumulative_total_time_ms += duration.count();
  total_frame_count++;
  
  LOG(INFO) << "=== Timing: estimate() took " << duration.count() << "ms ===";
  LOG(INFO) << "=== Total Cumulative: Frame " << total_frame_count 
            << ", Total time: " << cumulative_total_time_ms << "ms"
            << " (" << cumulative_total_time_ms/1000.0 << "s)";

  LOG(INFO) << "=== Timing: estimate() took " << duration.count() << "ms ===";
  
  // ===== 所有计时的汇总输出 =====
  LOG(INFO) << "=== TIMING BREAKDOWN SUMMARY ===";
  LOG(INFO) << "Constraints: " << cumulative_constraint_time_ms << "ms (" 
            << (constraint_count > 0 ? cumulative_constraint_time_ms/constraint_count : 0) << "ms avg)";
  LOG(INFO) << "Optimize: " << cumulative_optimize_time_ms << "ms (" 
            << (optimize_count > 0 ? cumulative_optimize_time_ms/optimize_count : 0) << "ms avg)";
  LOG(INFO) << "Outlier: " << cumulative_outlier_time_ms << "ms (" 
            << (outlier_count > 0 ? cumulative_outlier_time_ms/outlier_count : 0) << "ms avg)";
  LOG(INFO) << "Marginalization: " << cumulative_margin_time_ms << "ms (" 
            << (margin_count > 0 ? cumulative_margin_time_ms/margin_count : 0) << "ms avg)";
  LOG(INFO) << "VISUAL TOTAL: " << cumulative_visual_total_ms << "ms (" 
            << cumulative_visual_total_ms/1000.0 << "s)";
  LOG(INFO) << "SYSTEM TOTAL: " << cumulative_total_time_ms << "ms (" 
            << cumulative_total_time_ms/1000.0 << "s)";
  LOG(INFO) << "================================";

  // ===== CSV输出所有数据 =====
  static std::ofstream detail_csv("detailed_timing.csv", std::ios::app);
  if (detail_csv.is_open()) {
      detail_csv << total_frame_count << "," 
                << cumulative_constraint_time_ms << ","
                << cumulative_optimize_time_ms << ","
                << cumulative_outlier_time_ms << ","
                << cumulative_margin_time_ms << ","
                << cumulative_visual_total_ms << ","
                << cumulative_total_time_ms << ","
                << "Original" << std::endl;
  }

  return true;
}

// Set initializatin result
void RtkImuCameraRrrEstimator::setInitializationResult(
  const std::shared_ptr<MultisensorInitializerBase>& initializer)
{
  CHECK(initializer->finished());

  // Cast to desired initializer
  std::shared_ptr<GnssImuInitializer> gnss_imu_initializer = 
    std::static_pointer_cast<GnssImuInitializer>(initializer);
  CHECK_NOTNULL(gnss_imu_initializer);
  
  // Arrange to window length
  ImuMeasurements imu_measurements;
  std::deque<GnssSolution> gnss_measurement_temp;  // we do not use this
  gnss_imu_initializer->arrangeToEstimator(
    rrr_options_.max_gnss_window_length_minor, marginalization_error_, states_, 
    marginalization_residual_id_, gnss_extrinsics_id_, 
    gnss_measurement_temp, imu_measurements);
  for (auto it = imu_measurements.rbegin(); it != imu_measurements.rend(); it++) {
    imu_mutex_.lock();
    imu_measurements_.push_front(*it);
    imu_mutex_.unlock();
  }

  // Shift memory for states and measurements
  states_.push_back(State());
  gnss_measurement_pairs_.resize(states_.size());
  ambiguity_states_.resize(states_.size());
  frame_bundles_.push_back(nullptr);
}

// Marginalization
bool RtkImuCameraRrrEstimator::marginalization(const IdType& type)
{
  if (type == IdType::cPose) return frameMarginalization();
  else if (type == IdType::gPose) return gnssMarginalization();
  else return false;
}

// Marginalization when the new state is a frame state
bool RtkImuCameraRrrEstimator::frameMarginalization()
{
  // Check if we need marginalization
  if (isFirstEpoch()) return true;

  // Make sure that only current state can be non-keyframe
  for (int i = 0; i < states_.size() - 1; i++) {
    if (states_[i].id.type() != IdType::cPose) continue;
    if (!states_[i].is_keyframe) {

      // eraseReprojectionErrorResidualBlocks(states_[i]);
      
      // 清理DPO约束-非关键帧
      eraseDPOConstraints(states_[i]);

      eraseImuState(states_[i]);
      i--;

    }
  }

  // If current frame is a keyframe. Marginalize the oldest keyframe and corresponding 
  // IMU and GNSS states out. And sparsify GNSS states.
  if (states_[latest_state_index_].is_keyframe &&
      sizeOfKeyframeStates() > rrr_options_.max_keyframes) {

    // ===== 视觉边缘化计时 =====
    auto margin_start = std::chrono::high_resolution_clock::now();
        
    // Erase old marginalization item
    if (!eraseOldMarginalization()) return false;

    // Add marginalization items
    // marginalize oldest keyframe state
    bool reached_first_keyframe = false;
    State margin_keyframe_state;
    int n_state = 0, n_frame = 0, n_gnss = 0;
    for (auto it = states_.begin(); it != states_.end();) {
      State& state = *it;
      // reached the first keyframe
      if (state.is_keyframe) reached_first_keyframe = true;
      
      // mark marginalized keyframe state
      if (state.is_keyframe) margin_keyframe_state = state;

      // add margin blocks
      // Keyframe, we add the reprojection errors latter.
      if (state.is_keyframe) {
        addImuStateMarginBlock(state);
        addImuResidualMarginBlocks(state);
      }
      // GNSS state that not yet been marginalized by GNSS steps.
      else {
        auto it_ambiguity = ambiguityStateAt(state.timestamp);
        if (it_ambiguity != ambiguity_states_.end()) {
          addAmbiguityMarginBlocksWithResiduals(*it_ambiguity);
        }
        else {
          addGnssLooseResidualMarginBlocks(state);
        }
        addImuStateMarginBlock(state);
        addImuResidualMarginBlocks(state);
        addFrequencyMarginBlocksWithResiduals(state);
        addGnssResidualMarginBlocks(state);

        // erase ambiguity state and GNSS measurements
        if (it_ambiguity != ambiguity_states_.end()) {
          ambiguity_states_.erase(it_ambiguity);
        }
        if (gnssMeasurementPairAt(state.timestamp) != gnss_measurement_pairs_.end()) {
          gnss_measurement_pairs_.erase(gnssMeasurementPairAt(state.timestamp));
        }
      }

      // Erase state
      it = states_.erase(it);

      n_state++;

      if (reached_first_keyframe) break;
    }


    // DPO约束-关键帧边缘化
    eraseDPOConstraints(margin_keyframe_state);
    
    // Landmark边缘化
    // eraseEmptyLandmarks();
    // addLandmarkParameterMarginBlocksWithResiduals(margin_keyframe_state);

    auto margin_end = std::chrono::high_resolution_clock::now();
    auto margin_duration = std::chrono::duration_cast<std::chrono::milliseconds>(margin_end - margin_start);
    
    // ===== 边缘化加合统计 =====
    cumulative_margin_time_ms += margin_duration.count();
    cumulative_visual_total_ms += margin_duration.count();
    margin_count++;

    LOG(INFO) << "=== Visual Marginalization: took " << margin_duration.count() << "ms ===";
    LOG(INFO) << "=== Margin Cumulative: Total " << cumulative_margin_time_ms << "ms"
              << " (" << cumulative_margin_time_ms/1000.0 << "s)";
              
    // Apply marginalization and add the item into graph
    bool ret = applyMarginalization();

    return ret;
  }

  return true;
}

// Marginalization when the new state is a GNSS state
bool RtkImuCameraRrrEstimator::gnssMarginalization()
{
  // Check if we need marginalization
  if (isFirstEpoch()) return true;

  // Visual not initialzed, only handle GNSS and INS
  if (!visual_initialized_) {
    // Check if we need marginalization
    if (states_.size() < rrr_options_.max_gnss_window_length_minor) {
      return true;
    }

    // Erase old marginalization item
    if (!eraseOldMarginalization()) return false;

    // Add marginalization items
    // IMU states and residuals
    addImuStateMarginBlockWithResiduals(oldestState());
    // ambiguity
    addAmbiguityMarginBlocksWithResiduals(oldestAmbiguityState());
    // frequency
    addFrequencyMarginBlocksWithResiduals(oldestState());

    // Apply marginalization and add the item into graph
    bool ret = applyMarginalization();

    return ret;
  }
  // Visual initialized, handle all the sensors
  else {
    // Sparsify GNSS states
    sparsifyGnssStates();

    // Marginalize the GNSS states that in front of the oldest keyframe
    // We do this at both here and the visual margin step to smooth computational load
    for (auto it = states_.begin(); it != states_.end();) {
      State& state = *it;
      // reached the oldest keyframe
      if (state.is_keyframe) break;
      // check if GNSS state
      if (state.id.type() != IdType::gPose) continue;

      // Erase old marginalization item
      if (!eraseOldMarginalization()) return false;

      // margin
      auto it_ambiguity = ambiguityStateAt(state.timestamp);
      if (it_ambiguity != ambiguity_states_.end()) {
        addAmbiguityMarginBlocksWithResiduals(*it_ambiguity);
      }
      else {
        addGnssLooseResidualMarginBlocks(state);
      }
      addImuStateMarginBlock(state);
      addImuResidualMarginBlocks(state);
      addFrequencyMarginBlocksWithResiduals(state);
      addGnssResidualMarginBlocks(state);

      // erase ambiguity state and GNSS measurements
      if (it_ambiguity != ambiguity_states_.end()) {
        ambiguity_states_.erase(it_ambiguity);
      }
      if (gnssMeasurementPairAt(state.timestamp) != gnss_measurement_pairs_.end()) {
        gnss_measurement_pairs_.erase(gnssMeasurementPairAt(state.timestamp));
      }

      // Erase state
      it = states_.erase(it);

      // just handle one in one time
      bool ret = applyMarginalization();

      return ret;
    }
  }

  return true;
}

// Sparsify GNSS states to bound computational load
void RtkImuCameraRrrEstimator::sparsifyGnssStates()
{
  // Check if we need to sparsify
  std::vector<BackendId> gnss_ids;
  std::vector<int> num_neighbors;
  int max_num_neighbors = 0;
  for (int i = 0; i < states_.size(); i++) {
    if (states_[i].id.type() != IdType::gPose) continue;
    gnss_ids.push_back(states_[i].id);
    // find neighbors
    num_neighbors.push_back(0);
    int idx = num_neighbors.size() - 1;
    if (i - 1 >= 0) for (int j = i - 1; j >= 0; j--) {
      if (states_[j].id.type() != IdType::gPose) break;
      num_neighbors[idx]++;
    } 
    if (i + 1 < states_.size()) for (int j = i; j < states_.size(); j++) {
      if (states_[j].id.type() != IdType::gPose) break;
      num_neighbors[idx]++;
    } 
    if (max_num_neighbors < num_neighbors[idx]) {
      max_num_neighbors = num_neighbors[idx];
    }
  }
  if (gnss_ids.size() <= rrr_options_.max_keyframes) return;

  // Erase some GNSS states
  // The states with most neighbors will be erased first
  int num_to_erase = gnss_ids.size() - rrr_options_.max_keyframes;
  std::vector<BackendId> ids_to_erase;
  for (int m = max_num_neighbors; m >= 0; m--) {
    for (size_t i = 0; i < num_neighbors.size(); i++) {
      if (num_neighbors[i] != m) continue;
      if (i == 0) continue;  // in case the first one connects to margin block
      ids_to_erase.push_back(gnss_ids[i]);
      if (ids_to_erase.size() >= num_to_erase) break;
    }
    if (ids_to_erase.size() >= num_to_erase) break;
  }
  CHECK(ids_to_erase.size() >= num_to_erase);
  for (int i = 0; i < states_.size(); i++) {
    if (std::find(ids_to_erase.begin(), ids_to_erase.end(), states_[i].id) 
        == ids_to_erase.end()) {
      continue;
    }
    State& state = states_[i];

    auto it_ambiguity = ambiguityStateAt(state.timestamp);
    // the first one may be connected with margin error
    if (it_ambiguity == ambiguity_states_.begin()) continue;

    eraseGnssMeasurementResidualBlocks(state);
    eraseFrequencyParameterBlocks(state);

    // we may failed to find corresponding ambiguity state because the GNSS state can be 
    // loosely coupled state during initialization.
    if (it_ambiguity != ambiguity_states_.end()) {
      eraseAmbiguityParameterBlocks(*it_ambiguity);
      ambiguity_states_.erase(it_ambiguity);
    }
    else {
      eraseGnssLooseResidualBlocks(state);
    }

    eraseImuState(state);

    if (gnssMeasurementPairAt(state.timestamp) != gnss_measurement_pairs_.end()) {
      gnss_measurement_pairs_.erase(gnssMeasurementPairAt(state.timestamp));
    }

    i--;
  }
}

// Compute ambiguity covariance at current epoch
bool RtkImuCameraRrrEstimator::estimateAmbiguityCovariance(
  const State& state, Eigen::MatrixXd& covariance)
{
  // Computing the ambiguity covariance directly is very time-consuming, so we run a parallel
  // optimizer and forcely set some empirical constraints to estimate a coarse covariance.
  Transformation T_WS = getPoseEstimate(state);
  Eigen::Vector3d t_SR_S = getGnssExtrinsicsEstimate();
  Eigen::Vector3d t_WR_W = T_WS.getPosition() + T_WS.getRotationMatrix() * t_SR_S;
  Eigen::Vector3d position = 
    coordinate_->convert(t_WR_W, GeoType::ENU, GeoType::ECEF);
  std::shared_ptr<Graph> sub_graph = ambiguity_covariance_estimator_->getGraph();
  const State sub_state = ambiguity_covariance_estimator_->getState();

  // Compute covariance
  std::vector<uint64_t> parameter_block_ids;
  // the ids in curAmbiguityState() are as the same as in ambiguity_covariance_estimator_
  for (auto id : curAmbiguityState().ids) {
    if (!sub_graph->parameterBlockExists(id.asInteger())) return false;
    parameter_block_ids.push_back(id.asInteger());
  }
  sub_graph->computeCovariance(parameter_block_ids, covariance);

  return true;
}



// ===== 核心函数1: buildDPOConstraints =====
bool RtkImuCameraRrrEstimator::buildDPOConstraints(
    const State& current_state, const FrameBundlePtr& current_frame) {
    
  if (!current_frame->isKeyframe()) return false;

  // 调试特征轨迹状态
  int total_tracks = dpo_manager_.active_tracks.size();
  int sufficient_tracks = 0;
  for (const auto& [fid, track] : dpo_manager_.active_tracks) {
    if (track.pose_sequence.size() >= 3) sufficient_tracks++;
  }
  
  LOG(INFO) << "DPO tracks status: total=" << total_tracks 
            << ", sufficient=" << sufficient_tracks;
  
  int constraints_added = 0;
  
  // 遍历有足够观测的特征轨迹
  for (auto& [feature_id, track] : dpo_manager_.active_tracks) {
    if (track.pose_sequence.size() < 3) continue;
    
    // 为该特征点的所有三元组创建约束
    auto triplets = track.generateTriplets();
    for (const auto& triplet : triplets) {
      if (addDPOConstraint(triplet, track, feature_id)) {
        constraints_added++;
      }
    }
  }
  
  LOG(INFO) << "DPO constraints built: " << constraints_added;
  return constraints_added > 0;
}


bool RtkImuCameraRrrEstimator::addDPOConstraint(
  const std::tuple<BackendId, BackendId, BackendId>& triplet,
  const DPOFeatureTrack& track, 
  int feature_id) {

  auto [pose_i_id, pose_j_id, pose_k_id] = triplet;
  
  // 检查是否已存在约束
  if (dpo_manager_.hasConstraint(triplet)) return false;
  
  // 找到对应的State对象
  State* state_i = nullptr, *state_j = nullptr, *state_k = nullptr;
  for (auto& s : states_) {
      if (s.id == pose_i_id) state_i = &s;
      if (s.id == pose_j_id) state_j = &s;
      if (s.id == pose_k_id) state_k = &s;
  }
  
  if (!state_i || !state_j || !state_k) return false;
  
  // 验证参数块存在
  if (!graph_->parameterBlockExists(state_i->id_in_graph.asInteger()) ||
      !graph_->parameterBlockExists(state_j->id_in_graph.asInteger()) ||
      !graph_->parameterBlockExists(state_k->id_in_graph.asInteger())) {
      return false;
  }
  
  // 获取观测数据
  auto it_i = track.observations.find(pose_i_id);
  auto it_j = track.observations.find(pose_j_id);
  auto it_k = track.observations.find(pose_k_id);
  
  if (it_i == track.observations.end() || 
      it_j == track.observations.end() || 
      it_k == track.observations.end()) {
      return false;
  }
  
  // 转换为归一化坐标
  Eigen::Vector3d coord_i(it_i->second[0], it_i->second[1], 1.0);
  Eigen::Vector3d coord_j(it_j->second[0], it_j->second[1], 1.0);
  Eigen::Vector3d coord_k(it_k->second[0], it_k->second[1], 1.0);
    
  // 准备参数块
  std::vector<std::shared_ptr<ParameterBlock>> parameter_blocks;
  parameter_blocks.push_back(graph_->parameterBlockPtr(state_i->id_in_graph.asInteger()));
  parameter_blocks.push_back(graph_->parameterBlockPtr(state_j->id_in_graph.asInteger()));
  parameter_blocks.push_back(graph_->parameterBlockPtr(state_k->id_in_graph.asInteger()));
  
  
  // 权重
  double weight = sqrt(rrr_options_.dpo_constraint_weight);
  Eigen::Matrix3d sqrt_information = weight * Eigen::Matrix3d::Identity();

  // 调试LOG
  LOG(INFO) << "DPO weight setting: " << rrr_options_.dpo_constraint_weight;
  LOG(INFO) << "Actual sqrt_info diagonal: " << sqrt_information(0,0);
  LOG(INFO) << "Effective information weight: " << sqrt_information(0,0) * sqrt_information(0,0);

  // 创建DPO约束
  auto dpo_error = std::make_shared<DPOConstraintError>(coord_i, coord_j, coord_k, sqrt_information);
  
  ceres::ResidualBlockId residual_id = graph_->addResidualBlock(
      dpo_error,
      cauchy_loss_function_.get(),
      parameter_blocks);

  if (residual_id != nullptr) {
    LOG(INFO) << "Successfully added DPO constraint";
  } else {
    LOG(WARNING) << "Failed to add DPO constraint!";
  }

  if (residual_id == nullptr) return false;
  
  // 创建约束记录
  DPOTripletConstraint constraint(pose_i_id, pose_j_id, pose_k_id, feature_id);
  constraint.residual_id = BackendId(reinterpret_cast<uint64_t>(residual_id));
  
  // 添加到管理器
  dpo_manager_.addConstraint(constraint);
  
  return true;
}

// ===== 核心函数2: updateFeatureTracksWithNewFrame =====
void RtkImuCameraRrrEstimator::updateFeatureTracksWithNewFrame(
  const BackendId& pose_id, const FrameBundlePtr& frame) {
  
  if (!frame || frame->size() == 0) return;
  
  const FramePtr& first_frame = frame->at(0);
  if (!first_frame) return;
  
  // ===== 添加调试代码 - 位置1 =====
  LOG(INFO) << "=== Feature Track Analysis ===";
  LOG(INFO) << "Active tracks before update: " << dpo_manager_.active_tracks.size();
  LOG(INFO) << "Frame features count: " << first_frame->numFeatures();
  int new_tracks = 0, existing_tracks = 0, invalid_features = 0;

  // 遍历所有特征点
  for (size_t i = 0; i < first_frame->numFeatures(); ++i) {
      int track_id = first_frame->track_id_vec_[i];
      
      
      // 跳过无效特征
      if (track_id < 0 || first_frame->type_vec_[i] == FeatureType::kOutlier) {
          continue;
      }
      
      // ===== 统计新轨迹vs现有轨迹 =====
      if (dpo_manager_.active_tracks.find(track_id) == dpo_manager_.active_tracks.end()) {
          new_tracks++;
      } else {
          existing_tracks++;
      }

      // 获取归一化坐标
      Eigen::Vector2d px = first_frame->px_vec_.col(i);
      Eigen::Vector3d f_normalized;
      if (!first_frame->cam()->backProject3(px, &f_normalized)) {
          continue;
      }
      
      Eigen::Vector2d norm_coord(f_normalized[0]/f_normalized[2], 
                                f_normalized[1]/f_normalized[2]);
      
      // 更新或创建轨迹
      auto& track = dpo_manager_.active_tracks[track_id];
      if (track.feature_id == 0) {  // 新轨迹
          track.feature_id = track_id;
          track.is_active = true;
      }
      
      track.observations[pose_id] = norm_coord;
      if (std::find(track.pose_sequence.begin(), track.pose_sequence.end(), 
                   pose_id) == track.pose_sequence.end()) {
          track.pose_sequence.push_back(pose_id);
      }
  }
      // ===== 添加调试代码 - 位置2 =====
      LOG(INFO) << "Feature breakdown: Invalid=" << invalid_features 
      << ", New tracks=" << new_tracks 
      << ", Existing tracks=" << existing_tracks;
      
      // 统计有足够观测的轨迹
      int sufficient_tracks = 0;
      for (auto& [fid, track] : dpo_manager_.active_tracks) {
          if (track.pose_sequence.size() >= 3) {
              sufficient_tracks++;
          }
      }
      LOG(INFO) << "Tracks with >=3 observations: " << sufficient_tracks;
      LOG(INFO) << "Active tracks after update: " << dpo_manager_.active_tracks.size();
      LOG(INFO) << "===============================";
}


// ===== 核心函数3:eraseDPOConstraints =====
void RtkImuCameraRrrEstimator::eraseDPOConstraints(const State& state) {
    
  // 删除涉及该pose的所有DPO约束
  // 遍历所有活跃的DPO约束，需要在循环中删除元素
  for (auto it = dpo_manager_.active_constraints.begin(); 
       it != dpo_manager_.active_constraints.end();) {
      // 解构获取pose_i, pose_j, pose_k
      auto [pose_i, pose_j, pose_k] = it->first;
      // 如果三元组中任意一个pose是要删除的pose
      if (pose_i == state.id || pose_j == state.id || pose_k == state.id) {
          // 获取residual ID，然后从优化图中移除这个残差块
          ceres::ResidualBlockId residual_id = 
              ceres::ResidualBlockId(it->second.residual_id.asInteger());
          graph_->removeResidualBlock(residual_id);
          // erase返回下一个迭代器
          it = dpo_manager_.active_constraints.erase(it);
      } else {
          ++it;
      }
  }
  
  // 清理特征轨迹中该pose的观测数据
  for (auto it = dpo_manager_.active_tracks.begin(); 
       it != dpo_manager_.active_tracks.end();) {
      
      auto& track = it->second;
      
      // 从observations map中移除该pose对应的2D观测点
      track.observations.erase(state.id);
      track.pose_sequence.erase(
          std::remove(track.pose_sequence.begin(), track.pose_sequence.end(), state.id),
          track.pose_sequence.end());
      
      // 如果轨迹观测数不足，删除整个轨迹
      if (track.pose_sequence.size() == 0) {
          it = dpo_manager_.active_tracks.erase(it);
      } else {
          ++it;
      }
  }
  
  // 清理pose-to-constraints映射，清理该pose相关的所有约束ID记录
  dpo_manager_.pose_to_constraints.erase(state.id);
}


// ===== 核心函数4: numDPOConstraints =====
size_t RtkImuCameraRrrEstimator::numDPOConstraints(const State& state) {
    
  size_t num = 0;
  
  // 统计涉及该pose的所有三元组约束
  for (const auto& [triplet, constraint] : dpo_manager_.active_constraints) {
      auto [pose_i, pose_j, pose_k] = triplet;
      
      if ((pose_i == state.id || pose_j == state.id || pose_k == state.id) 
          && !constraint.is_outlier) {
          num++;
      }
  }
  
  return num;
}


// ===== 核心函数5: rejectDPOConstraintOutlier =====
void RtkImuCameraRrrEstimator::rejectDPOConstraintOutlier(const State& state) {
    
  auto it = dpo_manager_.pose_to_constraints.find(state.id);
  if (it == dpo_manager_.pose_to_constraints.end()) return;
  
  std::vector<ceres::ResidualBlockId> outliers_to_remove;
  
  // 评估该pose相关的所有约束
  for (const auto& constraint_id : it->second) {
      ceres::ResidualBlockId residual_id = 
          ceres::ResidualBlockId(constraint_id.asInteger());
      
      // 计算残差
      Eigen::Vector3d residuals;
      graph_->problem()->EvaluateResidualBlock(
          residual_id, false, nullptr, residuals.data(), nullptr);
      
        double residual_norm = residuals.norm();
        LOG(INFO) << "DPO final residual norm: " << residual_norm;  // 添加这行

      // 检查是否为outlier
      if (residuals.norm() > rrr_options_.dpo_outlier_threshold) {
          outliers_to_remove.push_back(residual_id);
      }
  }
  
  // 删除outlier约束并保持数据一致性
  for (auto residual_id : outliers_to_remove) {
    graph_->removeResidualBlock(residual_id);
    
    // 找到约束记录并获取涉及的所有pose
    BackendId residual_backend_id(reinterpret_cast<uint64_t>(residual_id));
    
    for (auto constraint_it = dpo_manager_.active_constraints.begin();
          constraint_it != dpo_manager_.active_constraints.end(); ++constraint_it) {
        
        if (constraint_it->second.residual_id.asInteger() == 
            reinterpret_cast<uint64_t>(residual_id)) {
            
            // 获取三元组中的所有pose
            auto [pose_i, pose_j, pose_k] = constraint_it->first;
            
            // 从所有三个pose的映射中删除该约束
            auto remove_constraint = [&](const BackendId& pose_id) {
                // 获取该pose对应的约束列表
                auto& constraints = dpo_manager_.pose_to_constraints[pose_id];
                constraints.erase(
                    std::remove(constraints.begin(), constraints.end(), residual_backend_id), 
                    constraints.end());
            };
            
            remove_constraint(pose_i);
            remove_constraint(pose_j);
            remove_constraint(pose_k);
            
            // 删除主约束记录
            dpo_manager_.active_constraints.erase(constraint_it);
            break;
        }
    }
  }
}

};