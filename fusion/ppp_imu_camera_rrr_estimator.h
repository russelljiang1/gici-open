/**
* @Function: PPP/IMU/Camera tightly couple estimator (GNSS raw (PPP formula) + IMU raw + camera raw)
*
* @Author  : Cheng Chi
* @Email   : chichengcn@sjtu.edu.cn
*
* Copyright (C) 2023 by Cheng Chi, All rights reserved.
**/
#pragma once

#include "gici/gnss/gnss_estimator_base.h"
#include "gici/imu/imu_estimator_base.h"
#include "gici/vision/visual_estimator_base.h"
#include "gici/gnss/ppp_estimator.h"
#include "gici/fusion/gnss_imu_initializer.h"
#include "gici/fusion/rtk_imu_camera_rrr_estimator.h"

namespace gici {

// PPP/IMU/Camera RRR couple options
struct PppImuCameraRrrEstimatorOptions {
  // Frame state window length
  // We only keep GNSS measurements near to keyframes (one-to-one) and throw the others
  // away after one optimization, because the GNSS measurement errors, especially for
  // the multipath, are highly correlated between epochs when we have a slow or zero motion.
  // Besides, we need at least 2 GNSS states in window. If current setting cannot ensure
  // this condition, we will ignore this option and extend the windows length.
  int max_keyframes = 5;

  // GNSS state window length before visual has been initialized
  int max_gnss_window_length_minor = 3;

  // Maximum yaw STD to start visual initialization (deg)
  double min_yaw_std_init_visual = 0.5;

  // ===== 新增：滑动窗口DPO相关参数 =====
  // 是否使用DPO方法替代传统landmark方法
  bool use_dpo_constraints = true;
  
  // 滑动窗口配置
  int dpo_window_size = 5;                         // DPO滑动窗口大小
  int dpo_max_constraints_per_feature = 10;        // 每个特征点最大约束数
  int dpo_min_observations_per_track = 3;          // 特征轨迹最小观测数

  // DPO约束质量控制
  double dpo_constraint_weight;              // DPO权重
  double dpo_min_parallax_angle = 0.5;             // 最小视差角(degree)
  double dpo_huber_loss = 1.0;                     // DPO约束Huber loss参数
  int dpo_min_features = 8;                        // DPO方法要求的最小特征数
  double dpo_outlier_threshold = 2.0;              // DPO约束outlier检测阈值
  double dpo_max_reject_ratio = 0.3;               // DPO约束的最大reject比例
  
  // 约束选择策略
  enum DPOConstraintSelectionStrategy {
    BEST_PARALLAX,     // 选择视差角最大的约束
    UNIFORM_SAMPLING,  // 均匀采样约束
    ADAPTIVE_QUALITY   // 自适应质量选择
  };
  DPOConstraintSelectionStrategy dpo_selection_strategy = BEST_PARALLAX;
  
  // 性能优化
  bool dpo_enable_constraint_caching = true;       // 启用约束缓存
  bool dpo_enable_parallel_processing = true;      // 启用并行处理
  bool dpo_relaxed_initialization = true;          // 使用更宽松的初始化条件
  bool dpo_use_information_matrix = false; 
  
  // 调试和日志
  bool dpo_verbose_statistics = false;             // 输出详细统计信息
  int dpo_log_frequency = 10;                      // 统计信息输出频率

};

// Estimator
class PppImuCameraRrrEstimator :
  public GnssEstimatorBase,
  public VisualEstimatorBase,
  public ImuEstimatorBase {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  PppImuCameraRrrEstimator(const PppImuCameraRrrEstimatorOptions& options,
               const GnssImuInitializerOptions& init_options,
               const PppEstimatorOptions ppp_options,
               const GnssEstimatorBaseOptions& gnss_base_options,
               const GnssLooseEstimatorBaseOptions& gnss_loose_base_options,
               const VisualEstimatorBaseOptions& visual_base_options,
               const ImuEstimatorBaseOptions& imu_base_options,
               const EstimatorBaseOptions& base_options,
               const AmbiguityResolutionOptions& ambiguity_options);
  ~PppImuCameraRrrEstimator();

  // Add measurement
  bool addMeasurement(const EstimatorDataCluster& measurement) override;

  // Estimate current graph
  bool estimate() override;

  // Set initializatin result
  void setInitializationResult(
    const std::shared_ptr<MultisensorInitializerBase>& initializer) override;

protected:
  // Add GNSS measurements and state
  bool addGnssMeasurementAndState(const GnssMeasurement& measurement);

  // Add image measurements and state
  bool addImageMeasurementAndState(const FrameBundlePtr& frame_bundle,
    const SpeedAndBias& speed_and_bias = SpeedAndBias::Zero());

  // Visual initialization
  bool visualInitialization(const FrameBundlePtr& frame_bundle);

  // Marginalization
  bool marginalization(const IdType& type);

  // Marginalization when the new state is a frame state
  bool frameMarginalization();

  // Marginalization when the new state is a GNSS state
  bool gnssMarginalization();

  // Sparsify GNSS states to bound computational load
  void sparsifyGnssStates();

  // Compute ambiguity covariance at current epoch
  bool estimateAmbiguityCovariance(const State& state, Eigen::MatrixXd& covariance);

  // Get latest state
  inline State& latestState() override { return states_[latest_state_index_]; }

  // 构建DPO约束
  bool buildDPOConstraints(const State& current_state, const FrameBundlePtr& current_frame);

  bool addDPOConstraint(const std::tuple<BackendId, BackendId, BackendId>& triplet,
    const DPOFeatureTrack& track, 
    int feature_id);
  
  // 更新特征点轨迹
  void updateFeatureTracksWithNewFrame(const BackendId& pose_id, const FrameBundlePtr& frame);
  
  // DPO约束清理
  void eraseDPOConstraints(const State& state);
  
  // 统计与指定frame相关的DPO约束数量
  size_t numDPOConstraints(const State& state);
  
  // DPO约束的outlier rejection
  void rejectDPOConstraintOutlier(const State& state);

protected:
  // Options
  PppImuCameraRrrEstimatorOptions rrr_options_;
  PppEstimatorOptions ppp_options_;

  // SPP estimator to get initial states
  std::unique_ptr<SppEstimator> spp_estimator_;

  // Phase wind-up handle
  PhaseWindupPtr phase_windup_;

  // Initialization control
  std::shared_ptr<GnssImuInitializer> gnss_imu_initializer_;
  std::shared_ptr<PppEstimator> initializer_sub_estimator_;
  bool visual_initialized_ = false;
  std::deque<FrameBundlePtr> init_keyframes_;
  std::deque<Solution> init_solution_store_;

  // PPP estimator used for ambiguity covariance estimation
  std::unique_ptr<PppEstimator> ambiguity_covariance_estimator_;
  bool ambiguity_covariance_coordinate_setted_ = false;

  // Status control
  int num_continuous_unfix_ = 0;
  int num_cotinuous_reject_gnss_ = 0;
  int num_cotinuous_reject_visual_ = 0;


  // DPO约束管理器
  DPOConstraintManager dpo_manager_;
  
};

}
