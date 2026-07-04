/**
* @Function: RTK/IMU/Camera tightly couple estimator (GNSS raw (RTK formula) + IMU raw + camera raw)
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
#include "gici/gnss/rtk_estimator.h"
#include "gici/fusion/gnss_imu_initializer.h"
#include <ceres/ceres.h>
#include <opencv2/opencv.hpp>
#include <map>
#include <deque>
#include <tuple>

namespace gici {

// 特征点跟踪轨迹
struct DPOFeatureTrack {
  int feature_id;                                      // 特征点唯一ID
  std::map<BackendId, Eigen::Vector2d> observations;   // pose_id -> 归一化坐标观测
  std::vector<BackendId> pose_sequence;                // 观测到该特征点的pose序列
  bool is_active = true;                               // 是否活跃跟踪中
  
  // 获取该特征点可以形成的所有三元组
  std::vector<std::tuple<BackendId, BackendId, BackendId>> generateTriplets() const {
    std::vector<std::tuple<BackendId, BackendId, BackendId>> triplets;
    if (pose_sequence.size() < 3) return triplets;
    
    // 生成所有可能的三元组合
    for (size_t i = 0; i < pose_sequence.size() - 2; i++) {
      for (size_t j = i + 1; j < pose_sequence.size() - 1; j++) {
        for (size_t k = j + 1; k < pose_sequence.size(); k++) {
          triplets.emplace_back(pose_sequence[i], pose_sequence[j], pose_sequence[k]);
        }
      }
    }
    return triplets;
  }
};

// 三元组约束信息
struct DPOTripletConstraint {
  BackendId pose_i, pose_j, pose_k;                 // 三个pose
  int feature_id;                                   // 对应的特征点ID
  BackendId residual_id;                            // 在图中的残差ID
  bool is_outlier = false;                          // 是否为outlier
  
  // 默认构造函数
  DPOTripletConstraint() : feature_id(-1) {}

  // 构造函数
  DPOTripletConstraint(const BackendId& pi, const BackendId& pj, const BackendId& pk, 
                       int fid) 
    : pose_i(pi), pose_j(pj), pose_k(pk), feature_id(fid) {}
  
  // 获取三元组标识
  std::tuple<BackendId, BackendId, BackendId> getTriplet() const {
    return std::make_tuple(pose_i, pose_j, pose_k);
  }
};

// DPO约束管理器
struct DPOConstraintManager {
  // 活跃的特征点轨迹
  std::map<int, DPOFeatureTrack> active_tracks;
  
  // 三元组约束记录 (triplet -> constraint_info)
  std::map<std::tuple<BackendId, BackendId, BackendId>, DPOTripletConstraint> active_constraints;
  
  // pose到约束的映射关系 (pose_id -> constraint_residual_ids)
  std::map<BackendId, std::vector<BackendId>> pose_to_constraints;

  // 检查三元组是否已存在约束
  bool hasConstraint(const std::tuple<BackendId, BackendId, BackendId>& triplet) const {
    return active_constraints.find(triplet) != active_constraints.end();
  }
  
  // 添加约束记录
  void addConstraint(const DPOTripletConstraint& constraint) {
    active_constraints[constraint.getTriplet()] = constraint;
    
    // 更新pose->constraints映射
    pose_to_constraints[constraint.pose_i].push_back(constraint.residual_id);
    pose_to_constraints[constraint.pose_j].push_back(constraint.residual_id);
    pose_to_constraints[constraint.pose_k].push_back(constraint.residual_id);
  }
  
  // 清空管理器
  void clear() {
    active_tracks.clear();
    active_constraints.clear();
    pose_to_constraints.clear();
  }
 
};


// ===== DPO约束误差类 =====
class DPOConstraintError :
    public ceres::SizedCostFunction<3, 7, 7, 7>,
    public ErrorInterface {
    
private:
    Eigen::Vector3d X_i_, X_j_, X_k_;
    Eigen::Matrix3d sqrt_information_;
    Eigen::Matrix3d sqrt_information_inverse_;

    struct DPOFunctor {
        Eigen::Vector3d X_i, X_j, X_k;

        DPOFunctor(const Eigen::Vector3d& X_i, const Eigen::Vector3d& X_j, const Eigen::Vector3d& X_k)
        : X_i(X_i), X_j(X_j), X_k(X_k) {}
        
        template <typename T>
        bool operator()(const T* const pose_i,
                        const T* const pose_j,
                        const T* const pose_k,
                        T* residuals) const {
            
            // 提取平移部分
            Eigen::Matrix<T, 3, 1> t_i(pose_i[0], pose_i[1], pose_i[2]);
            Eigen::Matrix<T, 3, 1> t_j(pose_j[0], pose_j[1], pose_j[2]);
            Eigen::Matrix<T, 3, 1> t_k(pose_k[0], pose_k[1], pose_k[2]);
            
            // 从四元数提取旋转矩阵
            Eigen::Quaternion<T> q_i(pose_i[6], pose_i[3], pose_i[4], pose_i[5]);
            Eigen::Quaternion<T> q_j(pose_j[6], pose_j[3], pose_j[4], pose_j[5]);
            Eigen::Quaternion<T> q_k(pose_k[6], pose_k[3], pose_k[4], pose_k[5]);
            
            Eigen::Matrix<T, 3, 3> R_i = q_i.toRotationMatrix();
            Eigen::Matrix<T, 3, 3> R_j = q_j.toRotationMatrix();
            Eigen::Matrix<T, 3, 3> R_k = q_k.toRotationMatrix();
            
            // 计算相对变换
            Eigen::Matrix<T, 3, 3> R_ij = R_j * R_i.transpose();
            Eigen::Matrix<T, 3, 3> R_ik = R_k * R_i.transpose();
            Eigen::Matrix<T, 3, 1> t_ij = t_j - t_i;
            Eigen::Matrix<T, 3, 1> t_ik = t_k - t_i;

            /*
            if (std::is_same<T, double>::value) {
              // 输出平移向量
              std::cout << "=== DPO Debug Info ===" << std::endl;
              std::cout << "t_ij: [" << t_ij[0] << ", " << t_ij[1] << ", " << t_ij[2] 
                        << "] norm=" << t_ij.norm() << std::endl;
              std::cout << "t_ik: [" << t_ik[0] << ", " << t_ik[1] << ", " << t_ik[2] 
                        << "] norm=" << t_ik.norm() << std::endl;
              
              // 将旋转矩阵转换为轴角表示
              Eigen::AngleAxis<T> aa_ij(R_ij);
              Eigen::AngleAxis<T> aa_ik(R_ik);
              
              T angle_ij = aa_ij.angle() * 180.0 / M_PI;  // 转为度数
              T angle_ik = aa_ik.angle() * 180.0 / M_PI;
              
              std::cout << "R_ij angle: " << angle_ij << "° axis: [" 
                        << aa_ij.axis()[0] << ", " << aa_ij.axis()[1] << ", " << aa_ij.axis()[2] << "]" << std::endl;
              std::cout << "R_ik angle: " << angle_ik << "° axis: [" 
                        << aa_ik.axis()[0] << ", " << aa_ik.axis()[1] << ", " << aa_ik.axis()[2] << "]" << std::endl;
            }
            */

            // 动态计算theta, a, b系数
            Eigen::Matrix<T, 3, 1> X_i_T = X_i.template cast<T>();
            Eigen::Matrix<T, 3, 1> X_j_T = X_j.template cast<T>();
            Eigen::Matrix<T, 3, 1> X_k_T = X_k.template cast<T>();
            
            // 计算theta_ij_sq, theta_ik_sq
            Eigen::Matrix<T, 3, 3> X_j_skew = skewMatrix(X_j_T);
            Eigen::Matrix<T, 3, 3> X_k_skew = skewMatrix(X_k_T);
            
            Eigen::Matrix<T, 3, 1> theta_ij = X_j_skew * (R_ij * X_i_T);
            Eigen::Matrix<T, 3, 1> theta_ik = X_k_skew * (R_ik * X_i_T);
            
            T theta_ij_sq = theta_ij.squaredNorm();
            T theta_ik_sq = theta_ik.squaredNorm();
            
            // 计算a_ik, b_ij系数
            Eigen::Matrix<T, 3, 1> a_ik = computeCoeffA(X_i_T, X_k_T, R_ik);
            Eigen::Matrix<T, 3, 1> b_ij = computeCoeffB(X_i_T, X_j_T, R_ij);
           
            /*
            // ===== 添加LOG验证 =====
            if (std::is_same<T, double>::value) {
              double theta_ij_val = *(reinterpret_cast<const double*>(&theta_ij_sq));
              double theta_ik_val = *(reinterpret_cast<const double*>(&theta_ik_sq));
              
              std::cout << "DPO Debug: theta_ij_sq=" << theta_ij_val 
                        << ", theta_ik_sq=" << theta_ik_val << std::endl;
              std::cout << "a_ik norm=" << a_ik.norm() << ", b_ij norm=" << b_ij.norm() << std::endl;
            }
            */

            // 第一项：θ²_{i,j} * R_{i,j} * X_i * (a^T_{i,k} * t_{i,k})
            T scalar_a_t_ik = a_ik.dot(t_ik);
            Eigen::Matrix<T, 3, 1> term1 = T(theta_ij_sq) * (R_ij * X_i_T) * scalar_a_t_ik;
            
            // 第二项：θ²_{i,k} * θ²_{i,j} * t_{i,j}  
            Eigen::Matrix<T, 3, 1> term2 = T(theta_ik_sq * theta_ij_sq) * t_ij; 
            
            // 第三项：θ²_{i,k} * X_j * (b^T_{i,j} * t_{i,j})
            T scalar_b_t_ij = b_ij.dot(t_ij);
            Eigen::Matrix<T, 3, 1> term3 = T(theta_ik_sq) * X_j_T * scalar_b_t_ij;

            /*
            // DPO约束：term1 + term2 - term3 = 0
            Eigen::Matrix<T, 3, 1> residual_vec = term1 + term2 - term3;

            if (std::is_same_v<T, double>) {
              double r0 = *(reinterpret_cast<const double*>(&residual_vec[0]));
              double r1 = *(reinterpret_cast<const double*>(&residual_vec[1]));
              double r2 = *(reinterpret_cast<const double*>(&residual_vec[2]));
              std::cout << "estimating residual: [" << r0 << ", " << r1 << ", " << r2 << "]" << std::endl;
            }
            
            residuals[0] = residual_vec[0];
            residuals[1] = residual_vec[1];
            residuals[2] = residual_vec[2];
            */

            const T epsilon = T(1e-10);  // 防除零
    
            residuals[0] = (term1[0] + term2[0]) / (abs(term3[0]) + epsilon) - T(1.0);
            residuals[1] = (term1[1] + term2[1]) / (abs(term3[1]) + epsilon) - T(1.0);
            residuals[2] = (term1[2] + term2[2]) / (abs(term3[2]) + epsilon) - T(1.0);

            if (std::is_same<T, double>::value) {
              double r0 = *(reinterpret_cast<const double*>(&residuals[0]));
              double r1 = *(reinterpret_cast<const double*>(&residuals[1]));
              double r2 = *(reinterpret_cast<const double*>(&residuals[2]));
              //std::cout << "Final residual: [" << r0 << ", " << r1 << ", " << r2 << "]" << std::endl;
              //std::cout << "Residual norm: " << sqrt(r0*r0 + r1*r1 + r2*r2) << std::endl;
            }
            

            return true;
        }

    private:
        template<typename T>
        Eigen::Matrix<T, 3, 3> skewMatrix(const Eigen::Matrix<T, 3, 1>& v) const {
            Eigen::Matrix<T, 3, 3> skew;
            skew << T(0), -v[2], v[1],
                    v[2], T(0), -v[0],
                    -v[1], v[0], T(0);
            return skew;
          }

        template<typename T>
        Eigen::Matrix<T, 3, 1> computeCoeffA(const Eigen::Matrix<T, 3, 1>& X_i, 
                                            const Eigen::Matrix<T, 3, 1>& X_k, 
                                            const Eigen::Matrix<T, 3, 3>& R_ik) const {
            Eigen::Matrix<T, 3, 1> R_ik_Xi = R_ik * X_i;
            Eigen::Matrix<T, 3, 3> R_ik_Xi_skew = skewMatrix(R_ik_Xi);
            Eigen::Matrix<T, 3, 1> cross_product = R_ik_Xi_skew * X_k;
            Eigen::Matrix<T, 3, 3> X_k_skew = skewMatrix(X_k);
            return X_k_skew * cross_product;
          }

        template<typename T>
        Eigen::Matrix<T, 3, 1> computeCoeffB(const Eigen::Matrix<T, 3, 1>& X_i, 
                                            const Eigen::Matrix<T, 3, 1>& X_j, 
                                            const Eigen::Matrix<T, 3, 3>& R_ij) const {
            Eigen::Matrix<T, 3, 1> R_ij_Xi = R_ij * X_i;
            Eigen::Matrix<T, 3, 3> R_ij_Xi_skew = skewMatrix(R_ij_Xi);
            Eigen::Matrix<T, 3, 1> cross_product = R_ij_Xi_skew * X_j;
            return R_ij_Xi_skew.transpose() * cross_product;
          }
    };
    
    std::unique_ptr<ceres::CostFunction> auto_diff_cost_function_;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
    typedef ceres::SizedCostFunction<3, 7, 7, 7> base_t;
    static const int kNumResiduals = 3;
    typedef Eigen::Matrix<double, 3, 3> information_t;

    DPOConstraintError(const Eigen::Vector3d& X_i,
      const Eigen::Vector3d& X_j,
      const Eigen::Vector3d& X_k,
      const Eigen::Matrix3d& sqrt_information = Eigen::Matrix3d::Identity()) 
        : X_i_(X_i), X_j_(X_j), X_k_(X_k), sqrt_information_(sqrt_information) {
        
        sqrt_information_inverse_ = sqrt_information_.inverse();
        
        // 自动微分
        auto_diff_cost_function_.reset(
            new ceres::AutoDiffCostFunction<DPOFunctor, 3, 7, 7, 7>(
                new DPOFunctor( X_i_, X_j_, X_k_)
            )
        );
    }
    
    // ========== ceres::CostFunction 接口 ==========
    virtual bool Evaluate(double const* const* parameters,
                          double* residuals,
                          double** jacobians) const override {
        
        if (!auto_diff_cost_function_->Evaluate(parameters, residuals, jacobians)) {
            return false;
        }
        
        Eigen::Map<Eigen::Vector3d> residual_map(residuals);
        residual_map = sqrt_information_ * residual_map;
        
        if (jacobians != nullptr) {
            for (int i = 0; i < 3; i++) { 
                if (jacobians[i] != nullptr) {
                    Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>> jac_map(jacobians[i]);
                    jac_map = sqrt_information_ * jac_map;
                }
            }
        }
        
        return true;
    }
    
    // ========== ErrorInterface 接口 ==========
    virtual bool EvaluateWithMinimalJacobians(
        double const* const* parameters,
        double* residuals,
        double** jacobians,
        double** jacobians_minimal) const override {
        
        return Evaluate(parameters, residuals, jacobians);
    }
    
    virtual size_t residualDim() const override { 
        return kNumResiduals; 
    }
    
    virtual size_t parameterBlocks() const override {
        return base_t::parameter_block_sizes().size();
    }
    
    virtual size_t parameterBlockDim(size_t parameter_block_idx) const override {
        return base_t::parameter_block_sizes().at(parameter_block_idx);
    }
    
    virtual ErrorType typeInfo() const override {
        return ErrorType::kDPOConstraintError;
    }
    
    virtual void deNormalizeResidual(double *residuals) const override {
        Eigen::Map<Eigen::Vector3d> residual_map(residuals);
        residual_map = sqrt_information_inverse_ * residual_map;
    }
    
    // 静态创建函数
    static std::shared_ptr<DPOConstraintError> Create(
        const Eigen::Vector3d& X_i,
        const Eigen::Vector3d& X_j,
        const Eigen::Vector3d& X_k,
        const Eigen::Matrix3d& sqrt_information = Eigen::Matrix3d::Identity()) {
        return std::make_shared<DPOConstraintError>(X_i, X_j, X_k, sqrt_information);
    }

};

  
// RTK/IMU/Camera RRR couple options
struct RtkImuCameraRrrEstimatorOptions {
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
class RtkImuCameraRrrEstimator : 
  public GnssEstimatorBase, 
  public VisualEstimatorBase, 
  public ImuEstimatorBase {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  RtkImuCameraRrrEstimator(const RtkImuCameraRrrEstimatorOptions& options, 
               const GnssImuInitializerOptions& init_options, 
               const RtkEstimatorOptions rtk_options,
               const GnssEstimatorBaseOptions& gnss_base_options, 
               const GnssLooseEstimatorBaseOptions& gnss_loose_base_options, 
               const VisualEstimatorBaseOptions& visual_base_options,
               const ImuEstimatorBaseOptions& imu_base_options,
               const EstimatorBaseOptions& base_options,
               const AmbiguityResolutionOptions& ambiguity_options);
  ~RtkImuCameraRrrEstimator();

  // Add measurement
  bool addMeasurement(const EstimatorDataCluster& measurement) override;

  // Estimate current graph
  bool estimate() override;

  // Set initializatin result
  void setInitializationResult(
    const std::shared_ptr<MultisensorInitializerBase>& initializer) override;

protected:

  // Add GNSS measurements and state
  bool addGnssMeasurementAndState(
    const GnssMeasurement& measurement_rov, 
    const GnssMeasurement& measurement_ref);

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
  RtkImuCameraRrrEstimatorOptions rrr_options_;
  RtkEstimatorOptions rtk_options_;

  // Initialization control
  std::shared_ptr<GnssImuInitializer> gnss_imu_initializer_;
  std::shared_ptr<RtkEstimator> initializer_sub_estimator_;
  bool visual_initialized_ = false;
  std::deque<FrameBundlePtr> init_keyframes_;
  std::deque<Solution> init_solution_store_;

  // Measurement alignment handle
  DifferentialMeasurementsAlign meausrement_align_;

  // RTK estimator used for ambiguity covariance estimation
  std::unique_ptr<RtkEstimator> ambiguity_covariance_estimator_;
  bool ambiguity_covariance_coordinate_setted_ = false;

  // Status control
  int num_continuous_unfix_ = 0;
  int num_cotinuous_reject_gnss_ = 0;
  int num_cotinuous_reject_visual_ = 0;


  // DPO约束管理器
  DPOConstraintManager dpo_manager_;
  
};

}