/**
* @Function: Stream functions
*
* @Author  : Cheng Chi
* @Email   : chichengcn@sjtu.edu.cn
*
* Copyright (C) 2023 by Cheng Chi, All rights reserved.
**/
#include "gici/stream/formator.h"

#include <sstream>
#include <mutex>
#include <glog/logging.h>
#include <vikit/timer.h>

#include "gici/gnss/gnss_common.h"
#include "gici/utility/transform.h"

namespace gici {

DataCluster::DataCluster(FormatorType type)
{
  if (type == FormatorType::RTCM2 || type == FormatorType::RTCM3 ||
    type == FormatorType::GnssRaw || type == FormatorType::RINEX ||
    type == FormatorType::DcbFile || type == FormatorType::AtxFile || 
    type == FormatorType::Sp3File || type == FormatorType::ClkFile) {
    gnss = std::make_shared<GNSS>();
    gnss->init();
    return;
  }
  if (type == FormatorType::IMUPack || type == FormatorType::IMUText) {
    imu = std::make_shared<IMU>();
    return;
  }
  if (type == FormatorType::OptionPack) {
    option = std::make_shared<Option>();
    return;
  }
  if (type == FormatorType::NMEA) {
    solution = std::make_shared<Solution>();
    return;
  }
  if (type == FormatorType::ImagePack || type == FormatorType::ImageV4L2) {
    LOG(FATAL) << "Cannot initialize DataCluster::Image: "
           << "Image length should be given!";
  }
  LOG(FATAL) << "Cannot initialize: Data format not recognized!";
}

DataCluster::DataCluster(FormatorType type, int _width, int _height, int _step)
{
  if (type == FormatorType::ImagePack || type == FormatorType::ImageV4L2) {
    image = std::make_shared<Image>();
    image->init(_width, _height, _step);
    return;
  }
  LOG(FATAL) << "Cannot initialize: Data format not recognized!";
}

DataCluster::~DataCluster()
{
  if (gnss != nullptr) gnss->free();
  if (image != nullptr) image->free();
}

void DataCluster::GNSS::init()
{
  released = false;
  if (!(observation = (obs_t *)malloc(sizeof(obs_t))) ||
    !(observation->data = (obsd_t *)malloc(sizeof(obsd_t) * MAXOBS)) ||
    !(ephemeris = (nav_t *)malloc(sizeof(nav_t))) ||
    !(antenna = (sta_t *)malloc(sizeof(sta_t))) ) {
    free();
  }

  memset(ephemeris, 0, sizeof(nav_t));
  memset(antenna, 0, sizeof(sta_t));
}

void DataCluster::GNSS::free()
{
  if (released) return;
  released = true;

  if (observation) {
    if (observation->data) {
      ::free(observation->data);
      observation->data = NULL;
    }
    ::free(observation);
    observation = NULL;
  }

  if (ephemeris) {
    if (ephemeris->eph) {
      ::free(ephemeris->eph);
      ephemeris->eph = NULL;
    }
    if (ephemeris->geph) {
      ::free(ephemeris->geph);
      ephemeris->geph = NULL;
    }
    if (ephemeris->seph) {
      ::free(ephemeris->seph);
      ephemeris->seph = NULL;
    }
    if (ephemeris->peph) {
      ::free(ephemeris->peph);
      ephemeris->peph = NULL;
    }
    if (ephemeris->pclk) {
      ::free(ephemeris->pclk);
      ephemeris->pclk = NULL;
    }
    ::free(ephemeris);
    ephemeris = NULL;
  }

  if (antenna) {
    ::free(antenna);
    antenna = NULL;
  }

  if (receiver_pcvs) {
    if (receiver_pcvs->pcv) {
      ::free(receiver_pcvs->pcv);
      receiver_pcvs->pcv = NULL;
    }
    ::free(receiver_pcvs);
    receiver_pcvs = NULL;
  }
}

void DataCluster::Image::init(int _width, int _height, int _step)
{
  released = false;
  width = _width;
  height = _height;
  step = _step;
  if (!(image = (uint8_t *)malloc(sizeof(uint8_t) * width * height * step)))
    free();
}

void DataCluster::Image::free()
{
  if (released) return;
  released = true;
  if (image) {
    ::free(image);
    image = NULL;
  }
}

// Load option with info
#define LOAD_COMMON(opt) \
  if (!option_tools::safeGet(node, #opt, &option.opt)) { \
  LOG(INFO) << __FUNCTION__ << ": Unable to load " << #opt \
         << ". Using default instead."; }
// Load option with fatal error
#define LOAD_REQUIRED(opt) \
  if (!option_tools::safeGet(node, #opt, &option.opt)) { \
  LOG(FATAL) << __FUNCTION__ << ": Unable to load " << #opt << "!"; }

// RTCM 2 ----------------------------------------------------
// Load date
inline bool loadStartTime(const YAML::Node& node, double& start_time) {
  std::string str;
  if (!option_tools::safeGet(node, "start_time", &str)) {
    LOG(INFO) << __FUNCTION__ << ": Unable to load " << "start_time"
         << ". Using default instead.";
    return false;
  }
  std::string strs[4];
  int index = 0;
  for (size_t i = 0; i < str.size(); i++) {
    if (str[i] == '.') {
      index++; continue;
    }
    strs[index] = strs[index] + str[i];
  }
  CHECK(index == 2) << "Start time format illegal!";
  int year = atoi(strs[0].data());
  int month = atoi(strs[1].data());
  int day = atoi(strs[2].data());
  int hour = 12;
  int min = 0, sec = 0;
  double ep[] = { (double)year, (double)month, (double)day,
                  (double)hour, (double)min, (double)sec };
  start_time = gnss_common::gtimeToDouble(epoch2time(ep));

  return true;
}

RTCM2Formator::RTCM2Formator(const Option& option)
{
  type_ = FormatorType::RTCM2;

  memset(&rtcm_, 0, sizeof(rtcm_t));
  init_rtcm(&rtcm_);
  sprintf(rtcm_.opt, "-EPHALL");
  rtcm_.time = gnss_common::doubleToGtime(option.start_time);
  for (int i = 0; i < MaxDataSize::RTCM2; i++) {
    data_.push_back(std::make_shared<DataCluster>(type_));
  }
}

RTCM2Formator::RTCM2Formator(const YAML::Node& node)
{
  Option option;
  option.start_time = vk::Timer::getCurrentTime();
  loadStartTime(node, option.start_time);

  type_ = FormatorType::RTCM2;

  memset(&rtcm_, 0, sizeof(rtcm_t));
  init_rtcm(&rtcm_);
  sprintf(rtcm_.opt, "-EPHALL");
  rtcm_.time = gnss_common::doubleToGtime(option.start_time);
  for (int i = 0; i < MaxDataSize::RTCM2; i++) {
    data_.push_back(std::make_shared<DataCluster>(type_));
  }
}

RTCM2Formator::~RTCM2Formator()
{
  free_rtcm(&rtcm_);
}

// Decode stream to data
int RTCM2Formator::decode(const uint8_t *buf, int size, 
    std::vector<std::shared_ptr<DataCluster>>& data)
{
  // Clear old informations and get GNSS data handle
  std::vector<std::shared_ptr<DataCluster::GNSS>> gnss_data;
  for (size_t i = 0; i < data_.size(); i++) {
    data_[i]->gnss->types.clear();
    gnss_data.push_back(data_[i]->gnss);
  }

  bool is_observation = false;
  bool has_others = false;
  int iobs = 0;
  for (int i = 0; i < size; i++) {
    int ret = input_rtcm2(&rtcm_, buf[i]);
    if (ret <= 0) continue;

    obs_t *obs = &rtcm_.obs;
    nav_t *nav = &rtcm_.nav;
    sta_t *sta = &rtcm_.sta;
    ssr_t *ssr = rtcm_.ssr;
    int sat = rtcm_.ephsat;
    gnss_common::updateStreamData(
        ret, obs, nav, sta, ssr, iobs, sat, gnss_data);
    GnssDataType type = static_cast<GnssDataType>(ret);
    std::shared_ptr<DataCluster::GNSS>& gnss = 
      type == GnssDataType::Observation ? gnss_data[iobs] : gnss_data[0];
    if (std::find(gnss->types.begin(), gnss->types.end(), type)
      == gnss->types.end()) {
      gnss->types.push_back(type);
    }
    
    if (type == GnssDataType::Observation) {
      if (iobs < MaxDataSize::RTCM2) iobs++;
      if (iobs >= MaxDataSize::RTCM2) {
        LOG(WARNING) << "Max data length surpassed!";
        break;
      }
      is_observation = true;
    }
    else {
      has_others = true;
    }
  }

  data = data_;

  return is_observation ? iobs : has_others;
}

// Encode data to stream
int RTCM2Formator::encode(const std::shared_ptr<DataCluster>& data, uint8_t *buf)
{
  LOG(ERROR) << "RTCM2 Encoding not supported!";
  return 0;
}

// RTCM 3 -------------------------------------------------
RTCM3Formator::RTCM3Formator(const Option& option)
{
  type_ = FormatorType::RTCM3;

  memset(&rtcm_, 0, sizeof(rtcm_t));
  init_rtcm(&rtcm_);
  sprintf(rtcm_.opt, "-EPHALL");
  rtcm_.time = gnss_common::doubleToGtime(option.start_time);
  for (int i = 0; i < MaxDataSize::RTCM3; i++) {
    data_.push_back(std::make_shared<DataCluster>(type_));
  }
}

RTCM3Formator::RTCM3Formator(const YAML::Node& node)
{
  Option option;
  option.start_time = vk::Timer::getCurrentTime();
  loadStartTime(node, option.start_time);

  type_ = FormatorType::RTCM3;

  memset(&rtcm_, 0, sizeof(rtcm_t));
  init_rtcm(&rtcm_);
  sprintf(rtcm_.opt, "-EPHALL");
  rtcm_.time = gnss_common::doubleToGtime(option.start_time);
  for (int i = 0; i < MaxDataSize::RTCM3; i++) {
    data_.push_back(std::make_shared<DataCluster>(type_));
  }
}

RTCM3Formator::~RTCM3Formator()
{
  free_rtcm(&rtcm_);
}

// Decode stream to data
int RTCM3Formator::decode(const uint8_t *buf, int size, 
    std::vector<std::shared_ptr<DataCluster>>& data)
{
  // Clear old informations and get GNSS data handle
  std::vector<std::shared_ptr<DataCluster::GNSS>> gnss_data;
  for (size_t i = 0; i < data_.size(); i++) {
    data_[i]->gnss->types.clear();
    gnss_data.push_back(data_[i]->gnss);
  }

  bool is_observation = false;
  bool has_others = false;
  int iobs = 0;
  for (int i = 0; i < size; i++) {
    int ret = input_rtcm3(&rtcm_, buf[i]);
    if (ret <= 0) continue;

    obs_t *obs = &rtcm_.obs;
    nav_t *nav = &rtcm_.nav;
    sta_t *sta = &rtcm_.sta;
    ssr_t *ssr = rtcm_.ssr;
    int sat = rtcm_.ephsat;
    gnss_common::updateStreamData(
        ret, obs, nav, sta, ssr, iobs, sat, gnss_data);
    GnssDataType type = static_cast<GnssDataType>(ret);
    std::shared_ptr<DataCluster::GNSS>& gnss = 
      type == GnssDataType::Observation ? gnss_data[iobs] : gnss_data[0];
    if (std::find(gnss->types.begin(), gnss->types.end(), type)
      == gnss->types.end()) {
      gnss->types.push_back(type);
    }

    if (type == GnssDataType::Observation) {
      if (iobs < MaxDataSize::RTCM3) iobs++;
      if (iobs >= MaxDataSize::RTCM3) {
        LOG(WARNING) << "Max data length surpassed!";
        break;
      }
      is_observation = true;
    }
    else {
      has_others = true;
    }
  }

  data = data_;

  return is_observation ? iobs : has_others;
}

// Encode data to stream
int RTCM3Formator::encode(
  const std::shared_ptr<DataCluster>& data, uint8_t *buf)
{
#if 0  // not finished yet
  // Check the control structure
  std::map<GnssDataType, bool> type_valid;
  type_valid.insert(std::make_pair(GnssDataType::Observation, false));
  type_valid.insert(std::make_pair(GnssDataType::Ephemeris, false));
  type_valid.insert(std::make_pair(GnssDataType::AntePos, false));
  type_valid.insert(std::make_pair(GnssDataType::IonAndUtcPara, false));
  type_valid.insert(std::make_pair(GnssDataType::SSR, false));
  type_valid.insert(std::make_pair(GnssDataType::PhaseCenter, false));
  for (auto it : data->gnss->types) {
    type_valid.at(it) = true;
  }

  // Encode data
  std::vector<int> msg_obs = {1077, 1087, 1097, 1127};
  std::vector<int> msg_eph = {1019, 1020, 1045, 1046, 1042};
  std::vector<int> msg_ant = {1005};
  int n = 0;
  rtcm_t rtcm = rtcm_;
  memcpy(&rtcm.obs, data->gnss->observation, sizeof(obs_t));
  memcpy(&rtcm.nav, data->gnss->ephemeris, sizeof(nav_t));
  memcpy(&rtcm.sta, data->gnss->antenna, sizeof(sta_t));
  memcpy(rtcm.ssr, data->gnss->ephemeris->ssr, sizeof(ssr_t) * MAXSAT);

  if (type_valid.at(GnssDataType::Observation)) {
    // Set time
    rtcm.time = rtcm.obs.data[0].time;
    if (fabs(timediff(rtcm.time, rtcm_.time)) > 30.0) rtcm_.time = rtcm.time;

    for (size_t i = 0; i < msg_obs.size(); i++) {
      gen_rtcm3(&rtcm, msg_obs[i], 0, i != 3);
      memcpy(buf+n, rtcm.buff, rtcm.nbyte);
      n += rtcm.nbyte;
    }
  }
  if (type_valid.at(GnssDataType::Ephemeris)) {
    for (size_t i = 0; i < msg_eph.size(); i++) {
      gen_rtcm3(&rtcm, msg_eph[i], 0, i != 4);
      memcpy(buf+n, rtcm.buff, rtcm.nbyte);
      n += rtcm.nbyte;
    }
  }
  if (type_valid.at(GnssDataType::AntePos)) {
    for (size_t i = 0; i < msg_ant.size(); i++) {
      gen_rtcm3(&rtcm, msg_ant[i], 0, i != 3);
      memcpy(buf+n, rtcm.buff, rtcm.nbyte);
      n += rtcm.nbyte;
    }
  }

  return n;
#else
  LOG(ERROR) << "RTCM3 Encoding not supported!";
  return 0;
#endif
}

// GNSS raw --------------------------------------------------------
GnssRawFormator::GnssRawFormator(const Option& option)
{
  type_ = FormatorType::GnssRaw;
  option_tools::convert(option.sub_type, format_);

  init_raw(&raw_, static_cast<int>(format_));
  raw_.time = gnss_common::doubleToGtime(option.start_time);
  for (int i = 0; i < MaxDataSize::GnssRaw; i++) {
    data_.push_back(std::make_shared<DataCluster>(type_));
  }
}
  
GnssRawFormator::GnssRawFormator(const YAML::Node& node)
{
  Option option;
  option.start_time = vk::Timer::getCurrentTime();
  loadStartTime(node, option.start_time);
  LOAD_REQUIRED(sub_type);

  type_ = FormatorType::GnssRaw;
  option_tools::convert(option.sub_type, format_);

  init_raw(&raw_, static_cast<int>(format_));
  raw_.time = gnss_common::doubleToGtime(option.start_time);
  for (int i = 0; i < MaxDataSize::GnssRaw; i++) {
    data_.push_back(std::make_shared<DataCluster>(type_));
  }
} 

GnssRawFormator::~GnssRawFormator()
{
  free_raw(&raw_);
}

// Decode stream to data
int GnssRawFormator::decode(const uint8_t *buf, int size, 
    std::vector<std::shared_ptr<DataCluster>>& data)
{
  // Clear old informations and get GNSS data handle
  std::vector<std::shared_ptr<DataCluster::GNSS>> gnss_data;
  for (size_t i = 0; i < data_.size(); i++) {
    data_[i]->gnss->types.clear();
    gnss_data.push_back(data_[i]->gnss);
  }

  bool is_observation = false;
  bool has_others = false;
  int iobs = 0;
  for (int i = 0; i < size; i++) {
    int ret = input_raw(&raw_, static_cast<int>(format_), buf[i]);
    if (ret <= 0) continue;

    obs_t *obs = &raw_.obs;
    nav_t *nav = &raw_.nav;
    sta_t *sta = &raw_.sta;
    int sat = raw_.ephsat;
    gnss_common::updateStreamData(
        ret, obs, nav, sta, NULL, iobs, sat, gnss_data);
    GnssDataType type = static_cast<GnssDataType>(ret);
    std::shared_ptr<DataCluster::GNSS>& gnss = 
      type == GnssDataType::Observation ? gnss_data[iobs] : gnss_data[0];
    if (std::find(gnss->types.begin(), gnss->types.end(), type)
      == gnss->types.end()) {
      gnss->types.push_back(type);
    }
    
    if (type == GnssDataType::Observation) {
      if (iobs < MaxDataSize::GnssRaw) iobs++;
      if (iobs >= MaxDataSize::GnssRaw) {
        LOG(WARNING) << "Max data length surpassed!";
        break;
      }
      is_observation = true;

      // modify some values
      if (format_ == GnssRawFormats::Tersus || format_ == GnssRawFormats::Novatel) {
        obs_t *obs = gnss->observation;
        for (int i = 0; i < obs->n; i++) {
          for (int j = 0; j < NFREQ+NEXOBS; j++) {
            obs->data[i].D[j] = -obs->data[i].D[j];
          }
        }
      }
    }
    else {
      has_others = true;
    }
  }

  data = data_;

  return is_observation ? iobs : has_others;
}

// Encode data to stream
int GnssRawFormator::encode(const std::shared_ptr<DataCluster>& data, uint8_t *buf)
{
  LOG(ERROR) << "GNSS-Raw Encoding not supported!";
  return 0;
}

// GNSS Rinex --------------------------------------------------------
RINEXFormator::RINEXFormator(const Option& option)
{
  type_ = FormatorType::RINEX;
  system_exclude_ = option.system_exclude;
  
  if (!(buf_memory_ = (char *)malloc(sizeof(char) * option.buffer_length))) {
    LOG(FATAL) << __FUNCTION__ << ": Buffer malloc error!";
    return;
  }
  p_memory_ = buf_memory_;
  fp_memory_ = fmemopen(buf_memory_, option.buffer_length, "w+");

  init_rnxctr(&rnx_);
  for (int i = 0; i < MaxDataSize::GnssRaw; i++) {
    data_.push_back(std::make_shared<DataCluster>(type_));
  }
}
  
RINEXFormator::RINEXFormator(const YAML::Node& node)
{
  type_ = FormatorType::RINEX;

  Option option;
  LOAD_COMMON(buffer_length);
  std::vector<std::string> system_excludes;
  if (option_tools::safeGet(node, "system_exclude", &system_excludes)) {
    for (const auto& system_exclude : system_excludes) {
      if (!system_exclude.empty()) {
        option.system_exclude.push_back(system_exclude[0]);
      }
    }
  }
  if (!(buf_memory_ = (char *)malloc(sizeof(char) * option.buffer_length))) {
    LOG(FATAL) << __FUNCTION__ << ": Buffer malloc error!";
    return;
  }
  system_exclude_ = option.system_exclude;
  p_memory_ = buf_memory_;
  fp_memory_ = fmemopen(buf_memory_, option.buffer_length, "r");

  init_rnxctr(&rnx_);
  for (int i = 0; i < MaxDataSize::GnssRaw; i++) {
    data_.push_back(std::make_shared<DataCluster>(type_));
  }
} 

RINEXFormator::~RINEXFormator()
{
  free_rnxctr(&rnx_);
  fclose(fp_memory_);
  free(buf_memory_);
}

// Decode stream to data
int RINEXFormator::decode(const uint8_t *buf, int size, 
    std::vector<std::shared_ptr<DataCluster>>& data)
{
  auto filter_observations = [this](obs_t* obs) {
    int kept = 0;
    for (int i = 0; i < obs->n; i++) {
      int prn = 0;
      char system = 0;
      switch (satsys(obs->data[i].sat, &prn)) {
        case SYS_GPS: system = 'G'; break;
        case SYS_GLO: system = 'R'; break;
        case SYS_GAL: system = 'E'; break;
        case SYS_CMP: system = 'C'; break;
        case SYS_QZS: system = 'J'; break;
        default: system = 0; break;
      }
      if (system != 0 && !gnss_common::useSystem(system_exclude_, system)) {
        continue;
      }
      if (kept != i) {
        obs->data[kept] = obs->data[i];
      }
      kept++;
    }
    obs->n = kept;
  };

  // Clear old informations and get GNSS data handle
  std::vector<std::shared_ptr<DataCluster::GNSS>> gnss_data;
  for (size_t i = 0; i < data_.size(); i++) {
    data_[i]->gnss->types.clear();
    gnss_data.push_back(data_[i]->gnss);
  }

  bool is_observation = false;
  bool has_others = false;
  int iobs = 0;
  obs_t *obs = &rnx_.obs;
  nav_t *nav = &rnx_.nav;
  sta_t *sta = &rnx_.sta;
  int *sat = &rnx_.ephsat;
  eph_t eph = {0};
  geph_t geph = {0};
  seph_t seph = {0};
  for (int i = 0; i < size; i++) {
    // form a line
    line_.push_back(buf[i]);
    if (buf[i] != '\n') continue;

    // add to memory
    p_memory_ += sprintf(p_memory_, "%s", line_.data());
    *p_memory_ = 0x00;
    line_.clear();

    // decode header
    char buff[1024];
    if (!header_decoded_)
    {
      if (readrnxh(fp_memory_, &rnx_.ver, &rnx_.type, 
          &rnx_.sys, &rnx_.tsys, rnx_.tobs, &rnx_.nav, &rnx_.sta)) {
        // handle new header data
        int ret[] = {static_cast<int>(GnssDataType::IonAndUtcPara), 
                     static_cast<int>(GnssDataType::AntePos)};
        gnss_common::updateStreamData(
          ret[0], obs, nav, sta, NULL, iobs, *sat, gnss_data);
        gnss_common::updateStreamData(
          ret[1], obs, nav, sta, NULL, iobs, *sat, gnss_data);
        GnssDataType type[] = {static_cast<GnssDataType>(ret[0]), 
                               static_cast<GnssDataType>(ret[1])};
        std::shared_ptr<DataCluster::GNSS>& gnss =  gnss_data[0];
        for (size_t k = 0; k < 2; k++)
        if (std::find(gnss->types.begin(), gnss->types.end(), type[k])
          == gnss->types.end()) {
          gnss->types.push_back(type[k]);
        }
        has_others = true;
        //set flag
        header_decoded_ = true;
        // clear memory
        memset(buf_memory_, 0x00, p_memory_ - buf_memory_);
        p_memory_ = buf_memory_;
      }
      rewind(fp_memory_);
      continue;
    }

    // decode observation body
    if (rnx_.type=='O') {
      int n, flag;
      if ((n = readrnxobsb(fp_memory_, rnx_.opt, rnx_.ver,
        &rnx_.tsys, rnx_.tobs, &flag, rnx_.obs.data, &rnx_.sta)) > 0) {
        rnx_.obs.n = n;
        filter_observations(&rnx_.obs);
        if (rnx_.obs.n <= 0) {
          memset(buf_memory_, 0x00, p_memory_ - buf_memory_);
          p_memory_ = buf_memory_;
          rewind(fp_memory_);
          continue;
        }
        // handle new epoch data
        rnx_.time = rnx_.obs.data[0].time;
        int ret = static_cast<int>(GnssDataType::Observation);
        gnss_common::updateStreamData(
          ret, obs, nav, sta, NULL, iobs, *sat, gnss_data);
        GnssDataType type = static_cast<GnssDataType>(ret);
        std::shared_ptr<DataCluster::GNSS>& gnss = gnss_data[iobs];
        if (std::find(gnss->types.begin(), gnss->types.end(), type)
          == gnss->types.end()) {
          gnss->types.push_back(type);
        }
        if (iobs < MaxDataSize::RINEX) iobs++;
        if (iobs >= MaxDataSize::RINEX) {
          LOG(WARNING) << "Max data length surpassed!";
          break;
        }
        is_observation = true;

        obs_t *obs = gnss->observation;
        for (int i = 0; i < obs->n; i++) {
          for (int j = 0; j < NFREQ+NEXOBS; j++) {
            obs->data[i].D[j] = -obs->data[i].D[j];
          }
        }

        // clear memory
        memset(buf_memory_, 0x00, p_memory_ - buf_memory_);
        p_memory_ = buf_memory_;
      }
      rewind(fp_memory_);
      continue;
    }
    // decode ephemeris body
    else {
      int sys, stat, type, prn, set;
      if ((stat = readrnxnavb(fp_memory_, rnx_.opt, 
        rnx_.ver, rnx_.sys, &type, &eph, &geph, &seph)) > 0) {
        // handle new ephemeris data
        bool valid = true;
        char system = 0;
        if (type == 1) { /* GLONASS ephemeris */
          sys = satsys(geph.sat, &prn);
          if (sys == SYS_GLO) system = 'R';
          if (system != 0 && !gnss_common::useSystem(system_exclude_, system)) {
            memset(buf_memory_, 0x00, p_memory_ - buf_memory_);
            p_memory_ = buf_memory_;
            rewind(fp_memory_);
            continue;
          }
          rnx_.nav.geph[prn-1] = geph;
          rnx_.time = geph.tof;
          rnx_.ephsat = geph.sat;
          rnx_.ephset = 0;
        }
        else { /* other ephemeris */
          sys = satsys(eph.sat, &prn);
          switch (sys) {
            case SYS_GPS: system = 'G'; break;
            case SYS_GAL: system = 'E'; break;
            case SYS_CMP: system = 'C'; break;
            case SYS_QZS: system = 'J'; break;
            default: system = 0; break;
          }
          if (system != 0 && !gnss_common::useSystem(system_exclude_, system)) {
            memset(buf_memory_, 0x00, p_memory_ - buf_memory_);
            p_memory_ = buf_memory_;
            rewind(fp_memory_);
            continue;
          }
          if (sys == SYS_GAL) {
            int sel = getseleph(sys);
            if (sel == 0 && !(eph.code&(1<<9))) valid = false;
            if (sel == 1 && !(eph.code&(1<<8))) valid = false;
          }
          if (valid) {
            rnx_.nav.eph[eph.sat-1] = eph;
            rnx_.time = eph.ttr;
            rnx_.ephsat = eph.sat;
            rnx_.ephset = 0;
          }
          else {
            // clear memory and continue
            memset(buf_memory_, 0x00, p_memory_ - buf_memory_);
            p_memory_ = buf_memory_;
            continue;
          }
        }

        if (valid)
        {
          int ret = static_cast<int>(GnssDataType::Ephemeris);
          gnss_common::updateStreamData(
            ret, obs, nav, sta, NULL, iobs, *sat, gnss_data);
          GnssDataType type = static_cast<GnssDataType>(ret);
          std::shared_ptr<DataCluster::GNSS>& gnss =  gnss_data[0];
          if (std::find(gnss->types.begin(), gnss->types.end(), type)
            == gnss->types.end()) {
            gnss->types.push_back(type);
          }
          has_others = true;
        }

        // clear memory
        memset(buf_memory_, 0x00, p_memory_ - buf_memory_);
        p_memory_ = buf_memory_;
      }
      rewind(fp_memory_);
      continue;
    }
  }

  data = data_;

  return is_observation ? iobs : has_others;
}

// Encode data to stream
int RINEXFormator::encode(const std::shared_ptr<DataCluster>& data, uint8_t *buf)
{
  LOG(ERROR) << "GNSS-Rinex Encoding not supported!";
  return 0;
}

// Image V4L2 ------------------------------------------------
ImageV4L2Formator::ImageV4L2Formator(const Option& option)
{
  type_ = FormatorType::ImageV4L2;

  init_img(&image_, option.width, option.height, option.step);
  data_.push_back(std::make_shared<DataCluster>(
    FormatorType::ImageV4L2, option.width, option.height, option.step));
}

ImageV4L2Formator::ImageV4L2Formator(const YAML::Node& node)
{
  Option option;
  LOAD_REQUIRED(width);
  LOAD_REQUIRED(height);
  LOAD_COMMON(step);

  type_ = FormatorType::ImageV4L2;

  init_img(&image_, option.width, option.height, option.step);
  data_.push_back(std::make_shared<DataCluster>(
    FormatorType::ImageV4L2, option.width, option.height, option.step));
}

ImageV4L2Formator::~ImageV4L2Formator()
{
  free_img(&image_);
}

// Decode stream to data
int ImageV4L2Formator::decode(const uint8_t *buf, int size, 
    std::vector<std::shared_ptr<DataCluster>>& data)
{
  int ret = input_image_v4l2(&image_, buf, size);
  if (ret <= 0) return 0;

  memcpy(data_[0]->image->image, image_.image, 
    sizeof(uint8_t) * image_.width * image_.height * image_.step);
  data = data_;

  return 1;
}

// Encode data to stream
int ImageV4L2Formator::encode(const std::shared_ptr<DataCluster>& data, uint8_t *buf)
{
  LOG(ERROR) << "Image-V4L2 Encoding not supported!";
  return 0;
}
  
// Image pack -------------------------------------------------
ImagePackFormator::ImagePackFormator(const Option& option)
{
  type_ = FormatorType::ImagePack;

  init_img(&image_, option.width, option.height, option.step);
  for (int i = 0; i < MaxDataSize::ImagePack; i++) {
    data_.push_back(std::make_shared<DataCluster>(
      FormatorType::ImagePack, option.width, option.height, option.step));
  }
}

ImagePackFormator::ImagePackFormator(const YAML::Node& node)
{
  Option option;
  LOAD_REQUIRED(width);
  LOAD_REQUIRED(height);
  LOAD_COMMON(step);

  type_ = FormatorType::ImagePack;

  init_img(&image_, option.width, option.height, option.step);
  for (int i = 0; i < MaxDataSize::ImagePack; i++) {
    data_.push_back(std::make_shared<DataCluster>(
      FormatorType::ImagePack, option.width, option.height, option.step));
  }
}

ImagePackFormator::~ImagePackFormator()
{
  free_img(&image_);
}

// Decode stream to data
int ImagePackFormator::decode(const uint8_t *buf, int size, 
    std::vector<std::shared_ptr<DataCluster>>& data)
{
  int iobs = 0;
  for (int i = 0; i < size; i++) {
    int ret = input_image(&image_, buf[i]);
    if (ret <= 0) continue;

    data_[iobs]->image->time = gnss_common::gtimeToDouble(image_.time);
    // TODO: This memcpy may increase memory occupation, according to
    // https://bbs.csdn.net/topics/390705325, maybe it is caused by the
    // convertion between physical and vitural memory, and this is not 
    // a memory leak.
    memcpy(data_[iobs]->image->image, image_.image, 
      sizeof(uint8_t) * image_.width * image_.height * image_.step);

    if (++iobs >= MaxDataSize::ImagePack) {
      LOG(WARNING) << "Max data length surpassed!";
      break;
    }
  }

  data = data_;

  return iobs;
}

// Encode data to stream
int ImagePackFormator::encode(
    const std::shared_ptr<DataCluster>& data, uint8_t *buf)
{
  img_t *image;
  init_img(image, data->image->width, data->image->height, data->image->step);
  image->time = gnss_common::doubleToGtime(data->image->time);
  memcpy(image->image, data->image->image, 
       data->image->width * data->image->height * data->image->step);

  if (!gen_img(image)) return 0;

  memcpy(buf, image->buff, image->nbyte);
  int nbyte = image->nbyte;
  free_img(image);

  return nbyte;
}

// IMU pack --------------------------------------------------
IMUPackFormator::IMUPackFormator(const Option& option)
{
  type_ = FormatorType::IMUPack;

  init_imu(&imu_);
}

IMUPackFormator::IMUPackFormator(const YAML::Node& node)
{
  type_ = FormatorType::IMUPack;

  init_imu(&imu_);
}

IMUPackFormator::~IMUPackFormator()
{
  free_imu(&imu_);
}

// Decode stream to data
int IMUPackFormator::decode(const uint8_t *buf, int size, 
    std::vector<std::shared_ptr<DataCluster>>& data)
{
  int n_data = 0;
  data.clear();
  for (int i = 0; i < size; i++) {
    int ret = input_imu(&imu_, buf[i]);
    if (ret <= 0) continue;

    std::shared_ptr<DataCluster> data_ptr;
    data_ptr = std::make_shared<DataCluster>(FormatorType::IMUPack);
    data_ptr->imu->time = gnss_common::gtimeToDouble(imu_.time);
    for (int k = 0; k < 3; k++) {
      data_ptr->imu->acceleration[k] = imu_.acc[k];
      data_ptr->imu->angular_velocity[k] = imu_.gyro[k];
    }
    data.push_back(data_ptr);

    if (++n_data >= MaxDataSize::IMUPack) {
      LOG(WARNING) << "Max data length surpassed!";
      break;
    }
  }

  return n_data;
}

// Encode data to stream
int IMUPackFormator::encode(const std::shared_ptr<DataCluster>& data, uint8_t *buf)
{
  imu_t *imu;
  init_imu(imu);
  imu->time = gnss_common::doubleToGtime(data->imu->time);
  for (int i = 0; i < 3; i++) {
    imu->acc[i] = data->imu->acceleration[i];
    imu->gyro[i] = data->imu->angular_velocity[i];
  }

  if (!gen_imu(imu)) return 0;

  memcpy(buf, imu->buff, imu->nbyte);
  int nbyte = imu->nbyte;
  free_imu(imu);
  return imu->nbyte;
}

// Option pack --------------------------------------------------
OptionFormator::OptionFormator(const Option& option)
{

}

OptionFormator::OptionFormator(const YAML::Node& node)
{

}

OptionFormator::~OptionFormator()
{

}

// Decode stream to data
int OptionFormator::decode(const uint8_t *buf, int size, 
    std::vector<std::shared_ptr<DataCluster>>& data)
{
  return 0;
}

// Encode data to stream
int OptionFormator::encode(const std::shared_ptr<DataCluster>& data, uint8_t *buf)
{
  return 0;
}

// NMEA ----------------------------------------------------------
NmeaFormator::NmeaFormator(const Option& option)
{
  type_ = FormatorType::NMEA;

  option_ = option;
}

NmeaFormator::NmeaFormator(const YAML::Node& node)
{
  type_ = FormatorType::NMEA;

  Option option;
  LOAD_COMMON(use_gga);
  LOAD_COMMON(use_rmc);
  LOAD_COMMON(use_esa);
  LOAD_COMMON(use_esd);
  LOAD_COMMON(talker_id);
  option_ = option;
}

NmeaFormator::~NmeaFormator()
{

}

// Decode stream to data
int NmeaFormator::decode(const uint8_t *buf, int size, 
    std::vector<std::shared_ptr<DataCluster>>& data)
{
  LOG(ERROR) << "NMEA decoding not supported!";

  return 0;
}

// Encode data to stream
int NmeaFormator::encode(const std::shared_ptr<DataCluster>& data, uint8_t *buf)
{
  if (data->solution == nullptr) return 0;

  uint8_t *p = buf;
  if (option_.use_rmc) {
    p += encodeRMC(*data->solution, p);
  }
  if (option_.use_gga) {
    p += encodeGGA(*data->solution, p);
  }
  if (option_.use_esa) {
    p += encodeESA(*data->solution, p);
  }
  if (option_.use_esd) {
    p += encodeESD(*data->solution, p);
  }

  return p - buf;
}

#define MAXFIELD   64           /* max number of fields in a record */
#define MAXNMEA    256          /* max length of nmea sentence */
#define KNOT2M     0.514444444  /* m/knot */
static const int nmea_solq[]={  /* NMEA GPS quality indicator */
    /* 0=Fix not available or invalid */
    /* 1=GPS SPS Mode, fix valid */
    /* 2=Differential GPS, SPS Mode, fix valid */
    /* 3=GPS PPS Mode, fix valid */
    /* 4=Real Time Kinematic. System used in RTK mode with fixed integers */
    /* 5=Float RTK. Satellite system used in RTK mode, floating integers */
    /* 6=Estimated (dead reckoning) Mode */
    /* 7=Manual Input Mode */
    /* 8=Simulation Mode */
    SOLQ_NONE ,SOLQ_SINGLE, SOLQ_DGPS, SOLQ_PPP , SOLQ_FIX,
    SOLQ_FLOAT,SOLQ_DR    , SOLQ_NONE, SOLQ_NONE, SOLQ_NONE
};

// Encode GNGGA message
int NmeaFormator::encodeGGA(const Solution& solution, uint8_t* buf)
{
  sol_t sol;
  convertSolution(solution, sol);

  gtime_t time;
  double h,ep[6],pos[3],dms1[3],dms2[3],dop=1.0;
  int solq,refid=0;
  char *p=(char *)buf,*q,sum;
  
  if (sol.stat<=SOLQ_NONE) {
    p+=sprintf(p,"$%sGGA,,,,,,,,,,,,,,",option_.talker_id.data());
    for (q=(char *)buf+1,sum=0;*q;q++) sum^=*q;
    p+=sprintf(p,"*%02X%c%c",sum,0x0D,0x0A);
    return p-(char *)buf;
  }
  for (solq=0;solq<8;solq++) if (nmea_solq[solq]==sol.stat) break;
  if (solq>=8) solq=0;
  time=gpst2utc(sol.time);
  time2epoch(time,ep);
  ecef2pos(sol.rr,pos);
  h=geoidh(pos);
  deg2dms(fabs(pos[0])*R2D,dms1,7);
  deg2dms(fabs(pos[1])*R2D,dms2,7);
  p+=sprintf(p,"$%sGGA,%02.0f%02.0f%06.3f,%02.0f%010.7f,%s,%03.0f%010.7f,%s,"
              "%d,%02d,%.1f,%.3f,M,%.3f,M,%.1f,%04d",
              option_.talker_id.data(),ep[3],ep[4],ep[5],dms1[0],dms1[1]+dms1[2]/60.0,
              pos[0]>=0?"N":"S",dms2[0],dms2[1]+dms2[2]/60.0,pos[1]>=0?"E":"W",
              solq,sol.ns,dop,pos[2]-h,h,sol.age,refid);
  for (q=(char *)buf+1,sum=0;*q;q++) sum^=*q; /* check-sum */
  p+=sprintf(p,"*%02X\r\n",sum);
  return p-(char *)buf;
}

// Encode GNRMC message
int NmeaFormator::encodeRMC(const Solution& solution, uint8_t* buf)
{
  sol_t sol;
  convertSolution(solution, sol);

  static double dirp=0.0;
  gtime_t time;
  double ep[6],pos[3],enuv[3],dms1[3],dms2[3],vel,dir,amag=0.0;
  char *p=(char *)buf,*q,sum;
  const char *emag="E",*mode="A",*status="V";
  
  trace(3,"outnmea_rmc:\n");
  
  if (sol.stat<=SOLQ_NONE) {
    p+=sprintf(p,"$%sRMC,,,,,,,,,,,,,",option_.talker_id.data());
    for (q=(char *)buf+1,sum=0;*q;q++) sum^=*q;
    p+=sprintf(p,"*%02X%c%c",sum,0x0D,0x0A);
    return p-(char *)buf;
  }
  time=gpst2utc(sol.time);
  time2epoch(time,ep);
  ecef2pos(sol.rr,pos);
  ecef2enu(pos,sol.rr+3,enuv);
  vel=norm(enuv,3);
  if (vel>=1.0) {
    dir=atan2(enuv[0],enuv[1])*R2D;
    if (dir<0.0) dir+=360.0;
    dirp=dir;
  }
  else {
    dir=dirp;
  }
  if      (sol.stat==SOLQ_DGPS ||sol.stat==SOLQ_SBAS) mode="D";
  else if (sol.stat==SOLQ_FLOAT||sol.stat==SOLQ_FIX ) mode="R";
  else if (sol.stat==SOLQ_PPP) mode="P";
  deg2dms(fabs(pos[0])*R2D,dms1,7);
  deg2dms(fabs(pos[1])*R2D,dms2,7);
  p+=sprintf(p,"$%sRMC,%02.0f%02.0f%06.3f,A,%02.0f%010.7f,%s,%03.0f%010.7f,"
              "%s,%4.2f,%4.2f,%02.0f%02.0f%02d,%.1f,%s,%s,%s",
              option_.talker_id.data(),ep[3],ep[4],ep[5],dms1[0],dms1[1]+dms1[2]/60.0,
              pos[0]>=0?"N":"S",dms2[0],dms2[1]+dms2[2]/60.0,pos[1]>=0?"E":"W",
              vel/KNOT2M,dir,ep[2],ep[1],(int)ep[0]%100,amag,emag,mode,status);
  for (q=(char *)buf+1,sum=0;*q;q++) sum^=*q; /* check-sum */
  p+=sprintf(p,"*%02X\r\n",sum);
  return p-(char *)buf;
}

// Encode GNESA (self-defined Extended Speed and Attitude) message
// Format: $GNESA,tod,Ve,Vn,Vu,Ar,Ap,Ay*checksum
int NmeaFormator::encodeESA(const Solution& solution, uint8_t* buf)
{
  sol_t sol;
  convertSolution(solution, sol);
  Eigen::Vector3d rpy = quaternionToEulerAngle(solution.pose.getEigenQuaternion());
  rpy *= R2D;

  gtime_t time;
  double ep[6];
  char *p=(char *)buf,*q,sum;
  
  if (sol.stat<=SOLQ_NONE) {
    p+=sprintf(p,"$%sESA,,,,,,,",option_.talker_id.data());
    for (q=(char *)buf+1,sum=0;*q;q++) sum^=*q;
    p+=sprintf(p,"*%02X%c%c",sum,0x0D,0x0A);
    return p-(char *)buf;
  }
  time=gpst2utc(sol.time);
  time2epoch(time,ep);
  p+=sprintf(p,"$%sESA,%02.0f%02.0f%06.3f,%+.3f,%+.3f,%+.3f,"
             "%+.3f,%+.3f,%+.3f",
             option_.talker_id.data(),ep[3],ep[4],ep[5],sol.rr[3],sol.rr[4],sol.rr[5],
             rpy[0],rpy[1],rpy[2]);
  for (q=(char *)buf+1,sum=0;*q;q++) sum^=*q; /* check-sum */
  p+=sprintf(p,"*%02X\r\n",sum);
  return p-(char *)buf;
}

// Encode GNESD (self-defined Extended STD) message
// Format: $GNESD,tod,STD_Pe,STD_Pn,STD_Pu,STD_Ve,STD_Vn,STD_Vu,
//         STD_Ar,STD_Ap,STD_Py*checksum
int NmeaFormator::encodeESD(const Solution& solution, uint8_t* buf)
{
  sol_t sol;
  convertSolution(solution, sol);

  gtime_t time;
  double ep[6], std_p[3], std_v[3], std_a[3];
  char *p=(char *)buf,*q,sum;

  for (size_t i = 0; i < 3; i++) {
    std_p[i] = sqrt(solution.covariance(i, i));
    std_v[i] = sqrt(solution.covariance(i + 6, i + 6));
    std_a[i] = sqrt(solution.covariance(i + 3, i + 3)) * R2D;
  }
  
  if (sol.stat<=SOLQ_NONE) {
    p+=sprintf(p,"$%sESD,,,,,,,",option_.talker_id.data());
    for (q=(char *)buf+1,sum=0;*q;q++) sum^=*q;
    p+=sprintf(p,"*%02X%c%c",sum,0x0D,0x0A);
    return p-(char *)buf;
  }
  time=gpst2utc(sol.time);
  time2epoch(time,ep);
  p+=sprintf(p,"$%sESD,%02.0f%02.0f%06.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f",
             option_.talker_id.data(),ep[3],ep[4],ep[5],std_p[0],std_p[1],std_p[2],
             std_v[0],std_v[1],std_v[2],std_a[0],std_a[1],std_a[2]);
  for (q=(char *)buf+1,sum=0;*q;q++) sum^=*q; /* check-sum */
  p+=sprintf(p,"*%02X\r\n",sum);
  return p-(char *)buf;
}

// Convert Solution to sol_t
void NmeaFormator::convertSolution(const Solution& solution, sol_t& sol)
{
  sol.type = 0;
  sol.time = gnss_common::doubleToGtime(solution.timestamp);
  sol.time = utc2gpst(sol.time);
  sol.age = solution.differential_age;
  sol.ns = solution.num_satellites;
  Eigen::Map<Eigen::Vector3d> rr(sol.rr);
  rr = solution.coordinate->convert(
    solution.pose.getPosition(), GeoType::ENU, GeoType::ECEF);
  Eigen::Map<Eigen::Vector3d> vv(sol.rr + 3);
  vv = solution.speed_and_bias.head<3>();
  if (solution.status == GnssSolutionStatus::Fixed) sol.stat = SOLQ_FIX;
  else if (solution.status == GnssSolutionStatus::Float) sol.stat = SOLQ_FLOAT;
  else if (solution.status == GnssSolutionStatus::DGNSS) sol.stat = SOLQ_DGPS;
  else if (solution.status == GnssSolutionStatus::Single) sol.stat = SOLQ_SINGLE;
  else if (solution.status == GnssSolutionStatus::DeadReckoning) sol.stat = SOLQ_DR;
  else sol.stat = SOLQ_NONE;
}

// DCB file pack --------------------------------------------------
DcbFileFormator::DcbFileFormator(const Option& option)
{
  type_ = FormatorType::DcbFile;
  line_.reserve(max_line_length_); 
}

DcbFileFormator::DcbFileFormator(const YAML::Node& node)
{
  type_ = FormatorType::DcbFile;
  line_.reserve(max_line_length_); 
}

DcbFileFormator::~DcbFileFormator()
{

}

// Decode stream to data
int DcbFileFormator::decode(const uint8_t *buf, int size, 
    std::vector<std::shared_ptr<DataCluster>>& data)
{
  auto parseLine = [&](const std::string& line) {
    if (line.empty()) return;

    if (line.substr(0, 5) == "*BIAS") {
      passed_header_ = true;
      return;
    }
    if (!passed_header_) return;
    if (line == "-BIAS/SOLUTION") {
      finished_reading_ = true;
      return;
    }

    std::istringstream iss(line);
    std::vector<std::string> fields;
    std::string field;
    while (iss >> field) {
      fields.push_back(field);
    }
    if (fields.size() < 9) return;

    const std::string& bias_type = fields[0];
    const std::string& str_prn = fields[2];
    if (str_prn.size() != 3) return;

    const char system = str_prn[0];
    if (std::find(getGnssSystemList().begin(), getGnssSystemList().end(),
        system) == getGnssSystemList().end()) {
      return;
    }

    if (bias_type == "DSB") {
      if (fields.size() < 10) return;
      const std::string& str_obs1 = fields[3];
      const std::string& str_obs2 = fields[4];
      if (str_obs1.size() != 3 || str_obs2.size() != 3 ||
          str_obs1[0] != 'C' || str_obs2[0] != 'C') {
        return;
      }

      Dcb dcb;
      dcb.code1 = gnss_common::rinexTypeToCodeType(system, str_obs1.substr(1, 2));
      dcb.code2 = gnss_common::rinexTypeToCodeType(system, str_obs2.substr(1, 2));
      dcb.value = atof(fields[8].data()) * 1e-9 * CLIGHT; // ns to m
      dcb.std = atof(fields[9].data()) * 1e-9 * CLIGHT;
      if (dcb.code1 > 0 && dcb.code2 > 0) {
        dcbs_.insert(std::make_pair(str_prn, dcb));
      }
      return;
    }

    if (bias_type == "OSB") {
      if (line.size() < 82) return;
      const std::string station = line.substr(15, 9);
      if (station.find_first_not_of(' ') != std::string::npos) {
        return;
      }
      const std::string obs = line.substr(25, 4);
      const size_t obs_end = obs.find_last_not_of(' ');
      if (obs_end == std::string::npos) return;
      const std::string obs_type = obs.substr(0, obs_end + 1);
      if (obs_type.size() != 3) return;
      OsbData osb;
      osb.time1 = line.substr(35, 14);
      osb.time2 = line.substr(50, 14);
      osb.value = atof(line.substr(70, 21).c_str());
      osb.std = atof(line.substr(92).c_str());

      const int prn_number = atoi(str_prn.substr(1).c_str());
      if (obs_type[0] == 'C') {
        code_osbs_[system][prn_number][obs_type] = osb;
      }
      else if (obs_type[0] == 'L') {
        phase_osbs_[system][prn_number][obs_type] = osb;
      }
    }
  };

  if (finished_reading_) return 0;

  for (int i = 0; i < size; i++) {
    if (buf[i] != '\n') {
      if ((int)line_.size() < max_line_length_) line_ = line_ + (char)buf[i];
    }
    else {
      parseLine(line_);
      line_.clear();
      if (finished_reading_) break;
    }
  }

  if (size == 0 && !line_.empty()) {
    parseLine(line_);
    line_.clear();
  }
  if (size == 0 && passed_header_ &&
      (!dcbs_.empty() || !code_osbs_.empty() || !phase_osbs_.empty())) {
    finished_reading_ = true;
  }

  // convert to DataCluster
  if (finished_reading_) {
    std::shared_ptr<DataCluster> data_cluster = 
      std::make_shared<DataCluster>(FormatorType::DcbFile);
    nav_t* ephemeris = data_cluster->gnss->ephemeris;
    size_t num_code_osb_entries = 0;
    size_t num_phase_osb_entries = 0;
    size_t num_code_osb_sats = 0;
    size_t num_phase_osb_sats = 0;
    size_t num_phase_same_band_conflicts = 0;
    std::vector<std::string> phase_same_band_conflict_samples;
    for (const auto& system_entry : code_osbs_) {
      num_code_osb_sats += system_entry.second.size();
      for (const auto& satellite_entry : system_entry.second) {
        num_code_osb_entries += satellite_entry.second.size();
      }
    }
    for (const auto& system_entry : phase_osbs_) {
      num_phase_osb_sats += system_entry.second.size();
      for (const auto& satellite_entry : system_entry.second) {
        num_phase_osb_entries += satellite_entry.second.size();
      }
    }

    for (const auto& system_entry : phase_osbs_) {
      char system = system_entry.first;
      for (const auto& satellite_entry : system_entry.second) {
        int prn_number = satellite_entry.first;
        std::map<int, std::vector<std::pair<std::string, double>>> phase_to_values;
        for (const auto& code_entry : satellite_entry.second) {
          const std::string& obs_name = code_entry.first;
          if (obs_name.size() != 3) continue;
          int code = gnss_common::rinexTypeToCodeType(system, obs_name.substr(1, 2));
          if (code <= 0) continue;
          int phase_id = gnss_common::getPhaseID(system, code);
          if (phase_id == 0) continue;
          phase_to_values[phase_id].push_back(std::make_pair(
            obs_name, code_entry.second.value * 1e-9 * CLIGHT));
        }

        for (const auto& phase_entry : phase_to_values) {
          const auto& values = phase_entry.second;
          if (values.size() <= 1) continue;

          double min_value = values.front().second;
          double max_value = values.front().second;
          for (const auto& item : values) {
            min_value = std::min(min_value, item.second);
            max_value = std::max(max_value, item.second);
          }
          if (std::fabs(max_value - min_value) <= 1e-4) continue;

          num_phase_same_band_conflicts++;
          if (phase_same_band_conflict_samples.size() < 8) {
            std::string prn = std::string(1, system);
            if (prn_number < 10) prn += "0";
            prn += std::to_string(prn_number);
            std::ostringstream oss;
            oss << prn << ":phase=" << phase_entry.first << ":";
            for (size_t i = 0; i < values.size(); ++i) {
              if (i > 0) oss << ",";
              oss << values[i].first << "=" << values[i].second;
            }
            phase_same_band_conflict_samples.push_back(oss.str());
          }
        }
      }
    }

    for (const auto& system_entry : code_osbs_) {
      char system = system_entry.first;
      for (const auto& satellite_entry : system_entry.second) {
        int prn_number = satellite_entry.first;
        std::string prn = std::string(1, system);
        if (prn_number < 10) prn = prn + "0";
        prn = prn + std::to_string(prn_number);

        int sat = gnss_common::prnToSat(prn);
        if (sat <= 0) continue;

        for (const auto& code_entry : satellite_entry.second) {
          const std::string& code_name = code_entry.first;
          const OsbData& osb = code_entry.second;
          int code = gnss_common::rinexTypeToCodeType(system, code_name.substr(1, 2));
          if (code <= 0) continue;

          double cbias = osb.value * 1e-9 * CLIGHT;
          if (cbias == 0.0) cbias = 1e-4;
          ephemeris->ssr[sat - 1].t0[4] = timeget();
          ephemeris->ssr[sat - 1].cbias[code - 1] = cbias;
          ephemeris->ssr[sat - 1].update = 1;
        }
      }
    }

    for (const auto& system_entry : phase_osbs_) {
      char system = system_entry.first;
      for (const auto& satellite_entry : system_entry.second) {
        int prn_number = satellite_entry.first;
        std::string prn = std::string(1, system);
        if (prn_number < 10) prn = prn + "0";
        prn = prn + std::to_string(prn_number);

        int sat = gnss_common::prnToSat(prn);
        if (sat <= 0) continue;

        for (const auto& code_entry : satellite_entry.second) {
          const std::string& obs_name = code_entry.first;
          const OsbData& osb = code_entry.second;
          int code = gnss_common::rinexTypeToCodeType(system, obs_name.substr(1, 2));
          if (code <= 0) continue;

          ephemeris->ssr[sat - 1].t0[5] = timeget();
          ephemeris->ssr[sat - 1].pbias[code - 1] = osb.value * 1e-9 * CLIGHT;
          ephemeris->ssr[sat - 1].stdpb[code - 1] = osb.std * 1e-9 * CLIGHT;
          ephemeris->ssr[sat - 1].update = 1;
        }
      }
    }

    // get all PRNs
    std::vector<std::string> prns;
    for (auto dcb : dcbs_) {
      if (prns.size() == 0 || prns.back() != dcb.first) {
        prns.push_back(dcb.first);
      }
    }
    // fill DCBs of every satellites
    for (auto prn : prns) {
      std::vector<int> codes;
      std::vector<Dcb> dcbs;
      for (auto dcb = dcbs_.lower_bound(prn); 
          dcb != dcbs_.upper_bound(prn); dcb++) {
        if (std::find(codes.begin(), codes.end(), dcb->second.code1) 
            == codes.end()) { 
          codes.push_back(dcb->second.code1);
        } 
        if (std::find(codes.begin(), codes.end(), dcb->second.code2) 
            == codes.end()) { 
          codes.push_back(dcb->second.code2);
        } 
        dcbs.push_back(dcb->second);
      }

      // apply a least-square to convert DCB to code biases
      Eigen::VectorXd x = Eigen::VectorXd::Zero(codes.size());
      Eigen::VectorXd z = Eigen::VectorXd::Zero(dcbs.size() + 1);
      Eigen::MatrixXd H = Eigen::MatrixXd::Zero(dcbs.size() + 1, codes.size());
      for (size_t i = 0; i < dcbs.size(); i++) {
        z(i) = dcbs[i].value;
        for (size_t j = 0; j < codes.size(); j++) {
          if (codes[j] == dcbs[i].code1) H(i, j) = -1.0;
          if (codes[j] == dcbs[i].code2) H(i, j) = 1.0;
        }
      }
      // set the first code of each satellite as zero
      z(dcbs.size()) = 0.0;
      H(dcbs.size(), 0) = 1.0;
      // add a further contraint for Galileo: C1C = C1X
      if (prn[0] == 'E') {
        z.conservativeResize(dcbs.size() + 2);
        H.conservativeResize(dcbs.size() + 2, Eigen::NoChange);
        size_t i = dcbs.size() + 1;
        z(i) = 0.0;
        for (size_t j = 0; j < codes.size(); j++) {
          if (codes[j] == CODE_L1C) H(i, j) = 1.0;
          else if (codes[j] == CODE_L1X) H(i, j) = -1.0;
          else H(i, j) = 0.0;
        }
      }

      // Check rank
      if (checkZero((H.transpose() * H).determinant())) {
        LOG(INFO) << "Input DCBs are not closed for " << prn;
        continue;
      }

      // solve
      x = (H.transpose() * H).inverse() * H.transpose() * z;

      // fill ephemeris
      int sat = gnss_common::prnToSat(prn);
      for (size_t i = 0; i < codes.size(); i++) {
        if (sat <= 0 || codes[i] <= 0) continue;
        if (x(i) == 0.0) x(i) = 1e-4;
        ephemeris->ssr[sat - 1].t0[4] = timeget();
        ephemeris->ssr[sat - 1].cbias[codes[i] - 1] = x(i);
        ephemeris->ssr[sat - 1].update = 1;
        ephemeris->ssr[sat - 1].isdcb = 1;
      }
    }
    
	    data_cluster->gnss->types.push_back(GnssDataType::SSR);
      LOG(INFO) << "DcbFileFormator: Loaded OSB/DCB from BIA. "
                << "DSB pairs=" << dcbs_.size()
                << ", code OSB entries=" << num_code_osb_entries
                << " on " << num_code_osb_sats << " satellites"
                << ", phase OSB entries=" << num_phase_osb_entries
                << " on " << num_phase_osb_sats << " satellites"
                << ", same-band phase conflicts=" << num_phase_same_band_conflicts << ".";
      if (!phase_same_band_conflict_samples.empty()) {
        std::ostringstream oss;
        oss << "DcbFileFormator: Phase OSB differs across tracking codes on the same "
            << "carrier band. Current PPP phase-bias handling is frequency-based, "
            << "so these products may need exact code-level mapping. samples=";
        for (size_t i = 0; i < phase_same_band_conflict_samples.size(); ++i) {
          if (i > 0) oss << " | ";
          oss << phase_same_band_conflict_samples[i];
        }
        LOG(WARNING) << oss.str();
      }
	    data.clear();
	    data.push_back(data_cluster);
	    return 1;
  }

  return 0;
}

// Encode data to stream
int DcbFileFormator::encode(const std::shared_ptr<DataCluster>& data, uint8_t *buf)
{
  LOG(ERROR) << "DCB file encoding not supported!";

  return 0;
}

// ATX file pack --------------------------------------------------
AtxFileFormator::AtxFileFormator(const Option& option)
{
  type_ = FormatorType::AtxFile;
  line_.reserve(max_line_length_); 
  if (!(pcvs_ = (pcvs_t *)malloc(sizeof(pcvs_t)))) {
    free(pcvs_); return;
  }
  memset(pcvs_, 0, sizeof(pcvs_t));
}

AtxFileFormator::AtxFileFormator(const YAML::Node& node)
{
  type_ = FormatorType::AtxFile;
  line_.reserve(max_line_length_); 
  if (!(pcvs_ = (pcvs_t *)malloc(sizeof(pcvs_t)))) {
    free(pcvs_); return;
  }
  memset(pcvs_, 0, sizeof(pcvs_t));
}

AtxFileFormator::~AtxFileFormator()
{
  free(pcvs_);
}

// Decode stream to data
int AtxFileFormator::decode(const uint8_t *buf, int size, 
    std::vector<std::shared_ptr<DataCluster>>& data)
{
  if (finished_reading_) return 0;

  const pcv_t pcv0={0};
  double neu[3];
  int i,f,freq=0,freqs[]={1,2,5,0};

  for (int k = 0; k < size; k++) {
    if (buf[k] != '\n') line_ = line_ + (char)buf[k];
    // decode line
    else {
      const char *buff = line_.data();
      if (strlen(buff)<60||strstr(buff+60,"COMMENT")) {
        line_.clear(); continue;
      }
      
      if (strstr(buff+60,"START OF ANTENNA")) {
        pcv_=pcv0;
        state_=1;
      }
      if (strstr(buff+60,"END OF ANTENNA")) {
        addpcv(&pcv_,pcvs_);
        state_=0;
      }
      if (!state_) {
        line_.clear(); continue;
      }
      
      if (strstr(buff+60,"TYPE / SERIAL NO")) {
        strncpy(pcv_.type,buff   ,20); pcv_.type[20]='\0';
        strncpy(pcv_.code,buff+20,20); pcv_.code[20]='\0';
        if (!strncmp(pcv_.code+3,"        ",8)) {
            pcv_.sat=satid2no(pcv_.code);
        }
      }
      else if (strstr(buff+60,"VALID FROM")) {
        if (!str2time(buff,0,43,&pcv_.ts)) {
          line_.clear(); continue;
        }
      }
      else if (strstr(buff+60,"VALID UNTIL")) {
        if (!str2time(buff,0,43,&pcv_.te)) {
          line_.clear(); continue;
        }
      }
      else if (strstr(buff+60,"START OF FREQUENCY")) {
        // Path B: also load receiver-antenna PCV for non-GPS systems
        // (E1/E5a etc.) so that PCV correction is available for all
        // frequencies actually used in the IF combination.
        if (sscanf(buff+4,"%d",&f)<1) {
          line_.clear(); continue;
        }
        for (i=0;freqs[i];i++) if (freqs[i]==f) break;
        if (freqs[i]) freq=i+1;
      }
      else if (strstr(buff+60,"END OF FREQUENCY")) {
        freq=0;
      }
      else if (strstr(buff+60,"NORTH / EAST / UP")) {
        if (freq<1||NFREQ<freq) {
          line_.clear(); continue;
        }
        if (decodef((char *)buff,3,neu)<3) {
          line_.clear(); continue;
        }
        pcv_.off[freq-1][0]=neu[pcv_.sat?0:1]; /* x or e */
        pcv_.off[freq-1][1]=neu[pcv_.sat?1:0]; /* y or n */
        pcv_.off[freq-1][2]=neu[2];           /* z or u */
      }
      else if (strstr(buff,"NOAZI")) {
        if (freq<1||NFREQ<freq) {
          line_.clear(); continue;
        }
        if ((i=decodef((char *)(buff+8),19,pcv_.var[freq-1]))<=0) {
          line_.clear(); continue;
        }
        for (;i<19;i++) pcv_.var[freq-1][i]=pcv_.var[freq-1][i-1];
      }
      line_.clear();
    }
  }

  // convert to DataCluster
  if (last_size_ != 0 && last_size_ > size) {
    std::shared_ptr<DataCluster> data_cluster = 
      std::make_shared<DataCluster>(FormatorType::AtxFile);
    nav_t* ephemeris = data_cluster->gnss->ephemeris;

    pcv_t *pcv;
    for (i=0;i<MAXSAT;i++) {
      pcv=searchpcv(i+1,"",timeget(),pcvs_);
      ephemeris->pcvs[i]=pcv?*pcv:pcv0;
    }

    // Path B: hand the receiver-antenna entries (pcv_t::sat == 0) over to
    // the data cluster so they can be looked up by antenna type at runtime.
    pcvs_t *rcv_pcvs = (pcvs_t *)malloc(sizeof(pcvs_t));
    memset(rcv_pcvs, 0, sizeof(pcvs_t));
    int nrcv = 0;
    for (int k = 0; k < pcvs_->n; k++) {
      if (pcvs_->pcv[k].sat == 0) nrcv++;
    }
    if (nrcv > 0) {
      rcv_pcvs->pcv = (pcv_t *)malloc(sizeof(pcv_t) * nrcv);
      rcv_pcvs->nmax = nrcv;
      rcv_pcvs->n = 0;
      for (int k = 0; k < pcvs_->n; k++) {
        if (pcvs_->pcv[k].sat == 0) {
          rcv_pcvs->pcv[rcv_pcvs->n++] = pcvs_->pcv[k];
        }
      }
    }
    data_cluster->gnss->receiver_pcvs = rcv_pcvs;

    free(pcvs_->pcv);
    pcvs_->pcv = NULL;
    pcvs_->n = 0;
    pcvs_->nmax = 0;
    finished_reading_ = true;

    data_cluster->gnss->types.push_back(GnssDataType::PhaseCenter);
    data.clear();
    data.push_back(data_cluster);
    return 1;
  }
  last_size_ = size;

  return 0;
}

// Encode data to stream
int AtxFileFormator::encode(const std::shared_ptr<DataCluster>& data, uint8_t *buf)
{
  LOG(ERROR) << "ATX file encoding not supported!";

  return 0;
}

// Add antenna parameter
void AtxFileFormator::addpcv(const pcv_t *pcv, pcvs_t *pcvs)
{
  pcv_t *pcvs_pcv;
  
  if (pcvs->nmax<=pcvs->n) {
    pcvs->nmax+=256;
    if (!(pcvs_pcv=(pcv_t *)realloc(pcvs->pcv,sizeof(pcv_t)*pcvs->nmax))) {
      free(pcvs->pcv); pcvs->pcv=NULL; pcvs->n=pcvs->nmax=0;
      return;
    }
    pcvs->pcv=pcvs_pcv;
  }
  pcvs->pcv[pcvs->n++]=*pcv;
}

// Decode antenna parameter field
int AtxFileFormator::decodef(char *p, int n, double *v)
{
  int i;
  
  for (i=0;i<n;i++) v[i]=0.0;
  for (i=0,p=strtok(p," ");p&&i<n;p=strtok(NULL," ")) {
    v[i++]=atof(p)*1E-3;
  }
  return i;
}

// -------------------------------------------------------------
// Get formator handle from yaml
#define MAP_FORMATOR(Type, Formator) \
  if (type == Type) { return std::make_shared<Formator>(node); }
#define LOG_UNSUPPORT LOG(FATAL) << "Formator type not supported!";
inline static FormatorType loadType(const YAML::Node& node)
{
  if (!node["type"].IsDefined()) {
    LOG(FATAL) << "Unable to load formator type!";
  }
  std::string type_str = node["type"].as<std::string>();
  FormatorType type;
  option_tools::convert(type_str, type);
  return type;
}
std::shared_ptr<FormatorBase> makeFormator(const YAML::Node& node)
{
  FormatorType type = loadType(node);
  MAP_FORMATOR(FormatorType::RTCM2, RTCM2Formator);
  MAP_FORMATOR(FormatorType::RTCM3, RTCM3Formator);
  MAP_FORMATOR(FormatorType::GnssRaw, GnssRawFormator);
  MAP_FORMATOR(FormatorType::RINEX, RINEXFormator);
  MAP_FORMATOR(FormatorType::ImagePack, ImagePackFormator);
  MAP_FORMATOR(FormatorType::ImageV4L2, ImageV4L2Formator);
  MAP_FORMATOR(FormatorType::IMUPack, IMUPackFormator);
  MAP_FORMATOR(FormatorType::OptionPack, OptionFormator);
  MAP_FORMATOR(FormatorType::NMEA, NmeaFormator);
  MAP_FORMATOR(FormatorType::DcbFile, DcbFileFormator);
  MAP_FORMATOR(FormatorType::AtxFile, AtxFileFormator);
  // LOG_UNSUPPORT;
  return nullptr;
}

}
