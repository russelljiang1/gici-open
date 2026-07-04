/**
* @Function: Phase center handler
*
* @Author  : Cheng Chi
* @Email   : chichengcn@sjtu.edu.cn
*
* Copyright (C) 2023 by Cheng Chi, All rights reserved.
**/
#pragma once

#include <memory>
#include <map>
#include <unordered_map>
#include <string>

#include "gici/utility/rtklib_safe.h"

namespace gici {

// Phase center handle
//
// Holds:
//  * a pointer to the satellite-PCV table (sat 1..MAXSAT) that lives in
//    nav_t::pcvs (used implicitly by rtklib's satantoff/peph2pos);
//  * an optional copy of the receiver-antenna PCV that matches the rover
//    RINEX header `ANT # / TYPE` string, used to apply receiver-side
//    PCO+PCV corrections in the phase/pseudorange residuals (path B).
class PhaseCenter {
public:
  PhaseCenter(pcv_t *sat_pcvs) : sat_pcvs_(sat_pcvs), has_receiver_pcv_(false) {}
  ~PhaseCenter() {}

  // Look up the receiver antenna type in the ATX receiver-PCV table and
  // cache the matching pcv_t. Returns true if a match was found.
  // `time` is GPS time (rtklib gtime_t).
  bool setReceiverAntenna(const std::string& type,
                          const pcvs_t* receiver_pcvs,
                          gtime_t time);

  // Whether a receiver PCV has been successfully bound.
  bool hasReceiverPcv() const { return has_receiver_pcv_; }

  // Compute the LOS receiver-antenna correction (PCO projected onto LOS
  // plus elevation-dependent PCV) for the given system, RINEX code id,
  // azimuth and elevation. Returns the correction in meters that should
  // be ADDED to the modeled geometric range (i.e. subtracted from the
  // measured range). Returns 0 if no PCV is available or the frequency
  // index cannot be resolved.
  // azimuth/elevation are in radians.
  double receiverAntennaDelay(int sys, uint8_t code,
                              double azimuth, double elevation) const;

  // Raw access (kept for backward compatibility with code that only uses
  // the satellite PCV table reference).
  const pcv_t* satellitePcvs() const { return sat_pcvs_; }

private:
  pcv_t *sat_pcvs_;            // Satellite PCV table (sat 1..MAXSAT)
  pcv_t  rcv_pcv_{};           // Receiver PCV (valid iff has_receiver_pcv_)
  bool   has_receiver_pcv_;    // True once setReceiverAntenna succeeded
};

using PhaseCenterPtr = std::shared_ptr<PhaseCenter>;

} // namespace gici
