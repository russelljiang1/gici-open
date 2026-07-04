/**
* @Function: Phase center handler implementation (path B)
*
* @Author  : Cheng Chi
* @Email   : chichengcn@sjtu.edu.cn
*
* Copyright (C) 2023 by Cheng Chi, All rights reserved.
**/
#include "gici/gnss/phase_center.h"

#include <cstring>
#include <glog/logging.h>

#include "gici/utility/rtklib_safe.h"

namespace gici {

bool PhaseCenter::setReceiverAntenna(const std::string& type,
                                     const pcvs_t* receiver_pcvs,
                                     gtime_t time)
{
  has_receiver_pcv_ = false;
  if (receiver_pcvs == nullptr || receiver_pcvs->n <= 0) {
    return false;
  }

  // searchpcv expects a fixed-width antenna type string. ATX entries store
  // 20 chars antenna model + 20 chars radome (e.g. "CHCI85          NONE").
  // The RINEX header `ANT # / TYPE` stores the same convention. Pass the
  // raw string through.
  pcv_t* found = searchpcv(0, type.c_str(), time,
                           const_cast<pcvs_t*>(receiver_pcvs));
  if (found == nullptr) {
    LOG(WARNING) << "PhaseCenter: receiver antenna \"" << type
                 << "\" not found in ATX. Receiver-side PCO/PCV correction "
                    "will be skipped.";
    return false;
  }

  rcv_pcv_ = *found;
  has_receiver_pcv_ = true;
  LOG(INFO) << "PhaseCenter: bound receiver antenna \"" << type
            << "\". Receiver-side PCO/PCV correction enabled.";
  return true;
}

double PhaseCenter::receiverAntennaDelay(int sys, uint8_t code,
                                         double azimuth, double elevation) const
{
  if (!has_receiver_pcv_) return 0.0;

  // Map RINEX code -> pcv frequency-slot index used by ATX:
  //   GPS:      L1=0, L2=1, L5=2
  //   Galileo:  E1=0, E5b=1, E5a=2
  //   BDS:      B1=0, B2=1, B2a=2, ...
  //   GLONASS:  G1=0, G2=1
  // (See rtklib code2idx_*).
  int idx = code2idx(sys, code);
  if (idx < 0 || idx >= NFREQ) return 0.0;

  // rtklib's antmodel():
  //   dant[i] = -dot(off + del, e_los_ENU) + interpvar(zenith_deg, var[i])
  // and the calling convention in ppp.c:corr_meas is to SUBTRACT this from
  // the observed range (or equivalently ADD it to the modeled range).
  double azel[2] = { azimuth, elevation };
  double del[3]  = { 0.0, 0.0, 0.0 };  // RINEX antenna delta is already
                                       // applied upstream
  double dant[NFREQ] = { 0.0 };
  // opt=1 -> include zenith-angle PCV interpolation
  antmodel(const_cast<pcv_t*>(&rcv_pcv_), del, azel, 1, dant);

  return dant[idx];
}

}  // namespace gici
