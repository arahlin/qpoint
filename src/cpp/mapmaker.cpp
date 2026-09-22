#include <cmath>

#include "error.hpp"
#include "map.hpp"
#include "pointing.hpp"

namespace qp {

std::size_t num_vec(VecMode mode) {
  switch (mode) {
    case VecMode::Temp:
      return 1;
    case VecMode::D1:
    case VecMode::Pol:
      return 3;
    case VecMode::VPol:
      return 4;
    case VecMode::D2:
      return 6;
    case VecMode::D1Pol:
      return 9;
    case VecMode::D2Pol:
      return 18;
    default:
      return 0;
  }
}

std::size_t num_proj(ProjMode mode) {
  switch (mode) {
    case ProjMode::Temp:
      return 1;
    case ProjMode::Pol:
      return 6;
    case ProjMode::VPol:
      return 10;
    default:
      return 0;
  }
}

VecMode vec_mode_for(std::size_t nrow, bool pol, bool vpol) {
  switch (nrow) {
    case 1:
      return VecMode::Temp;
    case 3:
      return pol ? VecMode::Pol : VecMode::D1;
    case 4:
      return vpol ? VecMode::VPol : VecMode::None;
    // Row counts follow what map2tod1 actually reads: D2 is T with first
    // and second derivatives (6), D1Pol is (T,Q,U) with first (9).
    case 6:
      return VecMode::D2;
    case 9:
      return VecMode::D1Pol;
    case 18:
      return VecMode::D2Pol;
    default:
      return VecMode::None;
  }
}

ProjMode proj_mode_for(std::size_t nrow) {
  switch (nrow) {
    case 1:
      return ProjMode::Temp;
    case 6:
      return ProjMode::Pol;
    case 10:
      return ProjMode::VPol;
    default:
      return ProjMode::None;
  }
}

std::size_t num_proj_for_vec(std::size_t nvec) {
  return nvec ? nvec * (nvec + 1) / 2 : 0;
}

std::size_t num_vec_for_proj(std::size_t nproj) {
  switch (proj_mode_for(nproj)) {
    case ProjMode::Temp:
      return 1;
    case ProjMode::Pol:
      return 3;
    case ProjMode::VPol:
      return 4;
    default:
      return 0;
  }
}

namespace {

// The C tests `det->flag_init || pair->flag_init` and then dereferences
// both flag arrays, which is a null read when only one detector has flags.
// Each is guarded separately here.
inline bool flagged(const Det &det, std::size_t i) {
  return det.flag && det.flag[i];
}

// Map a sky pixel into a partial map. Returns false if the sample should be
// skipped; throws when error_missing is set.
inline bool remap(const Map &map, long &ipix, bool error_missing,
                  const char *what) {
  if (!map.partial()) return true;
  ipix = map.pixhash->repixelize(ipix);
  if (ipix >= 0) return true;
  if (error_missing) throw QpMapError(std::string(what) + ": pixel out of bounds");
  return false;
}

}  // namespace

void tod2map1(Pointing &mem, const Det &det, const PointData &pnt, Map &map) {
  const Options &opt = mem.opt();
  if (!opt.mean_aber && !pnt.ctime)
    throw QpPointError("tod2map1: ctime required if not mean_aber");

  const double w0 = det.weight;
  const double g = det.gain;
  const Mueller &m = det.mueller;

  const bool do_pol =
      at_least(map.vec_mode, VecMode::Pol) || at_least(map.proj_mode, ProjMode::Pol);
  const bool do_vpol =
      map.vec_mode == VecMode::VPol || map.proj_mode == ProjMode::VPol;

  double mt = m[0], mq = 0, mu = 0, mv = m[3];
  double wmt = w0 * m[0], wmq = 0, wmu = 0, wmv = w0 * m[3];

  for (std::size_t ii = 0; ii < pnt.n(); ++ii) {
    if (flagged(det, ii)) continue;
    const double ctime = pnt.time(ii);

    const Quat q = pnt.q_hwp
                       ? mem.bore2det_hwp(det.q_off, ctime, pnt.q_bore[ii],
                                          pnt.q_hwp[ii])
                       : mem.bore2det(det.q_off, ctime, pnt.q_bore[ii]);

    long ipix;
    double spp, cpp;
    mem.quat2pix(q, static_cast<int>(map.nside), PolOut::SinCos, ipix, spp, &cpp);

    if (!remap(map, ipix, opt.error_missing, "tod2map1")) continue;

    double w1 = w0;
    if (det.weights) {
      w1 = w0 * det.weights[ii];
      wmt = w1 * m[0];
      if (do_vpol) wmv = w1 * m[3];
    }

    if (do_pol) {
      mq = m[1] * cpp - m[2] * spp;
      mu = m[2] * cpp + m[1] * spp;
      if (!opt.polconv) mu *= -1;
      wmq = w1 * mq;
      wmu = w1 * mu;
    }

    if (det.tod && map.has_vec()) {
      const double gd = g * det.tod[ii];
      switch (map.vec_mode) {
        case VecMode::VPol:
          map.vrow(3)[ipix] += wmv * gd;
          [[fallthrough]];
        case VecMode::Pol:
          map.vrow(1)[ipix] += wmq * gd;
          map.vrow(2)[ipix] += wmu * gd;
          [[fallthrough]];
        case VecMode::Temp:
          map.vrow(0)[ipix] += wmt * gd;
          break;
        default:
          break;
      }
    }

    if (map.has_proj()) {
      switch (map.proj_mode) {
        case ProjMode::VPol:
          map.prow(0)[ipix] += wmt * mt;
          map.prow(1)[ipix] += wmt * mq;
          map.prow(2)[ipix] += wmt * mu;
          map.prow(3)[ipix] += wmt * mv;
          map.prow(4)[ipix] += wmq * mq;
          map.prow(5)[ipix] += wmq * mu;
          map.prow(6)[ipix] += wmq * mv;
          map.prow(7)[ipix] += wmu * mu;
          map.prow(8)[ipix] += wmu * mv;
          map.prow(9)[ipix] += wmv * mv;
          break;
        case ProjMode::Pol:
          map.prow(1)[ipix] += wmt * mq;
          map.prow(2)[ipix] += wmt * mu;
          map.prow(3)[ipix] += wmq * mq;
          map.prow(4)[ipix] += wmq * mu;
          map.prow(5)[ipix] += wmu * mu;
          [[fallthrough]];
        case ProjMode::Temp:
          map.prow(0)[ipix] += wmt * mt;
          break;
        default:
          break;
      }
    }
  }
}

void tod2map1_diff(Pointing &mem, const Det &det, const Det &pair,
                   const PointData &pnt, Map &map) {
  const Options &opt = mem.opt();
  if (!opt.mean_aber && !pnt.ctime)
    throw QpPointError("tod2map1_diff: ctime required if not mean_aber");

  const double w0 = det.weight, g = det.gain;
  const double w0_p = pair.weight, g_p = pair.gain;
  const Mueller &m = det.mueller;
  const Mueller &m_p = pair.mueller;

  double w = w0, w_p = w0_p;
  double wd = 0.5 * (w + w_p);
  const double mtd = 0.5 * (m[0] + m_p[0]);

  const bool do_pol =
      at_least(map.vec_mode, VecMode::Pol) || at_least(map.proj_mode, ProjMode::Pol);
  const bool do_vpol =
      map.vec_mode == VecMode::VPol || map.proj_mode == ProjMode::VPol;

  double alpha = 0, beta = 0, gamma = 0;
  double walpha = 0, wbeta = 0, wgamma = 0;

  for (std::size_t ii = 0; ii < pnt.n(); ++ii) {
    if (flagged(det, ii) || flagged(pair, ii)) continue;
    const double ctime = pnt.time(ii);

    Quat q, q_p;
    if (pnt.q_hwp) {
      q = mem.bore2det_hwp(det.q_off, ctime, pnt.q_bore[ii], pnt.q_hwp[ii]);
      q_p = mem.bore2det_hwp(pair.q_off, ctime, pnt.q_bore[ii], pnt.q_hwp[ii]);
    } else {
      q = mem.bore2det(det.q_off, ctime, pnt.q_bore[ii]);
      q_p = mem.bore2det(pair.q_off, ctime, pnt.q_bore[ii]);
    }

    long ipix, ipix_p;
    double spp, cpp, spp_p, cpp_p;
    const int nside = static_cast<int>(map.nside);
    mem.quat2pix(q, nside, PolOut::SinCos, ipix, spp, &cpp);
    mem.quat2pix(q_p, nside, PolOut::SinCos, ipix_p, spp_p, &cpp_p);

    if (!remap(map, ipix, opt.error_missing, "tod2map1_diff")) continue;
    if (!remap(map, ipix_p, opt.error_missing, "tod2map1_diff: pair")) continue;

    if (det.weights) w = w0 * det.weights[ii];
    if (pair.weights) w_p = w0_p * pair.weights[ii];
    if (det.weights || pair.weights) wd = 0.5 * (w + w_p);

    if (do_pol) {
      alpha = m[1] * cpp - m[2] * spp - (m_p[1] * cpp_p - m_p[2] * spp_p);
      beta = m[2] * cpp + m[1] * spp - (m_p[2] * cpp_p + m_p[1] * spp_p);
      if (!opt.polconv) beta *= -1;
      walpha = 0.5 * wd * alpha;
      wbeta = 0.5 * wd * beta;
    }

    if (do_vpol) {
      gamma = m[3] * cpp - m_p[3] * cpp_p;
      wgamma = 0.5 * wd * gamma;
    }

    if (det.tod && pair.tod && map.has_vec()) {
      const double delta = g * det.tod[ii] - g_p * pair.tod[ii];
      switch (map.vec_mode) {
        case VecMode::VPol:
          map.vrow(3)[ipix] += wgamma * delta;
          [[fallthrough]];
        case VecMode::Pol:
          map.vrow(1)[ipix] += walpha * delta;
          map.vrow(2)[ipix] += wbeta * delta;
          [[fallthrough]];
        case VecMode::Temp:
          map.vrow(0)[ipix] +=
              0.5 * wd * (g * m[0] * det.tod[ii] + g_p * m_p[0] * pair.tod[ii]);
          break;
        default:
          break;
      }
    }

    if (map.has_proj()) {
      switch (map.proj_mode) {
        case ProjMode::VPol:
          map.prow(0)[ipix] += wd * mtd;
          map.prow(1)[ipix] += 0.;
          map.prow(2)[ipix] += 0.;
          map.prow(3)[ipix] += 0.;
          map.prow(4)[ipix] += walpha * alpha;
          map.prow(5)[ipix] += walpha * beta;
          map.prow(6)[ipix] += walpha * gamma;
          map.prow(7)[ipix] += wbeta * beta;
          map.prow(8)[ipix] += wbeta * gamma;
          map.prow(9)[ipix] += wgamma * gamma;
          break;
        case ProjMode::Pol:
          map.prow(1)[ipix] += 0.;
          map.prow(2)[ipix] += 0.;
          map.prow(3)[ipix] += walpha * alpha;
          map.prow(4)[ipix] += walpha * beta;
          map.prow(5)[ipix] += wbeta * beta;
          [[fallthrough]];
        case ProjMode::Temp:
          map.prow(0)[ipix] += wd * mtd;
          break;
        default:
          break;
      }
    }
  }
}

void map2tod1(Pointing &mem, const Det &det, const PointData &pnt,
              const Map &map) {
  const Options &opt = mem.opt();
  if (!det.tod) throw QpInitError("map2tod1: det.tod not initialized");
  if (!opt.mean_aber && !pnt.ctime)
    throw QpPointError("map2tod1: ctime required if not mean_aber");

  const double g = det.gain;
  const Mueller &m = det.mueller;
  double mt = m[0], mq = 0, mu = 0, mv = m[3];

  const bool do_interp =
      opt.interp_pix &&
      (map.vec_mode == VecMode::Temp || map.vec_mode == VecMode::Pol ||
       map.vec_mode == VecMode::VPol);
  const bool do_deriv = at_least(map.vec_mode, VecMode::D1);
  const bool do_pol =
      at_least(map.vec_mode, VecMode::Pol) || at_least(map.proj_mode, ProjMode::Pol);

  // Plain read at one pixel, and the interpolated read across four.
  auto datum = [&](std::size_t row, long ipix) { return map.vrow(row)[ipix]; };
  auto idatum = [&](std::size_t row, const long pix[4], const double wt[4]) {
    return map.vrow(row)[pix[0]] * wt[0] + map.vrow(row)[pix[1]] * wt[1] +
           map.vrow(row)[pix[2]] * wt[2] + map.vrow(row)[pix[3]] * wt[3];
  };

  for (std::size_t ii = 0; ii < pnt.n(); ++ii) {
    if (flagged(det, ii)) continue;
    const double ctime = pnt.time(ii);

    const Quat q = pnt.q_hwp
                       ? mem.bore2det_hwp(det.q_off, ctime, pnt.q_bore[ii],
                                          pnt.q_hwp[ii])
                       : mem.bore2det(det.q_off, ctime, pnt.q_bore[ii]);

    long ipix;
    double spp, cpp, dtheta = 0, dphi = 0;
    long pix[4] = {0, 0, 0, 0};
    double weight[4] = {0, 0, 0, 0};
    const int nside = static_cast<int>(map.nside);

    if (do_deriv || do_interp) {
      double ra, dec;
      mem.quat2radec(q, DecOut::Dec, PolOut::SinCos, ra, dec, spp, &cpp);
      ipix = mem.radec2pix(ra, dec, nside);
      mem.pixel_offset(nside, ipix, ra, dec, dtheta, dphi);
      if (do_interp)
        map.pixinfo->interpol(ra, dec, opt.pix_order == 1, pix, weight);
    } else {
      mem.quat2pix(q, nside, PolOut::SinCos, ipix, spp, &cpp);
    }

    if (map.partial()) {
      ipix = map.pixhash->repixelize(ipix);
      if (ipix < 0) {
        if (opt.error_missing)
          throw QpMapError("map2tod1: pixel out of bounds");
        if (opt.nan_missing) det.tod[ii] = 0.0 / 0.0;
        continue;
      }
      if (do_interp) {
        int bad_pix = 0;
        for (int jj = 0; jj < 4; ++jj) {
          pix[jj] = map.pixhash->repixelize(pix[jj]);
          if (pix[jj] >= 0) continue;
          if (opt.error_missing)
            throw QpMapError("map2tod1: neighbor pixel out of bounds");
          if (opt.interp_missing) {
            // drop this neighbour and renormalize the rest
            double norm1 = 0.0, norm2 = 0.0;
            for (int kk = 0; kk < 4; ++kk) {
              norm1 += weight[kk];
              if (kk != jj) norm2 += weight[kk];
            }
            pix[jj] = 0;
            weight[jj] = 0;
            for (int kk = 0; kk < 4; ++kk)
              if (kk != jj) weight[kk] *= norm1 / norm2;
            bad_pix += 1;
          } else {
            bad_pix = 1;
            break;
          }
        }
        // still usable as long as one neighbour survived
        if (opt.interp_missing && bad_pix < 4) bad_pix = 0;
        if (bad_pix) {
          if (opt.nan_missing) det.tod[ii] = 0.0 / 0.0;
          continue;
        }
      }
    }

    if (do_pol) {
      mq = m[1] * cpp - m[2] * spp;
      mu = m[2] * cpp + m[1] * spp;
      if (!opt.polconv) mu *= -1;
    }

    // do_interp is loop-invariant, and false for every derivative mode, so
    // one reader serves both paths.
    auto rd = [&](std::size_t r) {
      return do_interp ? idatum(r, pix, weight) : datum(r, ipix);
    };
    auto T = [&](std::size_t r) { return mt * rd(r); };
    auto P = [&](std::size_t r) {
      return mt * rd(r) + mq * rd(r + 1) + mu * rd(r + 2);
    };

    switch (map.vec_mode) {
      case VecMode::VPol:
        det.tod[ii] += g * (P(0) + mv * rd(3));
        break;
      case VecMode::D2Pol:
        det.tod[ii] += g * (dphi * dphi * P(15) + dtheta * dphi * P(12) +
                            dtheta * dtheta * P(9));
        [[fallthrough]];
      case VecMode::D1Pol:
        det.tod[ii] += g * (dphi * P(6) + dtheta * P(3));
        [[fallthrough]];
      case VecMode::Pol:
        det.tod[ii] += g * P(0);
        break;
      case VecMode::D2:
        det.tod[ii] += g * (dphi * dphi * T(5) + dtheta * dphi * T(4) +
                            dtheta * dtheta * T(3));
        [[fallthrough]];
      case VecMode::D1:
        det.tod[ii] += g * (dphi * T(2) + dtheta * T(1));
        [[fallthrough]];
      case VecMode::Temp:
        det.tod[ii] += g * T(0);
        break;
      default:
        break;
    }
  }
}

void add_map(Map &map, const Map &local) {
  if (map.vec_mode != local.vec_mode || map.proj_mode != local.proj_mode ||
      map.nside != local.nside || map.npix != local.npix)
    throw QpMapError("add_map: maps are not compatible");

  // The zero test is not an optimization that can be dropped: (-0.0) += 0.0
  // yields +0.0, so skipping zeros preserves signed zeros exactly as the C
  // reduction does.
  auto accum = [npix = map.npix](double *dst, const double *src) {
    for (std::size_t p = 0; p < npix; ++p)
      if (src[p] != 0) dst[p] += src[p];
  };

  if (map.has_vec() && local.has_vec())
    for (std::size_t i = 0; i < num_vec(map.vec_mode); ++i)
      accum(map.vrow(i), local.vrow(i));

  if (map.has_proj() && local.has_proj())
    for (std::size_t i = 0; i < num_proj(map.proj_mode); ++i)
      accum(map.prow(i), local.prow(i));
}

}  // namespace qp
