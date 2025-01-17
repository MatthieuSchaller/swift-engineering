/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Filip Husko (filip.husko@durham.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#ifndef SWIFT_SPIN_JET_BLACK_HOLES_SPIN_H
#define SWIFT_SPIN_JET_BLACK_HOLES_SPIN_H

/* Standard headers */
#include <float.h>

/* Local includes */
#include "black_holes_properties.h"
#include "black_holes_struct.h"
#include "inline.h"
#include "physical_constants.h"

/**
 * @brief Compute the radius of the horizon of a BH particle in gravitational
 * units.
 *
 * @param a Black hole spin, -1 < a < 1.
 */
__attribute__((always_inline)) INLINE static float r_hor(float a) {
  return 1. + sqrtf((1. - a) * (1. + a));
}

/**
 * @brief Compute the radius of the innermost stable circular orbit of a
 * BH particle in gravitational units.
 *
 * The expression is given in Appendix B of Fiacconi et al. (2018) or eqn. 4 in
 * Griffin et al. (2019).
 *
 * @param a Black hole spin, -1 < a < 1.
 */
__attribute__((always_inline)) INLINE static float r_isco(float a) {
  const float Z1 = 1. + (cbrtf((1. + fabsf(a)) * (1. - a * a)) +
                         cbrtf((1. - fabsf(a)) * (1. - a * a)));
  const float Z2 = sqrtf(3. * a * a + Z1 * Z1);

  const float R_ISCO =
      3. + Z2 - a / fabsf(a) * sqrtf((3. - Z1) * (3. + Z1 + 2. * Z2));

#ifdef SWIFT_DEBUG_CHECKS
  if (Z1 > 3.) {
    error(
        "Something went wrong with calculation of Z1 factor for r_isco of"
        " black holes. Z1 is %f instead of Z1 > 3.",
        Z1);
  }

  if ((3. + Z1 + 2. * Z2) < 0.) {
    error(
        "Something went wrong with calculation of (3. + Z1 + 2. * Z2 ) "
        "factor for r_isco of black holes. (3. + Z1 + 2. * Z2 ) is %f instead "
        "of"
        " (3. + Z1 + 2. * Z2 ) > 0.",
        3. + Z1 + 2. * Z2);
  }

  if (R_ISCO < 1.) {
    error(
        "Something went wrong with calculation of R_ISCO of black holes. "
        "R_ISCO is %f instead >= 1.",
        R_ISCO);
  }
#endif

  return R_ISCO;
}

/**
 * @brief Compute the magnitude of the angular momentum of the black hole
 * given its spin.
 *
 * @param a Black hole spin magnitude, 0 < a < 1.
 * @param constants Physical constants (in internal units).
 */
__attribute__((always_inline)) INLINE static float j_BH(
    struct bpart* bp, const struct phys_const* constants) {

  const float J_BH =
      fabs(bp->subgrid_mass * bp->subgrid_mass * bp->spin *
           constants->const_newton_G / constants->const_speed_light_c);

#ifdef SWIFT_DEBUG_CHECKS
  if (J_BH <= 0.) {
    error(
        "Something went wrong with calculation of j_BH of black holes. "
        " J_BH is %f instead of J_BH > 0.",
        J_BH);
  }
#endif

  return J_BH;
}

/**
 * @brief Compute the gravitational radius of a black hole.
 *
 * @param a Black hole mass.
 * @param constants Physical constants (in internal units).
 */
__attribute__((always_inline)) INLINE static float R_gravitational(
    float mass, const struct phys_const* constants) {

  const float r_G =
      mass * constants->const_newton_G /
      (constants->const_speed_light_c * constants->const_speed_light_c);

#ifdef SWIFT_DEBUG_CHECKS
  if (r_G <= 0.) {
    error(
        "Something went wrong with calculation of R_G of black holes. "
        " R_G is %f instead of R_G > 0.",
        r_G);
  }
#endif

  return r_G;
}

/**
 * @brief Compute the warp radius of a black hole particle.
 *
 * The result depends on bp->accretion_mode (thick disk, thin disk or
 * slim disk). For the thick disk and slim disk, the radius is calculated
 * from Lubow et al. (2002), eqn. 22 with x=1. The result will be different
 * only due to different aspect ratios H/R=h_0.
 *
 * For the thin disk, the result depends on props->TD_region (B - region b from
 * Shakura & Sunyaev 1973, C - region c from Shakura & Sunyaev 1973). The warp
 * radii are taken as eqns. 11 from Griffin et al. (2019) and A8 from Fiacconi
 * et al. (2018), respectively.
 *
 * For the thin disk we also have to include the possibility that the self-
 * gravity radius is smaller than the warp radius. In this case r_warp=r_sg
 * because the disk cannot be larger than the self-gravity radius, and the
 * entire disk is warped. The sg radius is taken as eqns. 16 in Griffin et al.
 * (2019) and A6 in Fiacconi et al. (2018), respectively.
 *
 * @param bp Pointer to the b-particle data.
 * @param constants Physical constants (in internal units).
 * @param props Properties of the black hole scheme.
 */
__attribute__((always_inline)) INLINE static float r_warp(
    struct bpart* bp, const struct phys_const* constants,
    const struct black_holes_props* props) {

  /* Define placeholder variable for the result */
  float Rw = -1.;

  /* Gravitational radius */
  const float R_G = R_gravitational(bp->subgrid_mass, constants);

  /* Start branching depending on which accretion mode the BH is in */
  if (bp->accretion_mode == BH_thick_disc) {

    /* Eqn. 22 from Lubow et al. (2002) with H/R=h_0_ADAF (thick disk) */
    const float base = 15.36 * fabsf(bp->spin) / props->h_0_ADAF_2;
    Rw = R_G * powf(base, 0.4);
  } else if (bp->accretion_mode == BH_slim_disc) {

    /* Eqn. 22 from Lubow et al. (2002) with H/R=1/gamma_SD (slim disk) */
    const float base = 15.36 * fabsf(bp->spin) * props->gamma_SD;
    Rw = R_G * powf(base, 0.4);
  } else if (bp->accretion_mode == BH_thin_disc) {

    /* Start branching depending on which region of the thin disk we wish to
       base the model upon (TD_region=B: region b from Shakura & Sunyaev 1973,
       or TD_region=C: region c) */
    if (props->TD_region == TD_region_B) {

      /* Calculate different factors in eqn. 11 (Griffin et al. 2019) for warp
          radius of region b in Shakura & Sunyaev (1973) */
      float mass_factor =
          powf(bp->subgrid_mass / (1e8 * constants->const_solar_mass), 0.2);
      float edd_factor = powf(bp->eddington_fraction, 0.4);

      /* Gather the factors and finalize calculation */
      const float base = mass_factor * fabsf(bp->spin) /
                         (props->xi_TD * props->alpha_factor_08 * edd_factor);
      const float rw = 3410. * 2. * R_G * powf(base, 0.625);

      /* Self-gravity radius in region b: eqn. 16 in Griffin et al. */
      mass_factor =
          powf(bp->subgrid_mass / (1e8 * constants->const_solar_mass), -0.961);
      edd_factor = powf(bp->eddington_fraction, -0.353);

      const float rs = 4790. * 2. * R_G * mass_factor *
                       props->alpha_factor_0549 * edd_factor;

      /* Take the minimum */
      Rw = min(rs, rw);
    }

    if (props->TD_region == TD_region_C) {

      /* Calculate different factors in eqn. A8 (Fiacconi et al. 2018) */
      float mass_factor =
          powf(bp->subgrid_mass / (1e6 * constants->const_solar_mass), 0.2);
      float edd_factor = powf(bp->eddington_fraction, 0.3);

      /* Gather the factors and finalize calculation */
      const float base = mass_factor * fabsf(bp->spin) /
                         (props->xi_TD * props->alpha_factor_02 * edd_factor);
      const float rw = 1553. * 2. * R_G * powf(base, 0.5714);

      /* Repeat the same for self-gravity radius - eqn. A6 in F2018 */
      mass_factor =
          powf(bp->subgrid_mass / (1e6 * constants->const_solar_mass), -1.1556);
      edd_factor = powf(bp->eddington_fraction, -0.48889);

      const float rs = 1.2 * 100000. * 2. * R_G * mass_factor *
                       props->alpha_factor_06222 * edd_factor;

      /* Take the minimum */
      Rw = min(rs, rw);
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (Rw < 0.) {
    error(
        "Something went wrong with calculation of Rw of black holes. "
        " Rw is %f instead of Rw >= 0.",
        Rw);
  }
#endif

  return Rw;
}

/**
 * @brief Compute the warp mass of a black hole particle.
 *
 * Calculated as the integral of the surface density of the disk up to R_warp.
 * The result again depends on type of accretion mode, both due to different
 * R_warp and different surface densities.
 *
 * The surface densities for the thick and slim disk take the same form
 * (eqn. 2.3 in Narayan & Yi 1995 for rho, and then sigma = rho * 2H =
 * dot(M_BH) / (2pi * R * abs(v_r))). They differ due to different radial
 * radial velocities in the disks: v_r = -alpha * v_0 * v_K (with v_K
 * the Keplerian velocity). These differences are encoded in the numerical
 * constant v_0, which depends on alpha in Narayan & Yi for the thick disk,
 * and is roughly constant for the slim disk (Wang & Zhou 1999).
 *
 * For the thin disk the surface densities are more complex, and again depend
 * on which region of the disk is chosen to be modelled (region b or c from
 * Shakura & Sunyaev 1973). Sigma for region b is given by eqn. 7 in Griffin
 * et al. (2019) and for region c, it is not given explicitly but can be
 * calculated based on Appendix A in Fiacconi et al. (2018).
 *
 * @param bp Pointer to the b-particle data.
 * @param constants Physical constants (in internal units).
 * @param props Properties of the black hole scheme.
 */
__attribute__((always_inline)) INLINE static float m_warp(
    struct bpart* bp, const struct phys_const* constants,
    const struct black_holes_props* props) {

  /* Define placeholder variable for the result */
  float Mw = -1.;

  /* Gravitational radius */
  const float R_G = R_gravitational(bp->subgrid_mass, constants);

  /* Start branching depending on which accretion mode the BH is in */
  if ((bp->accretion_mode == BH_thick_disc) ||
      (bp->accretion_mode == BH_slim_disc)) {

    /* Define v_0, the only factor which differs between thick and slim
       disc */
    float v_0;
    if (bp->accretion_mode == BH_thick_disc) {
      v_0 = props->v_0_ADAF;
    } else {
      v_0 = props->gamma_SD_inv;
    }

    /* Final result based on eqn. 2.3 in Narayan & Yi 1995*/
    Mw = 2. * bp->accretion_rate /
         (3. * props->alpha_acc * v_0 *
          sqrtf(bp->subgrid_mass * constants->const_newton_G)) *
         powf(r_warp(bp, constants, props), 1.5);
  } else {

    /* Start branching depending on which region of the thin disk we wish to
       base the model upon (TD_region=B: region b from Shakura & Sunyaev 1973,
       or TD_region=C: region c) */
    if (props->TD_region == TD_region_B) {

      /* Calculate different factors that appear in result for M_warp */
      const float mass_factor =
          powf(bp->subgrid_mass / (1e8 * constants->const_solar_mass), 2.2);
      const float edd_factor = powf(bp->eddington_fraction, 0.6);
      const float R_factor =
          powf(r_warp(bp, constants, props) / (2. * R_G), 1.4);

      /* Gather factors and finalize calculation */
      Mw = constants->const_solar_mass * 1.35 * mass_factor *
           props->alpha_factor_08_inv * edd_factor * R_factor;
    }
    if (props->TD_region == TD_region_C) {

      /* Same as above but for region c of disk */
      const float mass_factor =
          powf(bp->subgrid_mass / (1e6 * constants->const_solar_mass), 2.2);
      const float edd_factor = powf(bp->eddington_fraction, 0.7);
      const float R_factor =
          powf(r_warp(bp, constants, props) / (2. * R_G), 1.25);

      Mw = constants->const_solar_mass * 0.01 * mass_factor *
           props->alpha_factor_08_inv_10 * edd_factor * R_factor;
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (Mw < 0.) {
    error(
        "Something went wrong with calculation of Mw of black holes. "
        " Mw is %f instead of Mw >= 0.",
        Mw);
  }
#endif

  return Mw;
}

/**
 * @brief Compute the warp angular momentum of a black hole particle.
 *
 * Calculated as the integral of the surface density times the specific
 * angular momentum of the disk up to R_warp. The result depends on type
 * of accretion mode, due to different R_warp, surface densities and spec.
 * ang. momenta of the disks.
 *
 * The surface densities are the same as for M_warp. For the thin disk, the
 * spec. ang. mom. is L(R) = R * v_K(R), because orbits are perfectly circular.
 * For the thick and slim disk, this is replaced by L(R) = Omega_0 * R * v_K(R)
 * , with Omega_0 a numerical constant between 0 and 1 which encodes the fact
 * that rotation is slower in the two disks. The values for Omega_0 are given
 * in Narayan & Yi (1995) and Wang & Zhou (1999) for the thick and slim disk,
 * respectively.
 *
 * @param bp Pointer to the b-particle data.
 * @param constants Physical constants (in internal units).
 * @param props Properties of the black hole scheme.
 */
__attribute__((always_inline)) INLINE static float j_warp(
    struct bpart* bp, const struct phys_const* constants,
    const struct black_holes_props* props) {

  /* Define placeholder variable for the result */
  float Jw = -1.;

  /* Start branching depending on which accretion mode the BH is in */
  if ((bp->accretion_mode == BH_thick_disc) ||
      (bp->accretion_mode == BH_slim_disc)) {

    /* Get numerical constants for radial and tangential velocities for the
       thick and slim disk, which factor into the surface density and spec.
       ang. mom., respectively */
    float v_0 = 0.;
    float omega_0 = 0.;
    if (bp->accretion_mode == BH_thick_disc) {
      v_0 = props->v_0_ADAF;
      omega_0 = props->omega_0_ADAF;
    } else {
      v_0 = props->gamma_SD_inv;
      omega_0 = props->gamma_SD_inv;
    }

    /* Gather factors for the final result  */
    Jw = 2. * bp->accretion_rate * omega_0 / (2. * props->alpha_acc * v_0) *
         r_warp(bp, constants, props) * r_warp(bp, constants, props);
  } else {

    /* Start branching depending on which region of the thin disk we wish to
       base the model upon (TD_region=B: region b from Shakura & Sunyaev 1973,
       or TD_region=C: region c). The warp radius can generally be related to
       the warp mass and radius, if one assumes Keplerian rotation, with the
       following relation: J_warp = (c+2)/(c+5/2) * M_warp * sqrt(M_BH * G *
       R_warp), where c is the slope of the surface density profile: sigma~R^c.
       For region b, c=-3/5 (see Griffin et al. 2019), and for region c, c=-3/4
       (see Fiacconi et al. 2018). */
    if (props->TD_region == TD_region_B) {
      Jw = 0.737 * m_warp(bp, constants, props) *
           sqrtf(bp->subgrid_mass * constants->const_newton_G *
                 r_warp(bp, constants, props));
    }
    if (props->TD_region == TD_region_C) {
      Jw = 0.714 * m_warp(bp, constants, props) *
           sqrtf(bp->subgrid_mass * constants->const_newton_G *
                 r_warp(bp, constants, props));
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (Jw < 0.) {
    error(
        "Something went wrong with calculation of Jw of black holes. "
        " Jw is %f instead of Jw >= 0.",
        Jw);
  }
#endif

  return Jw;
}

/**
 * @brief Compute the spin-dependant radiative efficiency of a BH particle in
 * the radiatively efficient (thin disc) regime.
 *
 * This is eqn. 3 in Griffin et al. (2019), based on Novikov & Thorne (1973).
 *
 * @param a Black hole spin, -1 < a < 1.
 */
__attribute__((always_inline)) INLINE static float eps_NT(float a) {

#ifdef SWIFT_DEBUG_CHECKS
  if (r_isco(a) <= 0.6667) {
    error(
        "Something went wrong with calculation of eps_NT of black holes. "
        " r_isco is %f instead of r_isco > 1.",
        r_isco(a));
  }
#endif

  return 1. - sqrtf(1. - 2. / 3. / r_isco(a));
}

/**
 * @brief Compute the spin- and accretion rate-dependant radiative efficiency
 * of a BH particle in the super-Eddington (slim disk) regime.
 *
 * This is eqn. 3 in Madau et al. (2014), which is based on numerical GR
 * results by Sadowski (2009).
 *
 * @param a Black hole spin, -1 < a < 1.
 * @param m_dot Accretion rate normalized to the Eddington rate.
 */
__attribute__((always_inline)) INLINE static float eps_SD(float a, float mdot) {
  const float B = powf(4.627 - 4.445 * a, -0.5524);
  const float C = powf(827.3 - 718.1 * a, -0.706);
  const float A = powf(0.9663 - 0.9292 * a, -0.5693);

#ifdef SWIFT_DEBUG_CHECKS
  if (mdot <= 0.) {
    error(
        "The calculation of eps_SD was called even though mdot is %f. "
        " This function should not have been called if the accretion rate is "
        " not > 0.",
        mdot);
  }
#endif

  return 0.1 / mdot * (0.985 / (B + 1.6 / mdot) + 0.015 / (C + 1.6 / mdot)) * A;
}

/**
 * @brief Decide which regime (mode) of accretion the BH particle is in.
 *
 * The possible modes are the thick disk, thin disk and slim disk, in
 * order of increasing accretion rate. The transition from thick to thin disk
 * is at 0.4*alpha^2, based on current theory (Yuan & Narayan 2014). The
 * transition from thin to slim disk occurs when the slim disk efficiency
 * becomes sufficiently weak compared to the thin disk one. We parametrize
 * this transition as occuring at eps_SD = props->TD_SD_eps_r_threshold *
 * eps_TD, with props->TD_SD_eps_r_threshold < 1 of order 0.5
 *
 * @param bp Pointer to the b-particle data.
 * @param constants Physical constants (in internal units).
 * @param props Properties of the black hole scheme.
 */
__attribute__((always_inline)) INLINE static void decide_mode(
    struct bpart* bp, const struct black_holes_props* props) {
  if (bp->eddington_fraction < props->mdot_crit_ADAF) {
    bp->accretion_mode = BH_thick_disc;
  } else {
    if ((eps_SD(bp->spin, bp->eddington_fraction) <
         props->TD_SD_eps_r_threshold * eps_NT(bp->spin)) &&
        (props->include_slim_disk)) {
      bp->accretion_mode = BH_slim_disc;
    } else {
      bp->accretion_mode = BH_thin_disc;
    }
  }

  /* If we do not include radiative feedback, then we force the disk to be in
     the thick disk mode */
  if (props->turn_off_radiative_feedback) {
    bp->accretion_mode = BH_thick_disc;
  }

  /* similar for if we do not include jets - we force the disk to be thin */
  if (props->include_jets == 0) {
    bp->accretion_mode = BH_thin_disc;
  }
}

/**
 * @brief Compute the aspect ratio of the subgrid accretion disk.
 *
 * The result depends on bp->accretion_mode (thick disk, thin disk or
 * slim disk). For the thick disk and slim disk, the aspect ratio is
 * a constant, H/R = h_0.
 *
 * For the thin disk, the result depends on props->TD_region (B - region b from
 * Shakura & Sunyaev 1973, C - region c from Shakura & Sunyaev 1973). In region
 * b, we take H/R as eqn. 8 in Griffin et al. (2019), and in region c H/r is
 * taken directly as eqn. 2.19 from Shakura & Sunyaev (1973).
 *
 * @param bp Pointer to the b-particle data.
 * @param constants Physical constants (in internal units).
 * @param props Properties of the black hole scheme.
 */
__attribute__((always_inline)) INLINE static float aspect_ratio(
    struct bpart* bp, const struct phys_const* constants,
    const struct black_holes_props* props) {

  /* Define placeholder variable for the result */
  float h_0 = -1.;

  /* Start branching depending on which accretion mode the BH is in */
  if ((bp->accretion_mode == BH_thick_disc) ||
      (bp->accretion_mode == BH_slim_disc)) {
    if (bp->accretion_mode == BH_thick_disc) {
      h_0 = props->h_0_ADAF;
    } else {
      h_0 = 0.5 * props->gamma_SD_inv;
    }
  } else {

    /* Start branching depending on which region of the thin disk we wish to
       base the model upon (TD_region=B: region b from Shakura & Sunyaev 1973,
       or TD_region=C: region c). */
    if (props->TD_region == TD_region_B) {

      /* Compute factors for eqn. 8 in Griffin et al. (2019). */
      const float mass_factor =
          powf(bp->subgrid_mass / (1e8 * constants->const_solar_mass), -0.1);
      const float edd_factor = powf(bp->eddington_fraction, 0.2);
      const float R_G = R_gravitational(bp->subgrid_mass, constants);
      const float R_factor =
          powf(r_warp(bp, constants, props) / (2. * R_G), 0.05);

      /* Gather factors and finalize calculation. */
      h_0 = 1.25 * 0.001 * mass_factor * props->alpha_factor_01 * edd_factor *
            R_factor;
    }
    if (props->TD_region == TD_region_C) {

      /* Compute factors for eqn. 2.19 in Shakura & Sunyaev (1973). */
      const float mass_factor =
          powf(bp->subgrid_mass / (1e8 * constants->const_solar_mass), -0.1);
      const float edd_factor = powf(bp->eddington_fraction, 0.15);
      const float R_G = R_gravitational(bp->subgrid_mass, constants);
      const float R_factor =
          powf(r_warp(bp, constants, props) / (2. * R_G), 0.125);

      /* Gather factors and finalize calculation. */
      h_0 = 1.15 * 0.001 * mass_factor * props->alpha_factor_01 * edd_factor *
            R_factor;
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (h_0 <= 0.) {
    error(
        "Something went wrong with calculation of h_0 of black holes. "
        " h_0 is %f instead of h_0 > 0.",
        h_0);
  }
#endif

  return h_0;
}

/**
 * @brief Compute the jet efficiency of a BH particle.
 *
 * The result depends on bp->accretion_mode (thick disk, thin disk or
 * slim disk), through the varying H/R aspect ratios.
 *
 * The equation implemented is eqn. 9 from Tchekhovskoy et al. (2010), with the
 * dimensionless magnetic flux phi taken as eqn. 9 from Narayan et al. (2021).
 *
 * The dependence on the aspect ratio comes from results in Tchekhovskoy et al.
 * (2014) and the dependence in classical Blandford & Znajek (1979) jet theory.
 *
 * @param bp Pointer to the b-particle data.
 * @param constants Physical constants (in internal units).
 * @param props Properties of the black hole scheme.
 */
__attribute__((always_inline)) INLINE static float jet_efficiency(
    struct bpart* bp, const struct black_holes_props* props) {

  float jet_eff = -1.;
  if (props->fix_jet_efficiency) {
    jet_eff = props->jet_efficiency;
  } else {
    const float kappa = 0.05;
    const float horizon_ang_vel =
        bp->spin / (2. * (1. + sqrtf(1. - bp->spin * bp->spin)));
    const float phi = -20.2 * bp->spin * bp->spin * bp->spin -
                      14.9 * bp->spin * bp->spin + 34. * bp->spin + 52.6;
    jet_eff = kappa * 0.25 * M_1_PI * phi * phi *
              powf(bp->aspect_ratio * 3.333, props->jet_h_r_slope) *
              horizon_ang_vel * horizon_ang_vel *
              (1. + 1.38 * horizon_ang_vel * horizon_ang_vel -
               9.2 * horizon_ang_vel * horizon_ang_vel * horizon_ang_vel *
                   horizon_ang_vel);
  }

  /* Turn off jet feedback if we want to do that */
  if (props->include_jets == 0) {
    jet_eff = 0.;
  }

  /* Turn off jets in thin disk mode if we want to do that */
  if ((props->turn_off_secondary_feedback) &&
      (bp->accretion_mode == BH_thin_disc)) {
    jet_eff = 0.;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (jet_eff < 0.) {
    error(
        "Something went wrong with calculation of jet efficiency of black "
        "holes. jet_eff is %f instead of jet_eff >= 0.",
        jet_eff);
  }
#endif

  return jet_eff;
}

/**
 * @brief Compute the radiative efficiency of a BH particle.
 *
 * The result depends on bp->accretion_mode (thick disk, thin disk or
 * slim disk), since all modes have different radiative physics.
 *
 * For the thin disk, we assume the Novikov-Thorne (1973) radiative efficiency
 * based on general relativity. For the slim disk, we take the fit from Madau
 * et al. (2014), which is based on numerical GR results by Sadowski (2009).
 * For the thick disk, we assume radiative efficiencies from Mahadevan et al.
 * (1997).
 *
 * @param bp Pointer to the b-particle data.
 * @param constants Physical constants (in internal units).
 * @param props Properties of the black hole scheme.
 */
__attribute__((always_inline)) INLINE static float rad_efficiency(
    struct bpart* bp, const struct black_holes_props* props) {

  /* Calculate Novikov-Thorne efficiency, which will be needed twice. */
  const float eps_TD = eps_NT(bp->spin);

  /* Define placeholder variable for the result */
  float rad_eff = -1;

  if (props->fix_radiative_efficiency) {
    rad_eff = props->radiative_efficiency;
  } else {

    /* Start branching depending on which accretion mode the BH is in */
    if (bp->accretion_mode == BH_thin_disc) {

      /* Assign Novikov-Thorne efficiency to the thin disk. */
      rad_eff = eps_TD;
    } else if (bp->accretion_mode == BH_slim_disc) {

      /* Assign Madau 2014 efficiency to the slim disk. */
      rad_eff = eps_SD(bp->spin, bp->eddington_fraction);
    } else {

#ifdef SWIFT_DEBUG_CHECKS
      if (props->beta_acc > 1.) {
        error(
            "Something went wrong with calculation of radiative efficiency of "
            " black holes. beta_acc is %f instead of beta_acc < 1.",
            props->beta_acc);
      }
#endif

      /* Assign Mahadevan 1997 efficiency to the thick disk. */
      if (bp->eddington_fraction < props->edd_crit_thick) {
        rad_eff = 4.8 * eps_TD / r_isco(bp->spin) * (1. - props->beta_acc) *
                  props->delta_ADAF;
      } else {
        rad_eff = 2.4 * eps_TD / r_isco(bp->spin) * props->beta_acc *
                  bp->eddington_fraction * props->alpha_acc_2_inv;
      }
    }
  }

  /* Turn off radiative feedback if we want to do that */
  if (props->turn_off_radiative_feedback) {
    rad_eff = 0.;
  }

  /* Turn off radiation in the thick disk mode if we want to do that */
  if ((props->turn_off_secondary_feedback) &&
      (bp->accretion_mode == BH_thick_disc)) {
    rad_eff = 0.;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (rad_eff < 0.) {
    error(
        "Something went wrong with calculation of radiative efficiency of "
        " black holes. rad_eff is %f instead of rad_eff >= 0.",
        rad_eff);
  }
#endif

  return rad_eff;
}

/**
 * @brief Compute the spec. ang. mom. at the inner radius of a BH particle.
 *
 * The result depends on bp->accretion_mode (thick disk, thin disk or
 * slim disk), since advection-dominated modes (thick and slim disk)
 * have more radial orbits.
 *
 * For the thin disk, we assume that the spec. ang. mom. consumed matches that
 * of the innermost stable circular orbit (ISCO). For the other two modes, we
 * assume that the accreted ang. mom. at the event horizon is 45 per cent of
 * that at the ISCO, based on the fit from Benson & Babul (2009).
 *
 * @param bp Pointer to the b-particle data.
 * @param constants Physical constants (in internal units).
 * @param props Properties of the black hole scheme.
 */
__attribute__((always_inline)) INLINE static float l_acc(
    struct bpart* bp, const struct phys_const* constants,
    const struct black_holes_props* props) {

  /* Define placeholder variable for the result */
  float L = -1.;

#ifdef SWIFT_DEBUG_CHECKS
  if (r_isco(bp->spin) <= 0.6667) {
    error(
        "Something went wrong with calculation of l_acc of black holes. "
        " r_isco is %f instead of r_isco > 1.",
        r_isco(bp->spin));
  }
#endif

  /* Spec. ang. mom. at ISCO */
  const float L_ISCO = 0.385 * (1. + 2. * sqrtf(3. * r_isco(bp->spin) - 2.));

  /* Branch depending on which accretion mode the BH is in */
  if ((bp->accretion_mode == BH_thick_disc) ||
      (bp->accretion_mode == BH_slim_disc)) {
    L = 0.45 * L_ISCO;
  } else {
    L = L_ISCO;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (L <= 0.) {
    error(
        "Something went wrong with calculation of l_acc of black holes. "
        " l_acc is %f instead of l_acc > 0.",
        L);
  }
#endif

  return L;
}

/**
 * @brief Compute the evolution of the spin of a BH particle.
 *
 * The result depends on bp->accretion_mode (thick disk, thin disk or
 * slim disk), due to differing spec. ang. momenta as well as jet and
 * radiative efficiencies.
 *
 * This equation corresponds to eqn. 2 in Benson & Babul (2009), including
 * a jet spindown term.
 *
 * @param bp Pointer to the b-particle data.
 * @param constants Physical constants (in internal units).
 * @param props Properties of the black hole scheme.
 */
__attribute__((always_inline)) INLINE static float da_dln_mbh_0(
    struct bpart* bp, const struct phys_const* constants,
    const struct black_holes_props* props) {
  const float a = bp->spin;

  if ((a == 0.) || (a < -0.9981) || (a > 0.9981)) {
    error(
        "The da_dln_mbh_0 function was called and spin is %f. Spin should "
        " not be a = 0, a < -0.998 or a > 0.998.",
        a);
  }

  float spinup_rate = 0.;

  if (props->include_GRMHD_spindown) {
    if (bp->accretion_mode == BH_thin_disc) {
      spinup_rate = l_acc(bp, constants, props) -
                    2. * a * (1. - rad_efficiency(bp, props));
    } else {
      spinup_rate = 0.45 - 12.53 * a - 7.8 * a * a + 9.44 * a * a * a +
                    5.71 * a * a * a * a - 4.03 * a * a * a * a * a;
    }
  } else {
    spinup_rate =
        l_acc(bp, constants, props) -
        2. * a * (1. - rad_efficiency(bp, props)) -
        sqrtf(1. - a * a) / a *
            (a * a + (1. + sqrtf(1. - a * a)) * (1. + sqrtf(1. - a * a))) *
            jet_efficiency(bp, props);
  }

  return spinup_rate;
}

/**
 * @brief Compute the jet kick velocity to be used for jet feedback.
 *
 * @param bp The #bpart doing feedback.
 * @param props Properties of the BH scheme.
 * @param cosmo The current cosmological model.
 * @param constants The physical constants (in internal units).
 */
__attribute__((always_inline)) INLINE static float black_hole_feedback_dv_jet(
    const struct bpart* bp, const struct black_holes_props* props,
    const struct cosmology* cosmo, const struct phys_const* constants) {

  float v_jet = -1.;
  if (props->AGN_jet_velocity_model == AGN_jet_velocity_BH_mass) {

    /* Assign the halo mass according to an empirical relation given in the
       parameter file */
    const float halo_mass =
        powf(bp->subgrid_mass / props->v_jet_BH_mass_scaling_reference_mass,
             props->v_jet_BH_mass_scaling_slope);

    /* Get the critical density and virial overdensity at this redshift */
    const float critical_density = cosmo->critical_density;
    const float overdensity = cosmo->overdensity_BN98;

    /* Gather the previous factors and compute the virial radius, virial
       velocity and finally the sound speed in the hot gas */
    const float virial_radius =
        cbrtf(3. * halo_mass / (4. * M_PI * overdensity * critical_density));
    const float virial_velocity =
        sqrtf(bp->group_mass * constants->const_newton_G / virial_radius);
    const float sound_speed = sqrtf(5. / 3. * 0.5) * virial_velocity;

    /* Return the jet velocity as some factor times the sound speed */
    v_jet = fmaxf(props->v_jet_min, props->v_jet_cs_ratio * sound_speed);

  } else if (props->AGN_jet_velocity_model == AGN_jet_velocity_constant) {
    v_jet = props->v_jet;
  } else {
    error(
        "The scaling of jet velocities with halo mass is currently not "
        "supported.");
  }

  if (v_jet <= 0.) {
    error(
        "The black_hole_feedback_dv_jet returned a value less than 0. which "
        " is v_jet = %f.",
        v_jet);
  }

  return v_jet;
}

/**
 * @brief Compute the resultant spin of a black hole merger.
 *
 * This implements the fitting formula from Rezzolla et al. (2008).
 * The effects of gravitational waves are ignored.
 *
 * @param bp Pointer to the b-particle data.
 * @param constants Physical constants (in internal units).
 * @param props Properties of the black hole scheme.
 */
__attribute__((always_inline)) INLINE static void merger_spin_evolve(
    struct bpart* bpi, const struct bpart* bpj,
    const struct phys_const* constants) {

  if ((bpj->subgrid_mass <= 0.) || (bpi->subgrid_mass <= 0.)) {
    error(
        "Something went wrong with calculation of spin of a black hole "
        " merger remnant. The black hole masses are %f and %f, instead of  > "
        "0.",
        bpj->subgrid_mass, bpi->subgrid_mass);
  }

  const float m1 = bpi->subgrid_mass;
  const float m2 = bpj->subgrid_mass;
  const float mass_ratio = m2 / m1;
  const float sym_mass_ratio =
      mass_ratio / ((mass_ratio + 1.) * (mass_ratio + 1.));
  const float reduced_mass = m1 * m2 / (m1 + m2);

  const float spin1 = fabsf(bpi->spin);
  const float spin2 = fabsf(bpj->spin);

  if ((spin1 == 0.) || (spin2 == 0.)) {
    error(
        "Something went wrong with calculation of spin of a black hole "
        " merger remnant. The black hole spins are %f and %f, instead of  > 0.",
        spin1, spin2);
  }

  const float spin_vec1[3] = {spin1 * bpi->angular_momentum_direction[0],
                              spin1 * bpi->angular_momentum_direction[1],
                              spin1 * bpi->angular_momentum_direction[2]};
  const float spin_vec2[3] = {spin2 * bpj->angular_momentum_direction[0],
                              spin2 * bpj->angular_momentum_direction[1],
                              spin2 * bpj->angular_momentum_direction[2]};

  const float relative_coordinates[3] = {
      bpj->x[0] - bpi->x[0], bpj->x[1] - bpi->x[1], bpj->x[2] - bpi->x[2]};
  const float relative_velocities[3] = {
      bpj->v[0] - bpi->v[0], bpj->v[1] - bpi->v[1], bpj->v[2] - bpi->v[2]};

  float orbital_angular_momentum[3] = {
      reduced_mass * (relative_coordinates[1] * relative_velocities[2] -
                      relative_coordinates[2] * relative_velocities[1]),
      reduced_mass * (relative_coordinates[2] * relative_velocities[0] -
                      relative_coordinates[0] * relative_velocities[2]),
      reduced_mass * (relative_coordinates[0] * relative_velocities[1] -
                      relative_coordinates[1] * relative_velocities[0])};

  const float orbital_angular_momentum_magnitude =
      sqrtf(orbital_angular_momentum[0] * orbital_angular_momentum[0] +
            orbital_angular_momentum[1] * orbital_angular_momentum[1] +
            orbital_angular_momentum[2] * orbital_angular_momentum[2]);

  if (orbital_angular_momentum_magnitude > 0.) {
    orbital_angular_momentum[0] =
        orbital_angular_momentum[0] / orbital_angular_momentum_magnitude;
    orbital_angular_momentum[1] =
        orbital_angular_momentum[1] / orbital_angular_momentum_magnitude;
    orbital_angular_momentum[2] =
        orbital_angular_momentum[2] / orbital_angular_momentum_magnitude;
  } else {
    orbital_angular_momentum[0] = 0.;
    orbital_angular_momentum[1] = 0.;
    orbital_angular_momentum[2] = 0.;
  }

  const float angle_0 =
      (spin_vec1[0] * spin_vec2[0] + spin_vec1[1] * spin_vec2[1] +
       spin_vec1[2] * spin_vec2[2]) /
      (spin1 * spin2);
  const float angle_1 = (spin_vec1[0] * orbital_angular_momentum[0] +
                         spin_vec1[1] * orbital_angular_momentum[1] +
                         spin_vec1[2] * orbital_angular_momentum[2]) /
                        spin1;
  const float angle_2 = (spin_vec2[0] * orbital_angular_momentum[0] +
                         spin_vec2[1] * orbital_angular_momentum[1] +
                         spin_vec2[2] * orbital_angular_momentum[2]) /
                        spin2;

  const float l =
      -0.129 / (1. + mass_ratio * mass_ratio) * 1. /
          (1. + mass_ratio * mass_ratio) *
          (spin1 * spin1 +
           spin2 * spin2 * mass_ratio * mass_ratio * mass_ratio * mass_ratio +
           2. * spin1 * spin2 * mass_ratio * mass_ratio * angle_0) +
      ((-0.384 * sym_mass_ratio - 0.686) / (1. + mass_ratio * mass_ratio)) *
          (spin1 * angle_1 + spin2 * mass_ratio * mass_ratio * angle_2) +
      3.464 - 3.454 * sym_mass_ratio + 2.353 * sym_mass_ratio * sym_mass_ratio;

#ifdef SWIFT_DEBUG_CHECKS
  if (l < 0.) {
    error(
        "Something went wrong with calculation of spin of a black hole "
        " merger remnant. The l factor is %f, instead of  >= 0.",
        l);
  }
#endif

  float final_spin[3] = {
      1. / (1. + mass_ratio) / (1. + mass_ratio) *
          (spin_vec1[0] + mass_ratio * mass_ratio * spin_vec2[0] +
           mass_ratio * l * orbital_angular_momentum[0]),
      1. / (1. + mass_ratio) / (1. + mass_ratio) *
          (spin_vec1[1] + mass_ratio * mass_ratio * spin_vec2[1] +
           mass_ratio * l * orbital_angular_momentum[1]),
      1. / (1. + mass_ratio) / (1. + mass_ratio) *
          (spin_vec1[2] + mass_ratio * mass_ratio * spin_vec2[2] +
           mass_ratio * l * orbital_angular_momentum[2])};
  const float final_spin_magnitude =
      sqrtf(final_spin[0] * final_spin[0] + final_spin[1] * final_spin[1] +
            final_spin[2] * final_spin[2]);

#ifdef SWIFT_DEBUG_CHECKS
  if (final_spin_magnitude <= 0.) {
    error(
        "Something went wrong with calculation of spin of a black hole "
        " merger remnant. The final spin magnitude is %f, instead of > 0.",
        final_spin_magnitude);
  }
#endif

  final_spin[0] = final_spin[0] / final_spin_magnitude;
  final_spin[1] = final_spin[1] / final_spin_magnitude;
  final_spin[2] = final_spin[2] / final_spin_magnitude;

  bpi->spin = min(final_spin_magnitude, 0.998);
  if (fabsf(bpi->spin) < 0.001) {
    bpi->spin = 0.001;
  }

  bpi->angular_momentum_direction[0] = final_spin[0];
  bpi->angular_momentum_direction[1] = final_spin[1];
  bpi->angular_momentum_direction[2] = final_spin[2];
}

#endif /* SWIFT_SPIN_JET_BLACK_HOLES_SPIN_H */
