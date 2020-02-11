
#include <iostream>
#include "MElCalc.hpp"

const MElParticle MElCalc::omegaMeson_ =
    MElParticle(223);
const MElParticle MElCalc::piChMeson_ =
    MElParticle(211);
const MElParticle MElCalc::rho0Meson_ =
  MElParticle(113, &MElCalc::getRho770Width);
const MElParticle MElCalc::rhoChMeson_ =
    MElParticle(213);

double MElCalc::getRho770Width(
    double q2, double massRho770, double widthOnMassRho770, double*) {
  double mass2Rho770 = massRho770 * massRho770;
  return widthOnMassRho770 * mass2Rho770 / q2 *
    pow(p2_pi(q2) / p2_pi(mass2Rho770), 1.5);
}

double MElCalc::p2_pi(double q2) {
  return 0.25 * q2 - piChMeson_.getMass() * piChMeson_.getMass();
}

double MElCalc::getEta2PiMEl2(
    const CFourVector& p_eta,
    const CFourVector& p_pimi,
    const CFourVector& p_pipl) {
  double q2 = (p_pimi + p_pipl).getM2().real();
  std::complex<double> prop = rho0Meson_.getPropagator(q2);
  CFourVector conv = epsilon_conv(p_eta, p_pimi, p_pipl);
  cdouble energy = 0.5 * (p_eta + p_pimi + p_pipl)[3];
  CFourVector l1(0., 0., energy, energy);
  CFourVector l2(0., 0., -energy, energy);
  return electron_current_conv(l1, l2, conv) * norm(prop);
}

double MElCalc::getOmega2PiMEl2(
    const std::pair<CFourVector, CFourVector>& p_pimi,
    const std::pair<CFourVector, CFourVector>& p_pipl,
    const CFourVector& p_pi0) {
  CFourVector p_lepton =
      p_pimi.first + p_pimi.second +
      p_pipl.first + p_pipl.second +
      p_pi0;
  cdouble energy = 0.5 * p_lepton[3];
  CFourVector l1(0., 0., energy, energy);
  CFourVector l2(0., 0., -energy, energy);
  return abs(electron_current_conv(
      l1, l2, getOmega2Pi_A(p_lepton, p_pimi, p_pipl, p_pi0)));
}

CFourVector MElCalc::getOmega2Pi_A(
    const CFourVector& p_lepton,
    const std::pair<CFourVector, CFourVector>& p_pimi,
    const std::pair<CFourVector, CFourVector>& p_pipl,
    const CFourVector& p_pi0) {
  return
      getOmega2Pi_B(p_lepton, p_pimi.first, p_pipl.first, p_pi0) +
      getOmega2Pi_B(p_lepton, p_pimi.first, p_pipl.second, p_pi0) +
      getOmega2Pi_B(p_lepton, p_pimi.second, p_pipl.first, p_pi0) +
      getOmega2Pi_B(p_lepton, p_pimi.second, p_pipl.second, p_pi0);
}

CFourVector MElCalc::getOmega2Pi_B(
    const CFourVector& p_lepton,
    const CFourVector& p_pimi,
    const CFourVector& p_pipl,
    const CFourVector& p_pi0) {
  CFourVector p_omega = p_pipl + p_pimi + p_pi0;
  double p2_omega = p_omega.getM2().real();
  CFourVector vC = getOmega2Pi_C(
      p2_omega, p_pimi, p_pipl, p_pi0);
  return
      p_omega * (p_lepton * vC) +
      (p_lepton * p_omega) * vC;
}

CFourVector MElCalc::getOmega2Pi_C(
    double p2_omega,
    const CFourVector& p_pimi,
    const CFourVector& p_pipl,
    const CFourVector& p_pi0) {
  cdouble prop_omega =
      omegaMeson_.getPropagator(p2_omega);
  CFourVector q_rho0 = p_pipl + p_pimi;
  CFourVector q_rho_minus = p_pimi + p_pi0;
  CFourVector q_rho_plus = p_pipl + p_pi0;
  double q2_rho0 = q_rho0.getM2().real();
  double q2_rho_minus = q_rho_minus.getM2().real();
  double q2_rho_plus = q_rho_plus.getM2().real();
  cdouble prop_rho =
      rho0Meson_.getPropagator(q2_rho0) +
      rhoChMeson_.getPropagator(q2_rho_minus) +
      rhoChMeson_.getPropagator(q2_rho_plus);
  return epsilon_conv(p_pi0, p_pipl, p_pimi) *
      prop_omega * prop_rho;
}
