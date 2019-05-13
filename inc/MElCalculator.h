#ifndef MElCalculator_H
#define MElCalculator_H

#include <utility>
#include "CFourVector.h"
#include "MElParticle.h"

class MElCalculator {
 public:
  static double getEta2PiMEl2(
      const CFourVector&,
      const CFourVector&,
      const CFourVector&);
  static double getOmega2PiMEl2(
    const std::pair<CFourVector, CFourVector>&,
    const std::pair<CFourVector, CFourVector>&,
    const CFourVector&);

 private:
  static const MElParticle omegaMeson_;
  static const MElParticle piChMeson_;
  static const MElParticle rho0Meson_;
  static const MElParticle rhoChMeson_;
  static double getRho770Width(
      double, double, double, double*);
  static double p2_pi(double);
  static CFourVector getOmega2Pi_A(
      const CFourVector&,
      const std::pair<CFourVector, CFourVector>&,
      const std::pair<CFourVector, CFourVector>&,
      const CFourVector&);
  static CFourVector getOmega2Pi_B(
      const CFourVector&,
      const CFourVector&,
      const CFourVector&,
      const CFourVector&);
  static CFourVector getOmega2Pi_C(
      double,
      const CFourVector&,
      const CFourVector&,
      const CFourVector&);
};

#endif
