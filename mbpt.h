#ifndef MBPT_H
#define MBPT_H

#include "hamiltonian.h"
#include <libmints/mints.h>
#include <boost/shared_ptr.hpp>
#include <libchkpt/chkpt.h>

namespace psi { namespace ugamp {

class MBPT: public Wavefunction {
public:
  MBPT(std::string wfn, boost::shared_ptr<Wavefunction> reference, boost::shared_ptr<Hamiltonian> H,
                 Options& options, boost::shared_ptr<PSIO> psio, bool fvno);
  virtual ~MBPT();

protected:
  std::string wfn_;     // wfn type (MP2, MP3, etc.)

  int no_;  // Number of active occupied MOs
  int nv_;  // Number of active virtual MOs

  boost::shared_ptr<Hamiltonian> H_; // integrals and Fock matrix
  boost::shared_ptr<Wavefunction> reference_;
  
  // Energy denominators
  double **D1_;
  double ****D2_;

  // Ground-state T amplitudes
  double **t1_2_;    // Second-order T1
  double ****t2_1_;  // First-order T2
  double ****l2_1_;  // First-order L2
  double ****t2_2_;  // Second-order T2
  double ******t3_2_; // Second-order T3

public:
  double compute_energy();

  double mp2();
  double mp3();
  double mp4();

  SharedMatrix DVV(); // compute the VV block of the MP2 onepdm

  int nv() { return nv_; }
  int no() { return no_; }
}; // MBPT

}} // psi::ugamp

#endif // MBPT_H
