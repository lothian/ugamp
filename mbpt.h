#ifndef MBPT_H
#define MBPT_H

#include "hamiltonian.h"
#include <libmints/mints.h>
#include <boost/shared_ptr.hpp>
#include <libchkpt/chkpt.h>

namespace psi {

class MBPT: public Wavefunction {
public:
  MBPT(boost::shared_ptr<Wavefunction> reference,
                 boost::shared_ptr<Hamiltonian> H,
                 Options &options, boost::shared_ptr<PSIO> psio);
  virtual ~MBPT();

protected:
  std::string wfn_;     // wfn type (MP2, MP3, etc.)
  double convergence_;  // conv. on RMS residual change between iterations
  int maxiter_;         // maximum number of iterations
  bool do_diis_;        // use DIIS algorithms?
  bool ooc_;            // Use out-of-core algorithms?
  int dertype_;         // Gradient level
  bool fvno_;           // Frozen-virtual NO calculation?
  int num_frzv_;        // Number of frozen-virtuals

  int no_;  // Number of active occupied MOs
  int nv_;  // Number of active virtual MOs

  boost::shared_ptr<Hamiltonian> H_; // integrals and Fock matrix
  
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
  int maxiter() { return maxiter_; }
  double convergence() { return convergence_; }
  std::string wfn() { return wfn_; } 
  bool do_diis() { return do_diis_; }
  bool ooc() { return ooc_; }
  int dertype() { return dertype_; }
  bool fvno() { return fvno_; }
  int num_frzv() { return num_frzv_; }

  double compute_energy();

  double mp2(boost::shared_ptr<Chkpt>);
  double mp3();
  double mp4();
}; // MBPT

} // psi

#endif // MBPT_H
