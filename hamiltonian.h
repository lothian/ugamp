#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <boost/shared_ptr.hpp>
#include <libmints/mints.h>
#include <libtrans/integraltransform.h>

namespace psi { namespace ugamp {

class Hamiltonian {
public:
  Hamiltonian(boost::shared_ptr<PSIO>, boost::shared_ptr<Wavefunction>, std::vector<boost::shared_ptr<MOSpace> >, bool, bool);
  virtual ~Hamiltonian();
//  Hamiltonian(const boost::shared_ptr<Hamiltonian> &H);

  double ** fock_p() { return fock_; }
  double **** ints_p() { return ints_; }
  double **** L_p() { return L_; }
  int nact() { return nact_; }

protected:
  int nmo_;
  int nact_;
  int nfzc_;
  int nfzv_;
  int no_;
  int nv_;

  bool ovov_only_;

  double **fock_;
  double ****ints_;
  double ****L_;

}; // Hamiltonian

}} // psi::ugamp

#endif // HAMILTONIAN_H
