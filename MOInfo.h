#ifndef _psi_psi_ugamp_moinfo_h
#define _psi_psi_ugamp_moinfo_h

namespace psi { namespace ugamp {

struct MOInfo {
  int nmo;          /* # molecular orbitals */
  int nso;          /* # symmetry orbitals */
  int nact;         /* # active orbitals */
  int no;           /* # occupied orbitals */
  int nv;           /* # unoccupied orbitals */
  int nfzc;         /* # frozen core orbitals */
  int nfzv;         /* # frozen virtual orbitals */
  double enuc;      /* nuclear repulsion energy */
  double escf;      /* SCF energy */
  double efzc;      /* frozen core energy */
  double **fock;    /* f(p,q) = h(p,q) + sum_i ( 2<pi|qi> - <pi|iq> ) */
  double ****ints;  /* <pq|rs> */
  double ****L;     /* 2<pq|rs> - <pq|sr> */

  /* perturbed energies */
  double emp2;      /* MP2 energy */
  double emp3;      /* MP3 energy */
  double emp4;      /* MP3 energy */

  // T-amplitude quantities
  double **D1;      /* one-electron denominators */
  double ****D2;    /* two-electron denominators */
  double **t1_2;    /* second-order t1 amplitudes */
  double ****t2_1;  /* first-order t2 amplitudes */
  double ****l2_1;  /* first-order t2 amplitudes with L integrals */
  double ****t2_2;  /* second-order t2 amplitudes */
  double ******t3_2;  /* third-order t2 amplitudes */
};

}} // namespace psi::ugamp

#endif // _psi_psi_ugamp_moinfo_h
