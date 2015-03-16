#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <libmints/mints.h>
#include <libqt/qt.h>
#include "libparallel/ParallelPrinter.h"
#include <libchkpt/chkpt.h>

#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugamp {

void amp_write_T2(double ****T2, int length, std::string label, std::string OutFileRMR);

double mp2(boost::shared_ptr<Wavefunction> wfn, boost::shared_ptr<Chkpt> chkpt)
{
  int no = moinfo.no;
  int nv = moinfo.nv;
  int nso = moinfo.nso;
  int nfzc = moinfo.nfzv;
  double ****L = moinfo.L;
  double ****ints = moinfo.ints;
  double ****D2 = moinfo.D2;

  double ****t2_1 = init_4d_array(no, no, nv, nv);
  double ****l2_1 = init_4d_array(no, no, nv, nv);
  moinfo.t2_1 = t2_1;
  moinfo.l2_1 = l2_1;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
          t2_1[i][j][a][b] = ints[i][j][a+no][b+no]/D2[i][j][a][b];
          l2_1[i][j][a][b] = 2.0 * L[i][j][a+no][b+no]/D2[i][j][a][b];
        }

  double emp2 = 0.0;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          emp2 += t2_1[i][j][a][b] * L[i][j][a+no][b+no];

  // Compute virtual natural orbitals
  if(params.fvno == true) {
    if(chkpt->rd_nirreps() > 1)
        throw PSIEXCEPTION("You must add\n\n\tsymmetry c1\n\nto the molecule{} block to compute FVNOs.");

    SharedMatrix Pa(new Matrix("MP2 VV Density", nv, nv));
    double **Pap = Pa->pointer();
    for(int a=0; a < nv; a++)
      for(int b=0; b < nv; b++)
        for(int i=0; i < no; i++)
          for(int j=0; j < no; j++)
            for(int c=0; c < nv; c++)
              Pap[a][b] += l2_1[i][j][a][c] * t2_1[i][j][b][c];

//    Pa->print();
    SharedMatrix Pa_V(new Matrix("MP2 VV Density Eigenvectors", nv, nv));
    SharedVector Pa_v(new Vector("MP2 VV Density Eigenvalues", nv));
    Pa->diagonalize(Pa_V, Pa_v, descending);
//    Pa_V->print();
//    Pa_v->print();

    // Transform SCF MOs to frozen virtual NOs
    double **C = chkpt->rd_scf();

    SharedMatrix Cv(new Matrix("SCF Virtual MOs", nso, nv));
    double **Cvp = Cv->pointer();
    for(int a=0; a < nv; a++)
      for(int p=0; p < nso; p++)
        Cvp[p][a] = C[p][a+no+nfzc];

    SharedMatrix Z(new Matrix("Cv * Pa_V", nso, nv));
    Z->gemm(false, false, 1.0, Cv, Pa_V, 0.0);
//    Z->print();

    // Copy Z over old SCF virtual MOs
    double **Zp = Z->pointer();
    for(int a=0; a < nv; a++)
      for(int p=0; p < nso; p++)
        C[p][a+no+nfzc] = Zp[p][a];

    chkpt->wt_scf(C);
  }

  return emp2;
}

}} // namespace psi::ugamp
