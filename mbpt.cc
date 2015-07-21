#include "mbpt.h"
#include "globals.h"
#include "hamiltonian.h"
#include <boost/shared_ptr.hpp>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <cmath>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include "perturbation.h"

namespace psi {

MBPT::MBPT(boost::shared_ptr<Wavefunction> reference, boost::shared_ptr<Hamiltonian> H, Options &options, boost::shared_ptr<PSIO> psio) : Wavefunction(options, psio)
{
  outfile->Printf("\n");
  outfile->Printf("\t\t\t**************************\n");
  outfile->Printf("\t\t\t*                        *\n");
  outfile->Printf("\t\t\t*         UGA-MP         *\n");
  outfile->Printf("\t\t\t*                        *\n");
  outfile->Printf("\t\t\t**************************\n");
  outfile->Printf("\n");

  if(options.get_str("REFERENCE") != "RHF")
    throw PSIEXCEPTION("Only for use with RHF references determinants.");

  wfn_ = options.get_str("WFN");
  convergence_ = options.get_double("R_CONVERGENCE");
  maxiter_ = options.get_int("MAXITER");
  do_diis_ = options.get_bool("DIIS");
  ooc_ = options.get_bool("OOC");
  if(options.get_str("DERTYPE") == "NONE") dertype_ = 0;
  else if(options.get_str("DERTYPE") == "FIRST") dertype_ = 1;
  fvno_ = options.get_bool("FVNO");
  num_frzv_ = options.get_int("NUM_FRZV");
  occ_tol_ = options.get_double("OCC_TOL");
  spatial_tol_ = options.get_double("SPATIAL_TOL");

  outfile->Printf("\tWave function        = %s\n", wfn_.c_str());
  outfile->Printf("\tMaxiter              = %d\n", maxiter_);
  outfile->Printf("\tConvergence          = %3.1e\n", convergence_);
  outfile->Printf("\tDIIS                 = %s\n", do_diis_ ? "Yes" : "No");
  outfile->Printf("\tOut-of-core          = %s\n", ooc_ ? "Yes" : "No");
  outfile->Printf("\tDertype              = %d\n", dertype_);
  outfile->Printf("\tFrozen-Virtual NO    = %d\n", fvno_);
  outfile->Printf("\tNO Occupation Cutoff = %3.1e\n", occ_tol_);
  outfile->Printf("\tNO <r^2> Cutoff      = %3.1e\n", spatial_tol_);
  outfile->Printf("\tNo. FVNOs            = %d\n", num_frzv_);
  if(occ_tol_ >= 0.0 && spatial_tol_ >= 0.0) 
    outfile->Printf("\tDeleting FVNOs based on both occupation numbers and spatial extent.\n");
  else if(occ_tol_ >= 0.0 && spatial_tol_ < 0.0) 
    outfile->Printf("\tDeleting FVNOs based on occupation numbers.\n");
  else if(occ_tol_ < 0.0 && spatial_tol_ >= 0.0) 
    outfile->Printf("\tDeleting FVNOs based on user input and spatial extent.\n");
  else if(occ_tol_ < 0.0 && spatial_tol_ < 0.0) 
    outfile->Printf("\tDeleting FVNOs based on user input.\n");

  set_reference_wavefunction(reference);
  copy(reference);

  if(fvno_)  // need full virtual space, because all we're doing is computing VNOs
    for(int i=0; i < nirrep_; i++) frzvpi_[i] = 0;

  int nfrzv = 0;
  no_ = nv_ = 0;
  for(int i=0; i < nirrep_; i++) {
    no_ += doccpi_[i] - frzcpi_[i];
    nv_ += nmopi_[i] - doccpi_[i] - frzvpi_[i];
    nfrzv += frzvpi_[i];
  }
  char ** labels = molecule_->irrep_labels();

  outfile->Printf("\n\tReference Wfn Parameters:\n");
  outfile->Printf("\t---------------------------\n");
  outfile->Printf("\tNumber of irreps        = %d\n", nirrep_);
  outfile->Printf("\tNumber of MOs           = %d\n", nmo_);
  outfile->Printf("\tNumber of active MOs    = %d\n", no_+nv_);
  outfile->Printf("\tNumber of active occ    = %d\n", no_);
  outfile->Printf("\tNumber of active vir    = %d\n", nv_);
  outfile->Printf("\tNumber of frozen occ    = %d\n", nfrzc_);
  outfile->Printf("\tNumber of frozen vir    = %d\n\n", nfrzv);
  outfile->Printf("\tLabel\t# MOs\t# FZDC\t# DOCC\t# VIRT\t# FZVR\n");
  outfile->Printf("\t-----\t-----\t------\t------\t------\t------\n");
  for(int i=0; i < nirrep_; i++) {
      outfile->Printf("\t %s\t   %d\t    %d\t    %d\t    %d\t    %d\n",
              labels[i],nmopi_[i],frzcpi_[i],doccpi_[i],nmopi_[i]-doccpi_[i],frzvpi_[i]);
    }
  outfile->Printf("\n\tNuclear Repulsion Energy    = %20.15f\n", molecule_->nuclear_repulsion_energy());
  outfile->Printf( "\tFrozen Core Energy          = %20.15f\n", reference->efzc());
  outfile->Printf( "\tTotal SCF Energy (ref)      = %20.15f\n", reference_wavefunction_->reference_energy());

  for(int i=0; i < nirrep_; i++) free(labels[i]);
  free(labels);

  H_ = H;                 // does this copy properly?
  reference_ = reference; // does this copy properly?

  // Prepare energy denominators
  int no = no_;
  int nv = nv_;
  double **fock = H_->fock_p();
  double ****ints = H_->ints_p();
  double ****L = H_->L_p();

  D1_ = block_matrix(no,nv);
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      D1_[i][a] = fock[i][i] - fock[a+no][a+no];

  D2_ = init_4d_array(no,no,nv,nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          D2_[i][j][a][b] = fock[i][i] + fock[j][j] - fock[a+no][a+no] - fock[b+no][b+no];

  t2_1_ = init_4d_array(no,no,nv,nv);
  l2_1_ = init_4d_array(no,no,nv,nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
          t2_1_[i][j][a][b] = ints[i][j][a+no][b+no]/D2_[i][j][a][b];
          l2_1_[i][j][a][b] = 2.0 * L[i][j][a+no][b+no]/D2_[i][j][a][b];
        }
}

MBPT::~MBPT()
{
  int no = no_;
  int nv = nv_;

  free_block(D1_);
  free_4d_array(D2_, no, no, nv);
  free_4d_array(t2_1_, no, no, nv);
  free_4d_array(l2_1_, no, no, nv);
}

double MBPT::compute_energy() { return 0.0; }

double MBPT::mp2(boost::shared_ptr<Chkpt> chkpt)
{
  int no = no_;
  int nv = nv_;
  double ****L = H_->L_p();

  double ****t2_1 = t2_1_;
  double ****l2_1 = l2_1_;

  double emp2 = 0.0;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          emp2 += t2_1[i][j][a][b] * L[i][j][a+no][b+no];

  // Compute virtual natural orbitals
  if(fvno_ == true) {
    if(nirrep() > 1)
      throw PSIEXCEPTION("You must add\n\n\tsymmetry c1\n\nto the molecule{} block to compute FVNOs.");

    SharedMatrix Pa(new Matrix("MP2 VV Density", nv, nv));
    double **Pap = Pa->pointer();
    for(int a=0; a < nv; a++)
      for(int b=0; b < nv; b++)
        for(int i=0; i < no; i++)
          for(int j=0; j < no; j++)
            for(int c=0; c < nv; c++)
              Pap[a][b] += l2_1[i][j][a][c] * t2_1[i][j][b][c];

    SharedMatrix Pa_V(new Matrix("MP2 VV Density Eigenvectors", nv, nv));
    SharedVector Pa_v(new Vector("MP2 VV Density Eigenvalues", nv));
    Pa->diagonalize(Pa_V, Pa_v, descending);
    // Pa_V transforms from MOV to NOV
    Pa_v->print();

    // Decision time: delete virtuals based on occupation numbers, num_frzv, spatial extent, or
    // combination thereof?

    int nvno=0;
    SharedMatrix Y(new Matrix("FVNOs (SO, NO)", nso_, nv));
    if(spatial_tol_ >= 0.0 && occ_tol_ >= 0.0) { 
      // Need <r^2> integrals for spatial tolerances

      // Prepare the full SO->NO transformer, X
      SharedMatrix C = Ca();
      SharedMatrix Cv = Ca_subset("SO", "ACTIVE_VIR");
      SharedMatrix X(new Matrix("Cv * Pa_V", nso_, nv));
      X->gemm(false, false, 1.0, Cv, Pa_V, 0.0);

      // Create a new reference wfn replacing its MOs with the NOs
      boost::shared_ptr<Wavefunction> newref = reference_;
      SharedMatrix SCF = newref->Ca();
      double **SCFp = SCF->pointer();
      double **Xp = X->pointer();
      for(int p=0; p < nso_; p++)
        for(int a=0; a < nv; a++)
          SCFp[p][a+no+nfrzc_] = Xp[p][a];
      newref->Ca()->set(SCFp);
      newref->Cb()->set(SCFp);

      // Use the new reference to generate the NO-basis r^2 integrals
      boost::shared_ptr<Perturbation> RR(new Perturbation("RR", newref));

      double **xx = RR->prop_p(0,0);
      double **yy = RR->prop_p(1,1);
      double **zz = RR->prop_p(2,2);
      outfile->Printf("MO#  <r^2>\n");
      for(int p=no; p < H_->nact(); p++)
        outfile->Printf("%d %8.4f\n", p-no, -1.0*(xx[p][p]+yy[p][p]+zz[p][p]));

      // Re-organize the NOs to keep only active virtuals
      double *Pa_vp = Pa_v->pointer();
      double **Pa_Vp = Pa_V->pointer();
      double **Yp = Y->pointer();
      for(int p=no; p < H_->nact(); p++) {
        if((-1.0*(xx[p][p]+yy[p][p]+zz[p][p])) > spatial_tol_ || Pa_vp[p] > occ_tol_) {
          for(int q=0; q < nso_; q++) Yp[q][nvno] = Xp[q][p-no];
          nvno++;
        }
      }
    }

    outfile->Printf("\n# Active Virtual NOs  = %d\n", nvno);
    outfile->Printf(  "# Deleted Virtual NOs = %d\n", nv - nvno);

    // Transform VV Fock matrix to NO space
    SharedMatrix FVV_NO(new Matrix("Y^+ * FVV_MO * Y", nvno, nvno));
    double **FVV_NOp = FVV_NO->pointer();
    double **Pa_Vp = Pa_V->pointer();
    double **fock = H_->fock_p();
    for(int a=0; a < nvno; a++)
      for(int b=0; b < nvno; b++)
        for(int c=0; c < nv; c++)
          FVV_NOp[a][b] += fock[c+no][c+no] * Yp[c][a] * Yp[c][b];

    SharedMatrix FVV_V(new Matrix("VV NO Fock Matrix Eigenvectors", nvno, nvno));
    SharedVector FVV_v(new Vector("VV NO Fock Matrix Eigenvalues", nvno));
    FVV_NO->diagonalize(FVV_V, FVV_v);

    // FVV_V transforms from truncated NOV to truncated semicanonical NOV
    double **FVV_Vp = FVV_V->pointer();
    double **Cp = C->pointer();
    for(int p=0; p < nso_; p++)
      for(int a=0; a < nvno; a++) {
        Cp[p][a+no+nfrzc_] = 0.0;
        for(int b=0; b < nvno; b++)
          Cp[p][a+no+nfrzc_] += Yp[p][b] * FVV_Vp[b][a];
      }

    chkpt->wt_scf(Cp);
    Process::environment.wavefunction()->Ca()->set(Cp);
    Process::environment.wavefunction()->Cb()->set(Cp);
  }

  return emp2;
}

double MBPT::mp3()
{
  int no = no_;
  int nv = nv_;
  double ****L = H_->L_p();
  double ****ints = H_->ints_p();
  double ****t2_1 = t2_1_;
  double ****l2_1 = l2_1_;

  double ****X2 = init_4d_array(no, no, nv, nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {

          for(int c=0; c < nv; c++)
            for(int d=0; d < nv; d++)
              X2[i][j][a][b] += 0.5 * t2_1[i][j][c][d] * ints[a+no][b+no][c+no][d+no];

          for(int k=0; k < no; k++)
            for(int l=0; l < no; l++)
              X2[i][j][a][b] += 0.5 * t2_1[k][l][a][b] * ints[i][j][k][l];

          for(int k=0; k < no; k++)
            for(int c=0; c < nv; c++)
              X2[i][j][a][b] += ( t2_1[i][k][a][c] * L[b+no][k][j][c+no]
                                - t2_1[k][j][a][c] * ints[b+no][k][c+no][i]
                                - t2_1[k][i][a][c] * ints[b+no][k][j][c+no] );
        }

  double emp3 = 0.0;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          emp3 += l2_1[i][j][a][b] * X2[i][j][a][b];

  outfile->Printf("\tEMP3-a (corr)  = %20.15f\n", emp3);

  free_4d_array(X2, no, no, nv);

  X2 = init_4d_array(no, no, nv, nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {

          for(int e=0; e < nv; e++)
            for(int f=0; f < nv; f++)
              X2[i][j][a][b] += 0.5 * l2_1[i][j][e][f] * ints[e+no][f+no][a+no][b+no];

          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              X2[i][j][a][b] += 0.5 * l2_1[m][n][a][b] * ints[i][j][m][n];

          for(int m=0; m < no; m++)
            for(int e=0; e < nv; e++)
              X2[i][j][a][b] += ( l2_1[j][m][b][e] * L[e+no][i][m][a+no]
                                - l2_1[m][i][b][e] * ints[j][e+no][m][a+no]
                                - l2_1[m][i][e][b] * ints[e+no][j][m][a+no] );
        }

  emp3 = 0.0;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          emp3 += t2_1[i][j][a][b] * X2[i][j][a][b];

  outfile->Printf("\tEMP3-b (corr)  = %20.15f\n", emp3);

  return emp3;
}

double MBPT::mp4()
{
  int no = no_;
  int nv = nv_;
  double **fock = H_->fock_p();
  double ****L = H_->L_p();
  double ****ints = H_->ints_p();
  double **D1 = D1_;
  double ****D2 = D2_;
  double ****t2_1 = t2_1_;
  double ****l2_1 = l2_1_;

  // Second-order singles amplitudes
  double **t1_2 = block_matrix(no, nv);
  t1_2_ = t1_2;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {

      for(int d=0; d < nv; d++)
        for(int k=0; k < no; k++)
          for(int l=0; l < no; l++)
            t1_2[i][a] -= t2_1[k][l][a][d] * L[k][l][i][d+no];

      for(int c=0; c < nv; c++)
        for(int d=0; d < nv; d++)
          for(int k=0; k < no; k++)
            t1_2[i][a] += t2_1[k][i][c][d] * L[k][a+no][c+no][d+no];

      t1_2[i][a] /= D1[i][a];
    }

  // Second-order doubles amplitudes
  double ****t2_2 = init_4d_array(no, no, nv, nv);
  t2_2_ = t2_2;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {

          for(int c=0; c < nv; c++)
            for(int d=0; d < nv; d++)
              t2_2[i][j][a][b] += t2_1[i][j][c][d] * ints[a+no][b+no][c+no][d+no];

          for(int k=0; k < no; k++)
            for(int l=0; l < no; l++)
              t2_2[i][j][a][b] += t2_1[k][l][a][b] * ints[k][l][i][j];

          for(int k=0; k < no; k++)
            for(int c=0; c < nv; c++)
              t2_2[i][j][a][b] += ( t2_1[i][k][a][c] * L[b+no][k][j][c+no]
                                  + t2_1[j][k][b][c] * L[a+no][k][i][c+no]
                                  - t2_1[k][j][a][c] * ints[b+no][k][c+no][i]
                                  - t2_1[k][i][b][c] * ints[a+no][k][c+no][j]
                                  - t2_1[k][i][a][c] * ints[b+no][k][j][c+no]
                                  - t2_1[k][j][b][c] * ints[a+no][k][i][c+no] );

          t2_2[i][j][a][b] /= D2[i][j][a][b];
        }

  // Second-order triples amplitudes
  double ******t3_2 = init_6d_array(no, no, no, nv, nv, nv);
  t3_2_ = t3_2;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int k=0; k < no; k++)
        for(int a=0; a < nv; a++)
          for(int b=0; b < nv; b++)
            for(int c=0; c < nv; c++) {
              double value = 0.0;

              for(int e=0; e < nv; e++) {
                value +=
                  + ints[i][e+no][a+no][b+no] * t2_1[k][j][c][e]
                  + ints[i][e+no][a+no][c+no] * t2_1[j][k][b][e]
                  + ints[k][e+no][c+no][a+no] * t2_1[j][i][b][e]
                  + ints[k][e+no][c+no][b+no] * t2_1[i][j][a][e]
                  + ints[j][e+no][b+no][c+no] * t2_1[i][k][a][e]
                  + ints[j][e+no][b+no][a+no] * t2_1[k][i][c][e];
              }
              for(int m=0; m < no; m++) {
                value -=
                  + ints[j][k][m][c+no] * t2_1[i][m][a][b]
                  + ints[k][j][m][b+no] * t2_1[i][m][a][c]
                  + ints[i][j][m][b+no] * t2_1[k][m][c][a]
                  + ints[j][i][m][a+no] * t2_1[k][m][c][b]
                  + ints[k][i][m][a+no] * t2_1[j][m][b][c]
                  + ints[i][k][m][c+no] * t2_1[j][m][b][a];
              }

              value /= (fock[i][i] + fock[j][j] + fock[k][k] -
                        fock[a+no][a+no] - fock[b+no][b+no] - fock[c+no][c+no]);

              t3_2[i][j][k][a][b][c] = value;
            }

  // intermediates for MP4-Q
  double ****Ioooo = init_4d_array(no, no, no, no);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int k=0; k < no; k++)
        for(int l=0; l < no; l++)
          for(int c=0; c < nv; c++)
            for(int d=0; d < nv; d++)
              Ioooo[i][j][k][l] += t2_1[i][j][c][d] * ints[k][l][c+no][d+no];

  double ****Ioovv_a = init_4d_array(no, no, nv, nv);
  double ****Ioovv_b = init_4d_array(no, no, nv, nv);
  for(int j=0; j < no; j++)
    for(int k=0; k < no; k++)
      for(int b=0; b < nv; b++)
        for(int c=0; c < nv; c++)
          for(int l=0; l < no; l++)
            for(int d=0; d < nv; d++) {
              Ioovv_a[j][k][b][c] += (t2_1[j][l][b][d] - t2_1[l][j][b][d]) * L[k][l][c+no][d+no];
              Ioovv_b[j][k][b][c] += t2_1[l][j][b][d] * ints[k][l][c+no][d+no];
            }

  double ****Ioovv_c = init_4d_array(no, no, nv, nv);
  for(int i=0; i < no; i++)
    for(int k=0; k < no; k++)
      for(int b=0; b < nv; b++)
        for(int d=0; d < nv; d++)
          for(int c=0; c < nv; c++)
            for(int l=0; l < no; l++)
              Ioovv_c[i][k][b][d] += t2_1[l][i][b][c] * ints[k][l][c+no][d+no];

  double **Ioo = block_matrix(no, no);
  for(int j=0; j < no; j++)
    for(int k=0; k < no; k++)
      for(int l=0; l < no; l++)
        for(int c=0; c < nv; c++)
          for(int d=0; d < nv; d++)
            Ioo[j][k] += t2_1[l][j][c][d] * L[l][k][c+no][d+no];

  double **Ivv = block_matrix(nv, nv);
  for(int b=0; b < nv; b++)
    for(int c=0; c < nv; c++)
      for(int k=0; k < no; k++)
        for(int l=0; l < no; l++)
          for(int d=0; d < nv; d++)
            Ivv[b][c] += t2_1[k][l][b][d] * L[k][l][c+no][d+no];

  double ****S = init_4d_array(no, no, nv, nv);
  double ****D = init_4d_array(no, no, nv, nv);
  double ****T = init_4d_array(no, no, nv, nv);
  double ****Q = init_4d_array(no, no, nv, nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {

          for(int c=0; c < nv; c++)
            S[i][j][a][b] += t1_2[j][c] * ints[a+no][b+no][i][c+no];

          for(int k=0; k < no; k++)
            S[i][j][a][b] -= t1_2[k][b] * ints[a+no][k][i][j];

          for(int c=0; c < nv; c++)
            for(int d=0; d < nv; d++)
              D[i][j][a][b] += 0.5 * t2_2[i][j][c][d] * ints[a+no][b+no][c+no][d+no];

          for(int k=0; k < no; k++)
            for(int l=0; l < no; l++)
              D[i][j][a][b] += 0.5 * t2_2[k][l][a][b] * ints[k][l][i][j];

          for(int k=0; k < no; k++)
            for(int c=0; c < nv; c++)
              D[i][j][a][b] += ( t2_2[i][k][a][c] * L[b+no][k][j][c+no]
                               - t2_2[k][j][a][c] * ints[b+no][k][c+no][i]
                               - t2_2[k][i][a][c] * ints[b+no][k][j][c+no] );

          for(int c=0; c < nv; c++)
            for(int d=0; d < nv; d++)
              for(int k=0; k < no; k++)
                T[i][j][a][b] += ( t3_2[i][j][k][a][c][d] * L[b+no][k][c+no][d+no]
                                 - t3_2[k][j][i][a][c][d] * ints[k][b+no][d+no][c+no] );

          for(int c=0; c < nv; c++)
            for(int k=0; k < no; k++)
              for(int l=0; l < no; l++)
                T[i][j][a][b] -= ( t3_2[i][k][l][a][b][c] * L[k][l][j][c+no]
                                 - t3_2[l][k][i][a][b][c] * ints[k][l][j][c+no] );

          for(int k=0; k < no; k++)
            for(int l=0; l < no; l++)
              Q[i][j][a][b] += 0.5 * t2_1[k][l][a][b] * Ioooo[i][j][k][l];

          for(int c=0; c < nv; c++)
            for(int k=0; k < no; k++) {
              Q[i][j][a][b] += t2_1[i][k][a][c] * Ioovv_a[j][k][b][c];
              Q[i][j][a][b] += 0.5 * t2_1[k][i][a][c] * Ioovv_b[j][k][b][c];
            }

          for(int d=0; d < nv; d++)
            for(int k=0; k < no; k++)
              Q[i][j][a][b] += 0.5 * t2_1[k][j][a][d] * Ioovv_c[i][k][b][d];

          for(int k=0; k < no; k++)
            Q[i][j][a][b] -= t2_1[i][k][a][b] * Ioo[j][k];

          for(int c=0; c < nv; c++)
            Q[i][j][a][b] -= t2_1[i][j][a][c] * Ivv[b][c];
        }

  free_4d_array(Ioooo, no, no, no);
  free_4d_array(Ioovv_a, no, no, nv);
  free_4d_array(Ioovv_b, no, no, nv);
  free_4d_array(Ioovv_c, no, no, nv);
  free_block(Ioo);
  free_block(Ivv);

  double emp4_s = 0.0;
  double emp4_d = 0.0;
  double emp4_t = 0.0;
  double emp4_q = 0.0;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
          emp4_s += l2_1[i][j][a][b] * S[i][j][a][b];
          emp4_d += l2_1[i][j][a][b] * D[i][j][a][b];
          emp4_t += l2_1[i][j][a][b] * T[i][j][a][b];
          emp4_q += l2_1[i][j][a][b] * Q[i][j][a][b];
        }

  free_4d_array(S, no, no, nv);
  free_4d_array(D, no, no, nv);
  free_4d_array(T, no, no, nv);
  free_4d_array(Q, no, no, nv);

  outfile->Printf("\tEMP4-S (corr)  = %20.15f\n", emp4_s);
  outfile->Printf("\tEMP4-D (corr)  = %20.15f\n", emp4_d);
  outfile->Printf("\tEMP4-T (corr)  = %20.15f\n", emp4_t);
  outfile->Printf("\tEMP4-Q (corr)  = %20.15f\n", emp4_q);

  return emp4_s + emp4_d + emp4_t + emp4_q;
}

} // psi
