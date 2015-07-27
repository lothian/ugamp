#include "mbpt.h"
#include "perturbation.h"

#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <libtrans/integraltransform.h>
#include "globals.h"

INIT_PLUGIN

using namespace boost;

namespace psi { namespace ugamp {

extern "C" 
int read_options(std::string name, Options& options)
{
  if(name == "UGAMP" || options.read_globals()) {
    options.add_int("PRINT", 1);
    options.add_str("REFERENCE", "RHF");
    options.add_str("WFN", "CCSD", "MP2 MP3 MP4 CCSD CCSD_T");

    options.add_bool("FVNO", false); // compute MP2 frozen-virtual NOs
    options.add_str("FREEZE_TYPE", "MULTI_ORB", "MULTI_ORB SINGLE_ORB OCCUPATION SPATIAL HYBRID");
    options.add_int("NUM_FRZV", 0); // number of FVNOs or which single orb to delete
    options.add_double("OCC_TOL", 0.0); // delete FVNOs below cutoff
    options.add_double("SPATIAL_TOL", 0.0); // only delete FVNOs below R^2 cutoff
  }

  return true;
}

extern "C" 
PsiReturnType ugamp(Options& options)
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

  std::string wfn = options.get_str("WFN");
  std::string freeze_type = options.get_str("FREEZE_TYPE");
  bool fvno = options.get_bool("FVNO");
  double occ_tol = options.get_double("OCC_TOL");
  double spatial_tol = options.get_double("SPATIAL_TOL");
  int num_frzv = options.get_int("NUM_FRZV");

  outfile->Printf("\tWave function        = %s\n", wfn.c_str());
  outfile->Printf("\tFrozen-Virtual NO    = %s\n", fvno ? "Yes" : "No");
  if(fvno) {
    outfile->Printf("\tFreeze Type          = %s\n", freeze_type.c_str());
    outfile->Printf("\tNO Occupation Cutoff = %3.1e\n", occ_tol);
    outfile->Printf("\tNO <r^2> Cutoff      = %3.1e\n", spatial_tol);
    outfile->Printf("\t#/Index FVNO         = %d\n", num_frzv);
  }

  boost::shared_ptr<PSIO> psio(_default_psio_lib_);
  boost::shared_ptr<Wavefunction> ref = Process::environment.wavefunction();
  if(!ref) throw PSIEXCEPTION("SCF has not been run yet!");

  std::vector<boost::shared_ptr<MOSpace> > spaces;
  spaces.push_back(MOSpace::all);

  boost::shared_ptr<Hamiltonian> H(new Hamiltonian(psio, ref, spaces, options.get_bool("FVNO")));
  boost::shared_ptr<MBPT> mbpt(new MBPT(wfn, ref, H, options, psio, options.get_bool("FVNO")));

  boost::shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));

  double eref=0.0, emp2=0.0, emp3=0.0, emp4=0.0;
  eref = mbpt->reference_energy();
 
  if(options.get_str("WFN") == "MP2" || options.get_str("WFN") == "MP3" || options.get_str("WFN") == "MP4") {
    emp2 = mbpt->mp2();
    outfile->Printf("\tEMP2 (corr)    = %20.15f\n", emp2);
    outfile->Printf("\tEMP2           = %20.15f\n", emp2 + eref);
  }

  if(options.get_str("WFN") == "MP3" || options.get_str("WFN") == "MP4") {
    emp3 = mbpt->mp3();
    outfile->Printf("\tEMP3 (corr)    = %20.15f\n", emp3);
    outfile->Printf("\tEMP3           = %20.15f\n", emp2 + emp3 + eref);
  }

  if(options.get_str("WFN") == "MP4") {
    emp4 = mbpt->mp4();
    outfile->Printf("\tEMP4 (corr)    = %20.15f\n", emp4);
    outfile->Printf("\tEMP4           = %20.15f\n", emp2 + emp3 + emp4 + eref);
  }

  if(fvno) {
    if(ref->nirrep() > 1)
      throw PSIEXCEPTION("You must add\n\n\tsymmetry c1\n\nto the molecule{} block to compute FVNOs.");

    int nv = mbpt->nv();
    int no = mbpt->no();
    int nso = ref->nso();
    Dimension frzvpi = ref->frzvpi();

    SharedMatrix DVV = mbpt->DVV();
    SharedMatrix T_VMO_2_VNO(new Matrix("MP2 VV Density Eigenvectors", nv, nv));
    SharedVector DVV_eps(new Vector("MP2 VV Density Eigenvalues", nv));
    DVV->diagonalize(T_VMO_2_VNO, DVV_eps, descending);
    DVV_eps->print();

    // Compute the spatial extent of each NO, regardless of FREEZE_TYPE
    SharedMatrix C = ref->Ca();
    SharedMatrix Cv = ref->Ca_subset("SO", "ACTIVE_VIR");
    SharedMatrix T_SO_2_VNO(new Matrix("Cv * T_VMO_2_VNO", nso, nv));
    T_SO_2_VNO->gemm(false, false, 1.0, Cv, T_VMO_2_VNO, 0.0);

    // Create a new reference wfn replacing its MOs with the NOs
    boost::shared_ptr<Wavefunction> newref = ref;
    SharedMatrix SCF = newref->Ca();
    double **SCFp = SCF->pointer();
    double **T_SO_2_VNOp = T_SO_2_VNO->pointer();
    for(int p=0; p < nso; p++)
      for(int a=0; a < nv; a++)
        SCFp[p][a+no+ref->nfrzc()] = T_SO_2_VNOp[p][a];

    // Use the new reference to generate the NO-basis r^2 integrals
    boost::shared_ptr<Perturbation> RR(new Perturbation("RR", newref, true));

    SharedVector RR_NO(new Vector("NO <R^2> Values", nv));
    double **xx = RR->prop_p(0,0);
    double **yy = RR->prop_p(1,1);
    double **zz = RR->prop_p(2,2);
    for(int p=no; p < H->nact(); p++) RR_NO->set(p-no, -1.0*(xx[p][p]+yy[p][p]+zz[p][p]));
    RR_NO->print();

    // (a) Build boolean vector that identifies active VNOs
    int nvno=0;
    std::vector<bool> Frozen(nv);
    if(freeze_type == "SINGLE_ORB") {
      if(num_frzv > nv) throw PSIEXCEPTION("Chosen VNO index too large.");
      Frozen[nv-num_frzv] = true; // num_frzv = 1 means lowest occupation VNO
      nvno = (num_frzv == 0) ? nv : nv - 1; // in case the user set num_frzv_ = 0, so no frozen orbs
      frzvpi[0] = (num_frzv == 0) ? 0 : 1;
    }
    else if(freeze_type == "MULTI_ORB") {
      if(num_frzv > nv) throw PSIEXCEPTION("Too many frozen virtuals requested.");
      for(int p=nv-num_frzv; p < nv; p++) Frozen[p] = true;
      nvno = nv - num_frzv;
      frzvpi[0] = num_frzv;
    }
    else if(freeze_type == "OCCUPATION") {
      num_frzv = 0;
      for(int p=0; p < nv; p++) { if(DVV_eps->get(p) < occ_tol) { Frozen[p] = true; num_frzv++; } }
      nvno = nv - num_frzv;
      frzvpi[0] = num_frzv;
    }
    else if(freeze_type == "SPATIAL") {
      num_frzv = 0;
      for(int p=0; p < nv; p++) { if(RR_NO->get(p) < spatial_tol) { Frozen[p] = true; num_frzv++; } }
      nvno = nv - num_frzv;
      frzvpi[0] = num_frzv;
    }
    else if(freeze_type == "HYBRID") {
      num_frzv = 0;
      for(int p=0; p < nv; p++) {
        if(RR_NO->get(p) < spatial_tol && DVV_eps->get(p) < occ_tol) { Frozen[p] = true; num_frzv++; }
      }
      nvno = nv - num_frzv;
      frzvpi[0] = num_frzv;
    }

    if(freeze_type == "SINGLE_ORB")
      outfile->Printf("Index of frozen virtual orbital: %d\n", num_frzv);
    else
      outfile->Printf("Number of frozen virtual orbitals: %d\n", num_frzv);
    outfile->Printf("Number of active virtual orbitals: %d\n", nvno);

    // (b) Build T_VMO_2_FVNO matrix by copying over only active VNOs (columns) from full matrix
    SharedMatrix T_VMO_2_FVNO(new Matrix("VMO to FVNO Transform", nv, nvno));
    double **T_VMO_2_FVNOp = T_VMO_2_FVNO->pointer();
    double **T_VMO_2_VNOp = T_VMO_2_VNO->pointer();
    int FVNO_offset = 0;
    for(int a=0; a < nv; a++) {
      if(!Frozen[a]) {
        for(int p=0; p < nv; p++) T_VMO_2_FVNOp[p][FVNO_offset] = T_VMO_2_VNOp[p][a];
        FVNO_offset++;
      }
    }
    if(FVNO_offset != nvno) throw PSIEXCEPTION("Number of active VNOs not equal user-chosen number.");

    // (c) Build T_SO_2_FVNO = Cv * T_VMO_2_FVNO (SO x FVNO)
    SharedMatrix T_SO_2_FVNO(new Matrix("SO to FVNO Transform", nso, nvno));
    T_SO_2_FVNO->gemm(false, false, 1.0, Cv, T_VMO_2_FVNO, 0.0);

    // (d) Fock_FVNO = T_VMO_2_FVNO^+ Fock_MO T_VMO_2_FVNO (FVNO x FVNO)
    SharedMatrix Fock_FVNO(new Matrix("T_VMO_2_FVNO^+ Fock_MO T_VMO_2_FVNO", nvno, nvno));
    double **Fock_FVNOp = Fock_FVNO->pointer();
    double **fock = H->fock_p();
    for(int a=0; a < nvno; a++)
      for(int b=0; b < nvno; b++) {
        Fock_FVNOp[a][b] = 0.0;
        for(int c=0; c < nv; c++)
          Fock_FVNOp[a][b] += T_VMO_2_FVNOp[c][a] * fock[c+no][c+no] * T_VMO_2_FVNOp[c][b];
      }

    // (e) diagonalize Fock_FVNO matrix (evecs T_FVNO_2_semiFVNO)
    SharedMatrix T_FVNO_2_semiFVNO(new Matrix("FVNO to semicanon. FVNO Transform", nvno, nvno));
    SharedVector Fock_FVNO_eps(new Vector("semicanon. FVNO Eigenvalues", nvno));
    Fock_FVNO->diagonalize(T_FVNO_2_semiFVNO, Fock_FVNO_eps);

    // (f) Build T_SO_2_semiFVNO = T_SO_2_FVNO * T_FVNO_2_semiFVNO
    SharedMatrix T_SO_2_semiFVNO(new Matrix("SO to semicanon. FVNO Transform", nso, nvno));
    double **T_FVNO_2_semiFVNOp = T_FVNO_2_semiFVNO->pointer();
    double **T_SO_2_FVNOp = T_SO_2_FVNO->pointer();
    double **Cp = C->pointer();
    for(int p=0; p < nso; p++)
      for(int a=0; a < nvno; a++) {
        Cp[p][a+no+ref->nfrzc()] = 0.0;
        for(int b=0; b < nvno; b++)
          Cp[p][a+no+ref->nfrzc()] += T_SO_2_FVNOp[p][b] * T_FVNO_2_semiFVNOp[b][a];
      }

    // Copy new info into necessary locations for subsequent CC computations
    chkpt->wt_scf(Cp);
    Process::environment.wavefunction()->Ca()->set(Cp);
    Process::environment.wavefunction()->Cb()->set(Cp);

    chkpt->wt_frzvpi(frzvpi);
    Process::environment.wavefunction()->set_frzvpi(frzvpi);
    if(Process::environment.options["FROZEN_UOCC"].size() == 0)
      Process::environment.options["FROZEN_UOCC"].add(0);
    else
      Process::environment.options["FROZEN_UOCC"][0].assign(frzvpi[0]);
  }

  return Success;
}

}} // End psi::ugamp

