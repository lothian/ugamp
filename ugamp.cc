#include "mbpt.h"

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
    options.add_str("DERTYPE", "NONE");
    options.add_int("MAXITER", 100);
    options.add_bool("DIIS", true);
    options.add_double("R_CONVERGENCE", 1e-7);
    options.add_bool("OOC", false);
    options.add_bool("FVNO", false);
    options.add_int("NUM_FRZV", 0);
  }

  return true;
}

extern "C" 
PsiReturnType ugamp(Options& options)
{
  boost::shared_ptr<PSIO> psio(_default_psio_lib_);
  boost::shared_ptr<Wavefunction> ref = Process::environment.wavefunction();
  if(!ref) throw PSIEXCEPTION("SCF has not been run yet!");

  std::vector<boost::shared_ptr<MOSpace> > spaces;
  spaces.push_back(MOSpace::all);

  boost::shared_ptr<Hamiltonian> H(new Hamiltonian(psio, ref, spaces, options.get_bool("FVNO")));
  boost::shared_ptr<MBPT> mbpt(new MBPT(ref, H, options, psio));

  boost::shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));

  double eref=0.0, emp2=0.0, emp3=0.0, emp4=0.0;
  eref = mbpt->reference_energy();
 
  if(mbpt->wfn() == "MP2" || mbpt->wfn() == "MP3" || mbpt->wfn() == "MP4") {
    emp2 = mbpt->mp2(chkpt);
    outfile->Printf("\tEMP2 (corr)    = %20.15f\n", emp2);
    outfile->Printf("\tEMP2           = %20.15f\n", emp2 + eref);
  }

  if(mbpt->wfn() == "MP3" || mbpt->wfn() == "MP4") {
    emp3 = mbpt->mp3();
    outfile->Printf("\tEMP3 (corr)    = %20.15f\n", emp3);
    outfile->Printf("\tEMP3           = %20.15f\n", emp2 + emp3 + eref);
  }

  if(mbpt->wfn() == "MP4") {
    emp4 = mbpt->mp4();
    outfile->Printf("\tEMP4 (corr)    = %20.15f\n", emp4);
    outfile->Printf("\tEMP4           = %20.15f\n", emp2 + emp3 + emp4 + eref);
  }

  return Success;
}

}} // End namespaces

