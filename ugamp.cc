#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libciomr/libciomr.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include "MOInfo.h"
#include "Params.h"
#include "globals.h"
#include "libparallel/ParallelPrinter.h"

INIT_PLUGIN

using namespace boost;

namespace psi { namespace ugamp {

void title(void);
void get_moinfo(boost::shared_ptr<Wavefunction> wfn, boost::shared_ptr<Chkpt> chkpt);
void integrals(void);
void denom(void);
void cleanup(void);
double mp2(void);
double mp3(void);
double mp4(void);

extern "C" 
int read_options(std::string name, Options& options)
{
  if(name == "UGAMP" || options.read_globals()) {
    options.add_int("PRINT", 1);
    options.add_str("REFERENCE", "RHF");
    options.add_str("WFN", "CCSD");
    options.add_str("DERTYPE", "NONE");
    options.add_int("MAXITER", 100);
    options.add_bool("DIIS", true);
    options.add_double("R_CONVERGENCE", 1e-7);
    options.add_bool("OOC", false);
  }

  return true;
}

extern "C" 
PsiReturnType ugamp(Options& options)
{
  title();
  params.ref = options.get_str("REFERENCE");
  params.wfn = options.get_str("WFN");
  if(options.get_str("DERTYPE") == "NONE") params.dertype = 0;
  else if(options.get_str("DERTYPE") == "FIRST") params.dertype = 1;
  params.convergence = options.get_double("R_CONVERGENCE");
  params.do_diis = options.get_bool("DIIS");
  params.maxiter = options.get_int("MAXITER");
  params.ooc = options.get_bool("OOC");

  outfile->Printf("\tWave function  = %s\n", params.wfn.c_str());
  outfile->Printf("\tReference      = %s\n", params.ref.c_str());
  outfile->Printf("\tComputation    = %s\n", params.dertype ? "Gradient" : "Energy");
  outfile->Printf("\tMaxiter        = %d\n", params.maxiter);
  outfile->Printf("\tConvergence    = %3.1e\n", params.convergence);
  outfile->Printf("\tDIIS           = %s\n", params.do_diis ? "Yes" : "No");

  boost::shared_ptr<PSIO> psio(_default_psio_lib_);
  boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
  if(!wfn) throw PSIEXCEPTION("SCF has not been run yet!");
  boost::shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));

  get_moinfo(wfn, chkpt);
  integrals();
  denom();

  moinfo.emp2 = mp2();
  outfile->Printf("\tEMP2 (corr)    = %20.15f\n", moinfo.emp2);
  outfile->Printf("\tEMP2           = %20.15f\n", moinfo.emp2 + moinfo.escf);

  moinfo.emp3 = mp3();
  outfile->Printf("\tEMP3 (corr)    = %20.15f\n", moinfo.emp3);
  outfile->Printf("\tEMP3           = %20.15f\n", moinfo.emp2 + moinfo.emp3 + moinfo.escf);

  moinfo.emp4 = mp4();
  outfile->Printf("\tEMP4 (corr)    = %20.15f\n", moinfo.emp4);
  outfile->Printf("\tEMP3           = %20.15f\n", moinfo.emp2 + moinfo.emp3 + moinfo.emp4 + moinfo.escf);

  cleanup();

  return Success;
}

void title(void)
{
  outfile->Printf("\n");
  outfile->Printf("\t\t\t**************************\n");
  outfile->Printf("\t\t\t*                        *\n");
  outfile->Printf("\t\t\t*         UGA-MP         *\n");
  outfile->Printf("\t\t\t*                        *\n");
  outfile->Printf("\t\t\t**************************\n");
  outfile->Printf("\n");
}

}} // End namespaces

