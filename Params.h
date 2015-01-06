#ifndef _psi_psi_ugamp_params_h
#define _psi_psi_ugamp_params_h

#include <string>

namespace psi { namespace ugamp {

struct Params {
    std::string ref;
    std::string wfn;
    int dertype;
    int maxiter;
    double convergence;
    int do_diis;        // DIIS boolean
    int ooc;            // out-of-core boolean
};

}} // namespace devel::ugamp

#endif // _psi_psi_ugamp_params_h
