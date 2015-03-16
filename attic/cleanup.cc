#include <libciomr/libciomr.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugamp {

void cleanup(void)
{
  int no = moinfo.no; 
  int nv = moinfo.nv;
  int nact = moinfo.nact;

  free_block(moinfo.fock);
  free_4d_array(moinfo.ints, nact, nact, nact);
  free_4d_array(moinfo.L, nact, nact, nact);
  free_block(moinfo.D1);
  free_4d_array(moinfo.D2, no, no, nv);

  // MP2- and MP3-related quantities
  free_4d_array(moinfo.t2_1, no, no, nv);

  // MP4-related quantities
//  free_block(moinfo.t1_2);
//  free_4d_array(moinfo.t2_2, no, no, nv);
//  free_6d_array(moinfo.t3_2, no, no, no, nv, nv);
}

}} // namespace psi::ugamp
