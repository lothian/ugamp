/* global variables for newT */

#ifndef _psi_psi_ugamp_globals_h
#define _psi_psi_ugamp_globals_h

namespace psi { namespace ugamp {

double ***init_3d_array(int, int, int);
void free_3d_array(double ***, int, int);

double ****init_4d_array(int, int, int, int);
void free_4d_array(double ****, int, int, int);

double ******init_6d_array(int, int, int, int, int, int);
void free_6d_array(double ******, int, int, int, int, int);

}} // namespace psi::ugamp

#endif // _psi_psi_ugamp_globals_h
