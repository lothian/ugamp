#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>

#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugamp {

void amp_write_T2(double ****T2, int length, std::string label, std::string OutFileRMR);

double mp2(void)
{
  int no = moinfo.no;
  int nv = moinfo.nv;
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

  return emp2;
}

}} // namespace psi::ugamp
