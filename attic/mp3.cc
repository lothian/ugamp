#include <string>
#include <psi4-dec.h>
#include "libparallel/ParallelPrinter.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugamp {

double mp3(void)
{
  int no = moinfo.no;
  int nv = moinfo.nv;
  double ****L = moinfo.L;
  double ****ints = moinfo.ints;
  double ****t2_1 = moinfo.t2_1;
  double ****l2_1 = moinfo.l2_1;

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

  std::string out;
  boost::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            boost::shared_ptr<OutFile>(new OutFile(out)));

  printer->Printf("\tEMP3-a (corr)  = %20.15f\n", emp3);

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

  printer->Printf("\tEMP3-b (corr)  = %20.15f\n", emp3);

  return emp3;
}

}} // namespace psi::ugamp
