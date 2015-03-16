#include <string>
#include <psi4-dec.h>
#include <libciomr/libciomr.h>
#include "MOInfo.h"
#define EXTERN
#include "globals.h"
#include "libparallel/ParallelPrinter.h"

namespace psi { namespace ugamp {

double mp4(void)
{
  int no = moinfo.no;
  int nv = moinfo.nv;
  double **fock = moinfo.fock;
  double ****L = moinfo.L;
  double ****ints = moinfo.ints;
  double **D1 = moinfo.D1;
  double ****D2 = moinfo.D2;
  double ****t2_1 = moinfo.t2_1;
  double ****l2_1 = moinfo.l2_1;

  // Second-order singles amplitudes
  double **t1_2 = block_matrix(no, nv);
  moinfo.t1_2 = t1_2;
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
  moinfo.t2_2 = t2_2;
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
  moinfo.t3_2 = t3_2;
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

  std::string out;
  boost::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            boost::shared_ptr<OutFile>(new OutFile(out)));

  printer->Printf("\tEMP4-S (corr)  = %20.15f\n", emp4_s);
  printer->Printf("\tEMP4-D (corr)  = %20.15f\n", emp4_d);
  printer->Printf("\tEMP4-T (corr)  = %20.15f\n", emp4_t);
  printer->Printf("\tEMP4-Q (corr)  = %20.15f\n", emp4_q);

  return emp4_s + emp4_d + emp4_t + emp4_q;
}

}} // namespace psi::ugamp
