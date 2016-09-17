/* Produced by CVXGEN, 2016-09-14 11:06:51 -0400.  */
/* CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2012 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: matrix_support.c. */
/* Description: Support functions for matrix multiplication and vector filling. */
#include "solver.h"
void multbymA(double *lhs, double *rhs) {
}
void multbymAT(double *lhs, double *rhs) {
  lhs[0] = 0;
  lhs[1] = 0;
}
void multbymG(double *lhs, double *rhs) {
  lhs[0] = -rhs[1]*(-1);
  lhs[1] = -rhs[1]*(1);
  lhs[2] = -rhs[0]*(1)-rhs[1]*(1);
  lhs[3] = -rhs[0]*(-1)-rhs[1]*(params.A[0]);
  lhs[4] = -rhs[0]*(-1)-rhs[1]*(params.A[1]);
  lhs[5] = -rhs[0]*(-1)-rhs[1]*(params.A[2]);
  lhs[6] = -rhs[0]*(-1)-rhs[1]*(params.A[3]);
  lhs[7] = -rhs[0]*(-1)-rhs[1]*(params.A[4]);
  lhs[8] = -rhs[0]*(-1)-rhs[1]*(params.A[5]);
  lhs[9] = -rhs[0]*(-1)-rhs[1]*(params.A[6]);
  lhs[10] = -rhs[0]*(-1)-rhs[1]*(params.A[7]);
  lhs[11] = -rhs[0]*(-1)-rhs[1]*(params.A[8]);
  lhs[12] = -rhs[0]*(-1)-rhs[1]*(params.A[9]);
  lhs[13] = -rhs[0]*(-1)-rhs[1]*(params.A[10]);
  lhs[14] = -rhs[0]*(-1)-rhs[1]*(params.A[11]);
  lhs[15] = -rhs[0]*(-1)-rhs[1]*(params.A[12]);
  lhs[16] = -rhs[0]*(-1)-rhs[1]*(params.A[13]);
  lhs[17] = -rhs[0]*(-1)-rhs[1]*(params.A[14]);
  lhs[18] = -rhs[0]*(-1)-rhs[1]*(params.A[15]);
  lhs[19] = -rhs[0]*(-1)-rhs[1]*(params.A[16]);
  lhs[20] = -rhs[0]*(-1)-rhs[1]*(params.A[17]);
  lhs[21] = -rhs[0]*(-1)-rhs[1]*(params.A[18]);
  lhs[22] = -rhs[0]*(-1)-rhs[1]*(params.A[19]);
  lhs[23] = -rhs[0]*(-1)-rhs[1]*(params.A[20]);
  lhs[24] = -rhs[0]*(-1)-rhs[1]*(params.A[21]);
  lhs[25] = -rhs[0]*(-1)-rhs[1]*(params.A[22]);
  lhs[26] = -rhs[0]*(-1)-rhs[1]*(params.A[23]);
  lhs[27] = -rhs[0]*(-1)-rhs[1]*(params.A[24]);
  lhs[28] = -rhs[0]*(-1)-rhs[1]*(params.A[25]);
  lhs[29] = -rhs[0]*(-1)-rhs[1]*(params.A[26]);
  lhs[30] = -rhs[0]*(-1)-rhs[1]*(params.A[27]);
  lhs[31] = -rhs[0]*(-1)-rhs[1]*(params.A[28]);
  lhs[32] = -rhs[0]*(-1)-rhs[1]*(params.A[29]);
  lhs[33] = -rhs[0]*(-1)-rhs[1]*(params.A[30]);
  lhs[34] = -rhs[0]*(-1)-rhs[1]*(params.A[31]);
  lhs[35] = -rhs[0]*(-1)-rhs[1]*(params.A[32]);
  lhs[36] = -rhs[0]*(-1)-rhs[1]*(params.A[33]);
  lhs[37] = -rhs[0]*(-1)-rhs[1]*(params.A[34]);
  lhs[38] = -rhs[0]*(-1)-rhs[1]*(params.A[35]);
  lhs[39] = -rhs[0]*(-1)-rhs[1]*(params.A[36]);
  lhs[40] = -rhs[0]*(-1)-rhs[1]*(params.A[37]);
  lhs[41] = -rhs[0]*(-1)-rhs[1]*(params.A[38]);
  lhs[42] = -rhs[0]*(-1)-rhs[1]*(params.A[39]);
  lhs[43] = -rhs[0]*(-1)-rhs[1]*(params.A[40]);
  lhs[44] = -rhs[0]*(-1)-rhs[1]*(params.A[41]);
  lhs[45] = -rhs[0]*(-1)-rhs[1]*(params.A[42]);
  lhs[46] = -rhs[0]*(-1)-rhs[1]*(params.A[43]);
  lhs[47] = -rhs[0]*(-1)-rhs[1]*(params.A[44]);
  lhs[48] = -rhs[0]*(-1)-rhs[1]*(params.A[45]);
  lhs[49] = -rhs[0]*(-1)-rhs[1]*(params.A[46]);
  lhs[50] = -rhs[0]*(-1)-rhs[1]*(params.A[47]);
  lhs[51] = -rhs[0]*(-1)-rhs[1]*(params.A[48]);
  lhs[52] = -rhs[0]*(-1)-rhs[1]*(params.A[49]);
}
void multbymGT(double *lhs, double *rhs) {
  lhs[0] = -rhs[2]*(1)-rhs[3]*(-1)-rhs[4]*(-1)-rhs[5]*(-1)-rhs[6]*(-1)-rhs[7]*(-1)-rhs[8]*(-1)-rhs[9]*(-1)-rhs[10]*(-1)-rhs[11]*(-1)-rhs[12]*(-1)-rhs[13]*(-1)-rhs[14]*(-1)-rhs[15]*(-1)-rhs[16]*(-1)-rhs[17]*(-1)-rhs[18]*(-1)-rhs[19]*(-1)-rhs[20]*(-1)-rhs[21]*(-1)-rhs[22]*(-1)-rhs[23]*(-1)-rhs[24]*(-1)-rhs[25]*(-1)-rhs[26]*(-1)-rhs[27]*(-1)-rhs[28]*(-1)-rhs[29]*(-1)-rhs[30]*(-1)-rhs[31]*(-1)-rhs[32]*(-1)-rhs[33]*(-1)-rhs[34]*(-1)-rhs[35]*(-1)-rhs[36]*(-1)-rhs[37]*(-1)-rhs[38]*(-1)-rhs[39]*(-1)-rhs[40]*(-1)-rhs[41]*(-1)-rhs[42]*(-1)-rhs[43]*(-1)-rhs[44]*(-1)-rhs[45]*(-1)-rhs[46]*(-1)-rhs[47]*(-1)-rhs[48]*(-1)-rhs[49]*(-1)-rhs[50]*(-1)-rhs[51]*(-1)-rhs[52]*(-1);
  lhs[1] = -rhs[0]*(-1)-rhs[1]*(1)-rhs[2]*(1)-rhs[3]*(params.A[0])-rhs[4]*(params.A[1])-rhs[5]*(params.A[2])-rhs[6]*(params.A[3])-rhs[7]*(params.A[4])-rhs[8]*(params.A[5])-rhs[9]*(params.A[6])-rhs[10]*(params.A[7])-rhs[11]*(params.A[8])-rhs[12]*(params.A[9])-rhs[13]*(params.A[10])-rhs[14]*(params.A[11])-rhs[15]*(params.A[12])-rhs[16]*(params.A[13])-rhs[17]*(params.A[14])-rhs[18]*(params.A[15])-rhs[19]*(params.A[16])-rhs[20]*(params.A[17])-rhs[21]*(params.A[18])-rhs[22]*(params.A[19])-rhs[23]*(params.A[20])-rhs[24]*(params.A[21])-rhs[25]*(params.A[22])-rhs[26]*(params.A[23])-rhs[27]*(params.A[24])-rhs[28]*(params.A[25])-rhs[29]*(params.A[26])-rhs[30]*(params.A[27])-rhs[31]*(params.A[28])-rhs[32]*(params.A[29])-rhs[33]*(params.A[30])-rhs[34]*(params.A[31])-rhs[35]*(params.A[32])-rhs[36]*(params.A[33])-rhs[37]*(params.A[34])-rhs[38]*(params.A[35])-rhs[39]*(params.A[36])-rhs[40]*(params.A[37])-rhs[41]*(params.A[38])-rhs[42]*(params.A[39])-rhs[43]*(params.A[40])-rhs[44]*(params.A[41])-rhs[45]*(params.A[42])-rhs[46]*(params.A[43])-rhs[47]*(params.A[44])-rhs[48]*(params.A[45])-rhs[49]*(params.A[46])-rhs[50]*(params.A[47])-rhs[51]*(params.A[48])-rhs[52]*(params.A[49]);
}
void multbyP(double *lhs, double *rhs) {
  /* TODO use the fact that P is symmetric? */
  /* TODO check doubling / half factor etc. */
  lhs[0] = rhs[0]*(2);
  lhs[1] = rhs[1]*(2);
}
void fillq(void) {
  work.q[0] = -2*params.l_gz[0];
  work.q[1] = -2*params.Vg[0];
}
void fillh(void) {
  work.h[0] = 0;
  work.h[1] = params.k[0];
  work.h[2] = params.log_SIGMA[0];
  work.h[3] = -params.B[0];
  work.h[4] = -params.B[1];
  work.h[5] = -params.B[2];
  work.h[6] = -params.B[3];
  work.h[7] = -params.B[4];
  work.h[8] = -params.B[5];
  work.h[9] = -params.B[6];
  work.h[10] = -params.B[7];
  work.h[11] = -params.B[8];
  work.h[12] = -params.B[9];
  work.h[13] = -params.B[10];
  work.h[14] = -params.B[11];
  work.h[15] = -params.B[12];
  work.h[16] = -params.B[13];
  work.h[17] = -params.B[14];
  work.h[18] = -params.B[15];
  work.h[19] = -params.B[16];
  work.h[20] = -params.B[17];
  work.h[21] = -params.B[18];
  work.h[22] = -params.B[19];
  work.h[23] = -params.B[20];
  work.h[24] = -params.B[21];
  work.h[25] = -params.B[22];
  work.h[26] = -params.B[23];
  work.h[27] = -params.B[24];
  work.h[28] = -params.B[25];
  work.h[29] = -params.B[26];
  work.h[30] = -params.B[27];
  work.h[31] = -params.B[28];
  work.h[32] = -params.B[29];
  work.h[33] = -params.B[30];
  work.h[34] = -params.B[31];
  work.h[35] = -params.B[32];
  work.h[36] = -params.B[33];
  work.h[37] = -params.B[34];
  work.h[38] = -params.B[35];
  work.h[39] = -params.B[36];
  work.h[40] = -params.B[37];
  work.h[41] = -params.B[38];
  work.h[42] = -params.B[39];
  work.h[43] = -params.B[40];
  work.h[44] = -params.B[41];
  work.h[45] = -params.B[42];
  work.h[46] = -params.B[43];
  work.h[47] = -params.B[44];
  work.h[48] = -params.B[45];
  work.h[49] = -params.B[46];
  work.h[50] = -params.B[47];
  work.h[51] = -params.B[48];
  work.h[52] = -params.B[49];
}
void fillb(void) {
}
void pre_ops(void) {
  work.quad_772795076608[0] = params.l_gz[0]*params.l_gz[0];
  work.quad_890218008576[0] = params.Vg[0]*params.Vg[0];
}
