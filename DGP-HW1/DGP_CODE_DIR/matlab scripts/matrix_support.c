/* Produced by CVXGEN, 2016-12-28 05:27:42 -0500.  */
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
  lhs[3] = -rhs[0]*(-1)-rhs[1]*(params.m[0]);
}
void multbymGT(double *lhs, double *rhs) {
  lhs[0] = -rhs[2]*(1)-rhs[3]*(-1);
  lhs[1] = -rhs[0]*(-1)-rhs[1]*(1)-rhs[2]*(1)-rhs[3]*(params.m[0]);
}
void multbyP(double *lhs, double *rhs) {
  /* TODO use the fact that P is symmetric? */
  /* TODO check doubling / half factor etc. */
  lhs[0] = rhs[0]*(2);
  lhs[1] = rhs[1]*(2);
}
void fillq(void) {
  work.q[0] = -2*params.in_Real_log_gz[0];
  work.q[1] = -2*params.in_abs_Vg[0];
}
void fillh(void) {
  work.h[0] = 0;
  work.h[1] = params.k[0];
  work.h[2] = params.log_SIGMA[0];
  work.h[3] = -params.b[0];
}
void fillb(void) {
}
void pre_ops(void) {
  work.quad_382403260416[0] = params.in_Real_log_gz[0]*params.in_Real_log_gz[0];
  work.quad_800024121344[0] = params.in_abs_Vg[0]*params.in_abs_Vg[0];
}
