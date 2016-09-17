/* Produced by CVXGEN, 2016-09-14 11:06:52 -0400.  */
/* CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2012 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: testsolver.c. */
/* Description: Basic test harness for solver.c. */
#include "solver.h"
Vars vars;
Params params;
Workspace work;
Settings settings;
#define NUMTESTS 0
int main(int argc, char **argv) {
  int num_iters;
#if (NUMTESTS > 0)
  int i;
  double time;
  double time_per;
#endif
  set_defaults();
  setup_indexing();
  load_default_data();
  /* Solve problem instance for the record. */
  settings.verbose = 1;
  num_iters = solve();
#ifndef ZERO_LIBRARY_MODE
#if (NUMTESTS > 0)
  /* Now solve multiple problem instances for timing purposes. */
  settings.verbose = 0;
  tic();
  for (i = 0; i < NUMTESTS; i++) {
    solve();
  }
  time = tocq();
  printf("Timed %d solves over %.3f seconds.\n", NUMTESTS, time);
  time_per = time / NUMTESTS;
  if (time_per > 1) {
    printf("Actual time taken per solve: %.3g s.\n", time_per);
  } else if (time_per > 1e-3) {
    printf("Actual time taken per solve: %.3g ms.\n", 1e3*time_per);
  } else {
    printf("Actual time taken per solve: %.3g us.\n", 1e6*time_per);
  }
#endif
#endif
  return 0;
}
void load_default_data(void) {
  params.l_gz[0] = 0.20319161029830202;
  params.Vg[0] = 0.8325912904724193;
  params.k[0] = -0.8363810443482227;
  params.log_SIGMA[0] = 0.04331042079065206;
  params.A[0] = 1.5717878173906188;
  params.A[1] = 1.5851723557337523;
  params.A[2] = -1.497658758144655;
  params.A[3] = -1.171028487447253;
  params.A[4] = -1.7941311867966805;
  params.A[5] = -0.23676062539745413;
  params.A[6] = -1.8804951564857322;
  params.A[7] = -0.17266710242115568;
  params.A[8] = 0.596576190459043;
  params.A[9] = -0.8860508694080989;
  params.A[10] = 0.7050196079205251;
  params.A[11] = 0.3634512696654033;
  params.A[12] = -1.9040724704913385;
  params.A[13] = 0.23541635196352795;
  params.A[14] = -0.9629902123701384;
  params.A[15] = -0.3395952119597214;
  params.A[16] = -0.865899672914725;
  params.A[17] = 0.7725516732519853;
  params.A[18] = -0.23818512931704205;
  params.A[19] = -1.372529046100147;
  params.A[20] = 0.17859607212737894;
  params.A[21] = 1.1212590580454682;
  params.A[22] = -0.774545870495281;
  params.A[23] = -1.1121684642712744;
  params.A[24] = -0.44811496977740495;
  params.A[25] = 1.7455345994417217;
  params.A[26] = 1.9039816898917352;
  params.A[27] = 0.6895347036512547;
  params.A[28] = 1.6113364341535923;
  params.A[29] = 1.383003485172717;
  params.A[30] = -0.48802383468444344;
  params.A[31] = -1.631131964513103;
  params.A[32] = 0.6136436100941447;
  params.A[33] = 0.2313630495538037;
  params.A[34] = -0.5537409477496875;
  params.A[35] = -1.0997819806406723;
  params.A[36] = -0.3739203344950055;
  params.A[37] = -0.12423900520332376;
  params.A[38] = -0.923057686995755;
  params.A[39] = -0.8328289030982696;
  params.A[40] = -0.16925440270808823;
  params.A[41] = 1.442135651787706;
  params.A[42] = 0.34501161787128565;
  params.A[43] = -0.8660485502711608;
  params.A[44] = -0.8880899735055947;
  params.A[45] = -0.1815116979122129;
  params.A[46] = -1.17835862158005;
  params.A[47] = -1.1944851558277074;
  params.A[48] = 0.05614023926976763;
  params.A[49] = -1.6510825248767813;
  params.B[0] = -0.06565787059365391;
  params.B[1] = -0.5512951504486665;
  params.B[2] = 0.8307464872626844;
  params.B[3] = 0.9869848924080182;
  params.B[4] = 0.7643716874230573;
  params.B[5] = 0.7567216550196565;
  params.B[6] = -0.5055995034042868;
  params.B[7] = 0.6725392189410702;
  params.B[8] = -0.6406053441727284;
  params.B[9] = 0.29117547947550015;
  params.B[10] = -0.6967713677405021;
  params.B[11] = -0.21941980294587182;
  params.B[12] = -1.753884276680243;
  params.B[13] = -1.0292983112626475;
  params.B[14] = 1.8864104246942706;
  params.B[15] = -1.077663182579704;
  params.B[16] = 0.7659100437893209;
  params.B[17] = 0.6019074328549583;
  params.B[18] = 0.8957565577499285;
  params.B[19] = -0.09964555746227477;
  params.B[20] = 0.38665509840745127;
  params.B[21] = -1.7321223042686946;
  params.B[22] = -1.7097514487110663;
  params.B[23] = -1.2040958948116867;
  params.B[24] = -1.3925560119658358;
  params.B[25] = -1.5995826216742213;
  params.B[26] = -1.4828245415645833;
  params.B[27] = 0.21311092723061398;
  params.B[28] = -1.248740700304487;
  params.B[29] = 1.808404972124833;
  params.B[30] = 0.7264471152297065;
  params.B[31] = 0.16407869343908477;
  params.B[32] = 0.8287224032315907;
  params.B[33] = -0.9444533161899464;
  params.B[34] = 1.7069027370149112;
  params.B[35] = 1.3567722311998827;
  params.B[36] = 0.9052779937121489;
  params.B[37] = -0.07904017565835986;
  params.B[38] = 1.3684127435065871;
  params.B[39] = 0.979009293697437;
  params.B[40] = 0.6413036255984501;
  params.B[41] = 1.6559010680237511;
  params.B[42] = 0.5346622551502991;
  params.B[43] = -0.5362376605895625;
  params.B[44] = 0.2113782926017822;
  params.B[45] = -1.2144776931994525;
  params.B[46] = -1.2317108144255875;
  params.B[47] = 0.9026784957312834;
  params.B[48] = 1.1397468137245244;
  params.B[49] = 1.8883934547350631;
}
