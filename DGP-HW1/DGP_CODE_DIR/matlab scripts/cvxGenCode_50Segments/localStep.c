//this is a MEX file for Matlab.
#define ZERO_LIBRARY_MODE
#include "mex.h"
#include "solver.h"
//#include <iostream>
Vars vars;
Params params;
Workspace work;
Settings settings;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	size_t m,n,i;
	double *src;
	double *dest;
	int steps;
	int j;
	double* global_abs_vg;
	double * global_real_lg;
	double* local_abs_vg ;
	double* local_real_lg ;
	if(nrhs != 6)
	{
		mexErrMsgIdAndTxt("MyToolbox:localStep:nrhs","Five inputs required.");
	}
	if(nlhs!=2) {
		mexErrMsgIdAndTxt("MyToolbox:localStep:nlhs","Two output required.");
	}
	m=mxGetM(prhs[0]);
	if(!(m > 1 &&  mxGetN(prhs[0]) == 1))
	{
		mexErrMsgIdAndTxt("MyToolbox:localStep:wrongInput", "First argument must be a column vector.");//A
	}

	if(!(mxGetM(prhs[1]) > 1 &&  mxGetN(prhs[1]) == 1) && m==mxGetM(prhs[1]))
	{
		mexErrMsgIdAndTxt("MyToolbox:localStep:wrongInput", "Second argument must be a column vector, and the same size as the first argument.");//B
	}
	n=mxGetM(prhs[2]);
	if(!(n > 1 && mxGetN(prhs[2]) == 1))
	{
		mexErrMsgIdAndTxt("MyToolbox:localStep:wrongInput", "Third argument must be a column vector.");//abs(Vg)
	}
	if(!(mxGetM(prhs[3]) > 1 && mxGetN(prhs[3]) == 1) && n==mxGetM(prhs[3]))
	{
		mexErrMsgIdAndTxt("MyToolbox:localStep:wrongInput", "Forth argument must be a column vector, and the same size as the third argument.");//real(l_gz)
	}
	if(!(mxGetM(prhs[4]) == 1 &&  mxGetN(prhs[4]) == 1) )
	{
		mexErrMsgIdAndTxt("MyToolbox:localStep:wrongInput", "Fifth argument must be a scalar.");//k
	}
	if(!(mxGetM(prhs[5]) == 1 &&  mxGetN(prhs[5]) == 1) )
	{
		mexErrMsgIdAndTxt("MyToolbox:localStep:wrongInput", "Sixth argument must be a scalar.");//log(SIGMA)
	}
	


	set_defaults();  // Set basic algorithm parameters.
    setup_indexing();
//	params.A = mxGetPr(prhs[0]);
	dest = params.A;
	src = mxGetPr(prhs[0]);
	for (j = 0; j < m; j++)
		*dest++ = *src++;

//	params.B = mxGetPr(prhs[1]);
	dest = params.B;
	src = mxGetPr(prhs[1]);
	for (j = 0; j < m; j++)
		*dest++ = *src++;

	*params.k = mxGetScalar(prhs[4]);
	*params.log_SIGMA = mxGetScalar(prhs[5]);

	global_abs_vg= mxGetPr(prhs[2]);
	global_real_lg= mxGetPr(prhs[3]);

	plhs[0] = mxCreateDoubleMatrix((mwSize)(n), 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix((mwSize)(n), 1, mxREAL);

	local_abs_vg = mxGetPr(plhs[0]);
	local_real_lg = mxGetPr(plhs[1]);

	for (i = 0; i < n; i++) {
		*params.in_Real_log_gz = global_real_lg[i];
		*params.in_abs_Vg = global_abs_vg[i];
		steps = solve();
	/*	if(!work.converged)
			mexErrMsgIdAndTxt("MyToolbox:localStep:error", "the probles did not converged.");*/
			/* For profiling purposes, allow extra silent solves if desired. */
		settings.verbose = 0;
		local_abs_vg[i] = *vars.out_abs_Vg;
		local_real_lg[i] = *vars.out_Real_log_gz;
	}
}
