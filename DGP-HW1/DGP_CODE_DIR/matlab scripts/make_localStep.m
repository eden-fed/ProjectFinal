% Produced by CVXGEN, 2016-09-14 11:06:51 -0400.
% CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com.
% The code in this file is Copyright (C) 2006-2012 Jacob Mattingley.
% CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial
% applications without prior written permission from Jacob Mattingley.

% Filename: make_csolve.m.
% Description: Calls mex to generate the csolve mex file.
%mex -v localStep.cpp ldl.c matrix_support.c solver.c util.c
mex localStep.c ldl.c matrix_support.c solver.c util.c