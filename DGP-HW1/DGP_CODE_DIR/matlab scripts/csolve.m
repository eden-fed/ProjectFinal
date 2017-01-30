% csolve  Solves a custom quadratic program very rapidly.
%
% [vars, status] = csolve(params, settings)
%
% solves the convex optimization problem
%
%   minimize(square(in_Real_log_gz - out_Real_log_gz) + square(in_abs_Vg - out_abs_Vg))
%   subject to
%     out_abs_Vg >= 0
%     out_abs_Vg <= k
%     out_abs_Vg <= log_SIGMA - out_Real_log_gz
%     out_Real_log_gz >= m*out_abs_Vg + b
%
% with variables
% out_Real_log_gz   1 x 1
% out_abs_Vg   1 x 1
%
% and parameters
%        b   1 x 1
% in_Real_log_gz   1 x 1
% in_abs_Vg   1 x 1
%        k   1 x 1
% log_SIGMA   1 x 1
%        m   1 x 1
%
% Note:
%   - Check status.converged, which will be 1 if optimization succeeded.
%   - You don't have to specify settings if you don't want to.
%   - To hide output, use settings.verbose = 0.
%   - To change iterations, use settings.max_iters = 20.
%   - You may wish to compare with cvxsolve to check the solver is correct.
%
% Specify params.b, ..., params.m, then run
%   [vars, status] = csolve(params, settings)
% Produced by CVXGEN, 2016-12-28 05:27:42 -0500.
% CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com.
% The code in this file is Copyright (C) 2006-2012 Jacob Mattingley.
% CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial
% applications without prior written permission from Jacob Mattingley.

% Filename: csolve.m.
% Description: Help file for the Matlab solver interface.
