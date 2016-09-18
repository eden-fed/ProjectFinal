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
%     out_Real_log_gz >= A(1)*out_abs_Vg + B(1)
%     out_Real_log_gz >= A(2)*out_abs_Vg + B(2)
%     out_Real_log_gz >= A(3)*out_abs_Vg + B(3)
%     out_Real_log_gz >= A(4)*out_abs_Vg + B(4)
%     out_Real_log_gz >= A(5)*out_abs_Vg + B(5)
%
% with variables
% out_Real_log_gz   1 x 1
% out_abs_Vg   1 x 1
%
% and parameters
%        A   5 x 1
%        B   5 x 1
% in_Real_log_gz   1 x 1
% in_abs_Vg   1 x 1
%        k   1 x 1
% log_SIGMA   1 x 1
%
% Note:
%   - Check status.converged, which will be 1 if optimization succeeded.
%   - You don't have to specify settings if you don't want to.
%   - To hide output, use settings.verbose = 0.
%   - To change iterations, use settings.max_iters = 20.
%   - You may wish to compare with cvxsolve to check the solver is correct.
%
% Specify params.A, ..., params.log_SIGMA, then run
%   [vars, status] = csolve(params, settings)
% Produced by CVXGEN, 2016-09-18 02:13:45 -0400.
% CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com.
% The code in this file is Copyright (C) 2006-2012 Jacob Mattingley.
% CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial
% applications without prior written permission from Jacob Mattingley.

% Filename: csolve.m.
% Description: Help file for the Matlab solver interface.
