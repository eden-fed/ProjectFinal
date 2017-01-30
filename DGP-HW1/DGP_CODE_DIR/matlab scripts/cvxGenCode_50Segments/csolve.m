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
%     out_Real_log_gz >= A(6)*out_abs_Vg + B(6)
%     out_Real_log_gz >= A(7)*out_abs_Vg + B(7)
%     out_Real_log_gz >= A(8)*out_abs_Vg + B(8)
%     out_Real_log_gz >= A(9)*out_abs_Vg + B(9)
%     out_Real_log_gz >= A(10)*out_abs_Vg + B(10)
%     out_Real_log_gz >= A(11)*out_abs_Vg + B(11)
%     out_Real_log_gz >= A(12)*out_abs_Vg + B(12)
%     out_Real_log_gz >= A(13)*out_abs_Vg + B(13)
%     out_Real_log_gz >= A(14)*out_abs_Vg + B(14)
%     out_Real_log_gz >= A(15)*out_abs_Vg + B(15)
%     out_Real_log_gz >= A(16)*out_abs_Vg + B(16)
%     out_Real_log_gz >= A(17)*out_abs_Vg + B(17)
%     out_Real_log_gz >= A(18)*out_abs_Vg + B(18)
%     out_Real_log_gz >= A(19)*out_abs_Vg + B(19)
%     out_Real_log_gz >= A(20)*out_abs_Vg + B(20)
%     out_Real_log_gz >= A(21)*out_abs_Vg + B(21)
%     out_Real_log_gz >= A(22)*out_abs_Vg + B(22)
%     out_Real_log_gz >= A(23)*out_abs_Vg + B(23)
%     out_Real_log_gz >= A(24)*out_abs_Vg + B(24)
%     out_Real_log_gz >= A(25)*out_abs_Vg + B(25)
%     out_Real_log_gz >= A(26)*out_abs_Vg + B(26)
%     out_Real_log_gz >= A(27)*out_abs_Vg + B(27)
%     out_Real_log_gz >= A(28)*out_abs_Vg + B(28)
%     out_Real_log_gz >= A(29)*out_abs_Vg + B(29)
%     out_Real_log_gz >= A(30)*out_abs_Vg + B(30)
%     out_Real_log_gz >= A(31)*out_abs_Vg + B(31)
%     out_Real_log_gz >= A(32)*out_abs_Vg + B(32)
%     out_Real_log_gz >= A(33)*out_abs_Vg + B(33)
%     out_Real_log_gz >= A(34)*out_abs_Vg + B(34)
%     out_Real_log_gz >= A(35)*out_abs_Vg + B(35)
%     out_Real_log_gz >= A(36)*out_abs_Vg + B(36)
%     out_Real_log_gz >= A(37)*out_abs_Vg + B(37)
%     out_Real_log_gz >= A(38)*out_abs_Vg + B(38)
%     out_Real_log_gz >= A(39)*out_abs_Vg + B(39)
%     out_Real_log_gz >= A(40)*out_abs_Vg + B(40)
%     out_Real_log_gz >= A(41)*out_abs_Vg + B(41)
%     out_Real_log_gz >= A(42)*out_abs_Vg + B(42)
%     out_Real_log_gz >= A(43)*out_abs_Vg + B(43)
%     out_Real_log_gz >= A(44)*out_abs_Vg + B(44)
%     out_Real_log_gz >= A(45)*out_abs_Vg + B(45)
%     out_Real_log_gz >= A(46)*out_abs_Vg + B(46)
%     out_Real_log_gz >= A(47)*out_abs_Vg + B(47)
%     out_Real_log_gz >= A(48)*out_abs_Vg + B(48)
%     out_Real_log_gz >= A(49)*out_abs_Vg + B(49)
%     out_Real_log_gz >= A(50)*out_abs_Vg + B(50)
%
% with variables
% out_Real_log_gz   1 x 1
% out_abs_Vg   1 x 1
%
% and parameters
%        A  50 x 1
%        B  50 x 1
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
% Produced by CVXGEN, 2016-11-24 03:33:00 -0500.
% CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com.
% The code in this file is Copyright (C) 2006-2012 Jacob Mattingley.
% CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial
% applications without prior written permission from Jacob Mattingley.

% Filename: csolve.m.
% Description: Help file for the Matlab solver interface.
