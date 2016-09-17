% csolve  Solves a custom quadratic program very rapidly.
%
% [vars, status] = csolve(params, settings)
%
% solves the convex optimization problem
%
%   minimize(quad_form(l_gz - R_l_gz, eye(1)) + quad_form(Vg - abs_Vg, eye(1)))
%   subject to
%     abs_Vg >= 0
%     abs_Vg <= k
%     abs_Vg <= log_SIGMA - R_l_gz
%     R_l_gz >= A(1)*abs_Vg + B(1)
%     R_l_gz >= A(2)*abs_Vg + B(2)
%     R_l_gz >= A(3)*abs_Vg + B(3)
%     R_l_gz >= A(4)*abs_Vg + B(4)
%     R_l_gz >= A(5)*abs_Vg + B(5)
%     R_l_gz >= A(6)*abs_Vg + B(6)
%     R_l_gz >= A(7)*abs_Vg + B(7)
%     R_l_gz >= A(8)*abs_Vg + B(8)
%     R_l_gz >= A(9)*abs_Vg + B(9)
%     R_l_gz >= A(10)*abs_Vg + B(10)
%     R_l_gz >= A(11)*abs_Vg + B(11)
%     R_l_gz >= A(12)*abs_Vg + B(12)
%     R_l_gz >= A(13)*abs_Vg + B(13)
%     R_l_gz >= A(14)*abs_Vg + B(14)
%     R_l_gz >= A(15)*abs_Vg + B(15)
%     R_l_gz >= A(16)*abs_Vg + B(16)
%     R_l_gz >= A(17)*abs_Vg + B(17)
%     R_l_gz >= A(18)*abs_Vg + B(18)
%     R_l_gz >= A(19)*abs_Vg + B(19)
%     R_l_gz >= A(20)*abs_Vg + B(20)
%     R_l_gz >= A(21)*abs_Vg + B(21)
%     R_l_gz >= A(22)*abs_Vg + B(22)
%     R_l_gz >= A(23)*abs_Vg + B(23)
%     R_l_gz >= A(24)*abs_Vg + B(24)
%     R_l_gz >= A(25)*abs_Vg + B(25)
%     R_l_gz >= A(26)*abs_Vg + B(26)
%     R_l_gz >= A(27)*abs_Vg + B(27)
%     R_l_gz >= A(28)*abs_Vg + B(28)
%     R_l_gz >= A(29)*abs_Vg + B(29)
%     R_l_gz >= A(30)*abs_Vg + B(30)
%     R_l_gz >= A(31)*abs_Vg + B(31)
%     R_l_gz >= A(32)*abs_Vg + B(32)
%     R_l_gz >= A(33)*abs_Vg + B(33)
%     R_l_gz >= A(34)*abs_Vg + B(34)
%     R_l_gz >= A(35)*abs_Vg + B(35)
%     R_l_gz >= A(36)*abs_Vg + B(36)
%     R_l_gz >= A(37)*abs_Vg + B(37)
%     R_l_gz >= A(38)*abs_Vg + B(38)
%     R_l_gz >= A(39)*abs_Vg + B(39)
%     R_l_gz >= A(40)*abs_Vg + B(40)
%     R_l_gz >= A(41)*abs_Vg + B(41)
%     R_l_gz >= A(42)*abs_Vg + B(42)
%     R_l_gz >= A(43)*abs_Vg + B(43)
%     R_l_gz >= A(44)*abs_Vg + B(44)
%     R_l_gz >= A(45)*abs_Vg + B(45)
%     R_l_gz >= A(46)*abs_Vg + B(46)
%     R_l_gz >= A(47)*abs_Vg + B(47)
%     R_l_gz >= A(48)*abs_Vg + B(48)
%     R_l_gz >= A(49)*abs_Vg + B(49)
%     R_l_gz >= A(50)*abs_Vg + B(50)
%
% with variables
%   R_l_gz   1 x 1
%   abs_Vg   1 x 1
%
% and parameters
%        A  50 x 1
%        B  50 x 1
%       Vg   1 x 1
%        k   1 x 1
%     l_gz   1 x 1
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
% Produced by CVXGEN, 2016-09-14 11:06:50 -0400.
% CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com.
% The code in this file is Copyright (C) 2006-2012 Jacob Mattingley.
% CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial
% applications without prior written permission from Jacob Mattingley.

% Filename: csolve.m.
% Description: Help file for the Matlab solver interface.
