
%we need to get from the code:
%a(size of A =set of the sampled boundary), n, R(invanted values), Ctag(the derivative of cauchy coords)
% z0, k, Sigma,sigma
% b(the number of barycenters in tha mesh - num of faces)

cvx_begin
    variable  psi(n) complex;
    variable r(a);
    minimize((r-R)'*(r-R)+sum(abs(Ctag*psi)^2));
    subject to
        abs(Ctag*psi)<=k*r;
        abs(Ctag*psi)<=Sigma-r;
        abs(Ctag*psi)<=k-sigma;
cvx_end