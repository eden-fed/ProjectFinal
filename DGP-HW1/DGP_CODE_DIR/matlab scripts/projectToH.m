n=size(C,2); %row length

%we need to get from the code:
%a(size of A =set of the sampled boundary), n, R(invanted values), C(the derivative of cauchy coords)
% z0, k, Sigma,sigma
% b(the number of barycenters in tha mesh - num of faces)

cvx_begin
    variable  psi(n) complex;
    variable r(a);
    minimize(sum((r-R)'*(r-R))+sum(abs(C*psi)^2));
    subject to
        abs(C*psi)<=k*r;
        abs(C*psi)<=Sigma-r;
        abs(C*psi)<=k-sigma;
cvx_end