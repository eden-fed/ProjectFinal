
%we need to get from the code:
%a(size of A =set of the sampled boundary), n, R(invanted values), Ctag(the derivative of cauchy coords)
% z0, k, Sigma,sigma
% b(the number of barycenters in the mesh - num of faces)
a=size(Ctag,1);
n=size(Ctag,2);
R=2.*ones(a,1);

cvx_begin
    variable  psai(n) complex;
    variable r(a);
    minimize((r-R)'*(r-R));
    subject to
        abs(Ctag*psai)<=k*r;
        abs(Ctag*psai)<=Sigma-r;
        abs(Ctag*psai)<=k-sigma;
cvx_end


%we need to get from the code:
%a(size of A =set of the sampled boundary), n, r(from the previous cvx), C(the cauchy coords)
% z0, indZ0(the index of z0)

cvx_begin
    variable  l(n) complex;
    minimize((real(C*l)-log(r))'*(real(C*l)-log(r)));
    subject to
        imag(C(indZ0,:)*l)=0;
cvx_end