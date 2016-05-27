
%*************step 3:solve 13 - obtain r, psi(z)******************
%we need to get from the code: Cz0=cauchy coordinates for z0 (the size is 1xn)
%b(the number of barycenters in the mesh = num of faces)

a=size(Ctag,1);
n=size(Ctag,2);

Cz0=C(Z0index,:);

R=2.*ones(a,1);

cvx_begin
    variable  psai(n) complex;
    variable r(a);
    minimize((r-R)'*(r-R));
    subject to
		Cz0*psai=0;
        abs(Ctag*psai)<=k*r;
        abs(Ctag*psai)<=SIGMA-r;
        abs(Ctag*psai)<=k-sigma;
cvx_end


%*************step 4:solve 15 - obtain l(z)******************
%we need to get from the code: z0 / indZ0(the index of z0)
cvx_begin
    variable  l(n) complex;
    minimize((real(C*l)-log(r))'*(real(C*l)-log(r)));
    subject to
        imag(Cz0*l)=0;
cvx_end

%*************step 5:find derivative of phi(z)******************
PHItag=exp(C*l);

%*************step 6:find phi(z) - integral******************
%spanning tree build : [disc, pred, closed] = graphtraverse(adjacencyGraph, anchorVertexIndex, 'Directed', false, 'Method', 'BFS');


%*************step 7:find f******************
