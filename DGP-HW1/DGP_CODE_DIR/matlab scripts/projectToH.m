% clearvars -except s edgeVectors_gpu
% clc
%*************step 3:solve 13 - obtain r, psi(z)******************

a=size(Ctag,1);
n=size(Ctag,2);

%Z0=0.495826+1i*0.741235;
% Z0=0.407+0.426*1i;

Cz0=C_sizeM(Z0index,:);

%R=2.*ones(a,1);
%R=abs(C_tempAN*cageVerteciesAfterMap);
R=ones(a,1);

cvx_begin
    variable  psai(n) complex;
    variable r(a);
    minimize((r-R)'*(r-R)+(Ctag*psai)'*(Ctag*psai));
    subject to
		Cz0*psai == zeros(size(Cz0*psai));
        abs(Ctag*psai)<=k*r;
        abs(Ctag*psai)<=SIGMA-r;
        abs(Ctag*psai)<=k-sigma;
cvx_end

%*************step 4:solve 15 - obtain l(z)******************

cvx_begin
    variable  l(n) complex;
    minimize((real(C_sizeA*l)-log(r))'*(real(C_sizeA*l)-log(r)));
    subject to
        %imag(Cz0*l)==mean(angle(cageVerteciesB4Map-mean(cageVerteciesB4Map))-angle(cageVerteciesAfterMap-mean(cageVerteciesAfterMap)));
        imag(Cz0*l)==0;
cvx_end

%*************step 5:find derivative of phi(z)******************
PHItag=exp(C_sizeM*l);

%*************step 6:find phi(z) - integral******************
%we need to send :cageVerteciesAfterMap,cageVerteciesB4Map, startIndices, endIndices, integral_on_edges(or edgeVectors_gpu)
if(exist('treeCumSum', 'file') ~= 3)
    fullPathName = which('treeCumSum.cpp');
    [folderName, fileName, ext] = fileparts(fullPathName);
    mex(fullPathName, '-outdir', folderName);
end

PHI_Z0=Z0+mean(cageVerteciesAfterMap)-mean(cageVerteciesB4Map);

%calc the integral on the edges
partialCalc_gpu=(gpuArray(PHItag(endIndices)) + gpuArray(PHItag(startIndices)))./2;
integral_on_edges_gpu=partialCalc_gpu.*edgeVectors_gpu;
integral_on_edges=gather(integral_on_edges_gpu);
% % 
% partialCalc=(PHItag(endIndices) + PHItag(startIndices))./2;
% integral_on_edges=partialCalc.*edgeVectors;

% find the integral on all the spanning tree
PHI = treeCumSum(uint32(Z0index), PHI_Z0, integral_on_edges, startIndices, endIndices);

%*************step 7:find f******************
PSI=C_sizeM*psai;
f=PHI+conj(PSI);
