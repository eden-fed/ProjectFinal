% clc
a=size(Ctag,1);
n=size(Ctag,2);
k=(SIGMA-sigma)/(SIGMA+sigma);
%*************step 2:Evaluate gz and gz_gag******************
% cageVerteciesAfterMap
% cageVerteciesB4Map
% 
deltaS=cageVerteciesB4Map-circshift(cageVerteciesB4Map,1);
deltaD=cageVerteciesAfterMap-circshift(cageVerteciesAfterMap,1);

deltaS=circshift(deltaS,size(deltaS,1)-1);
deltaD=circshift(deltaD,size(deltaD,1)-1);

%gpuArray(deltaS);
%gpuArray(deltaD);

%[ out1,out2 ] = arrayfun(@calcGzAndGzGag,deltaS,deltaD);
gz=0.5*((deltaD.*(abs(deltaS)+abs(deltaD)))./(deltaS.*abs(deltaD)));
gz_gag=(-0.5)*((deltaD.*(abs(deltaS)-abs(deltaD)))./(deltaS.*abs(deltaD)));

%gz=gather(out1);
%gz_gag=gather(out2);

%this shouldnt be done on the gpu
gz_enc=repelem(gz,NumOfVerticesInEdges);
gz_gag_enc=repelem(gz_gag,NumOfVerticesInEdges);



%*************step 3:solve 13 - obtain r, psi(z)******************
Cz0=C_sizeM(Z0index,:);

%R=abs(Ctag_tempAN*cageVerteciesAfterMap);
abs_gz=abs(gz_enc);
conj_gs_gag=conj(gz_gag_enc);


% abs_gz=ones(a,1);
% conj_gs_gag=zeros(a,1);

tic
cvx_begin
    variable  psai(n) complex;
    variable r(a);
    minimize((r-abs_gz)'*(r-abs_gz)+(Ctag*psai-conj_gs_gag)'*(Ctag*psai-conj_gs_gag));
    subject to
		Cz0*psai == zeros(size(Cz0*psai));
        abs(Ctag*psai)<=k*r;
        abs(Ctag*psai)<=SIGMA-r;
        abs(Ctag*psai)<=r-sigma;
cvx_end

time1=toc;
%*************step 4:solve 15 - obtain l(z)******************
tic
cvx_begin
    variable  l(n) complex;
    minimize((real(C_sizeA*l)-log(r))'*(real(C_sizeA*l)-log(r)));
    subject to
        imag(Cz0*l)==mean(angle(cageVerteciesB4Map-mean(cageVerteciesB4Map))-angle(cageVerteciesAfterMap-mean(cageVerteciesAfterMap)));
      %  imag(Cz0*l)==0;
cvx_end
time2=toc;
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
partialCalc_gpu=gpuArray(PHItag(endIndices)) + gpuArray(PHItag(startIndices));
integral_on_edges_gpu=partialCalc_gpu.*edgeVectors_gpu;
integral_on_edges=gather(integral_on_edges_gpu);
% % 
% partialCalc=(PHItag(endIndices) + PHItag(startIndices))./2;
% integral_on_edges=partialCalc.*edgeVectors;

% find the integral on all the spanning tree
PHI_Z0=PHI_Z0+0.000000000000000001i;
PHI = treeCumSum(uint32(Z0index), PHI_Z0, integral_on_edges, startIndices, endIndices);

%*************step 7:find f******************
PSI=C_sizeM*psai;
f=PHI+conj(PSI);
