% clc
a=size(Ctag,1);
n=size(Ctag,2);
%k=(SIGMA-sigma)/(SIGMA+sigma);
%lambda=1;

%*************step 2:Evaluate gz and gz_gag******************
deltaS=cageVerteciesB4Map-circshift(cageVerteciesB4Map,1);
deltaD=cageVerteciesAfterMap-circshift(cageVerteciesAfterMap,1);

deltaS=circshift(deltaS,size(deltaS,1)-1);
deltaD=circshift(deltaD,size(deltaD,1)-1);

gz=0.5*((deltaD.*(abs(deltaS)+abs(deltaD)))./(deltaS.*abs(deltaD)));
gz_gag=(-0.5)*((deltaD.*(abs(deltaS)-abs(deltaD)))./(deltaS.*abs(deltaD)));

gz_enc=repelem(gz,NumOfVerticesInEdges);
gz_gag_enc=repelem(gz_gag,NumOfVerticesInEdges);

%*************step 3,4:extract argument from gz, and evaluate log(gz), Vg on A******************
%ln_gz=log(abs(gz_enc));
log_gz=0.000000000001*1i+zeros(size(gz_enc));%get code from weber
Vg=(conj(gz_gag_enc))./gz_enc;

%*************step 5:solve 22 - obtain l(z), V(z)******************
cvx_begin
    variable  l(n) complex;
    variable v(n) complex;
    minimize norm(C_sizeA*l-log_gz,2)+lambda*norm(C_sizeA*v-Vg,2);
    subject to
        abs(C_sizeA*v)<=k;
        abs(C_sizeA*v)<=log(SIGMA)-real(C_sizeA*l);
        sigma*exp(-real(C_sizeA*l))+abs(C_sizeA*v)<=1;
cvx_end

Vz=C_sizeM*l;
lz=C_sizeM*v;

%*************step 6:find phi(z) - integral******************
PHItag=exp(lz);

if(exist('treeCumSum', 'file') ~= 3)
    fullPathName = which('treeCumSum.cpp');
    [folderName, fileName, ext] = fileparts(fullPathName);
    mex(fullPathName, '-outdir', folderName);
end

PHI_Z0=Z0+mean(cageVerteciesAfterMap)-mean(cageVerteciesB4Map);
%PHI_Z0=Z0;

%calc the integral on the edges
partialCalc_gpu=gpuArray(PHItag(endIndices)) + gpuArray(PHItag(startIndices));
integral_on_edges_gpu=partialCalc_gpu.*edgeVectors_gpu;
integral_on_edges=gather(integral_on_edges_gpu);

% find the integral on all the spanning tree
PHI_Z0=PHI_Z0+0.000000000000000001i;
PHI = treeCumSum(uint32(Z0index), PHI_Z0, integral_on_edges, startIndices, endIndices);

%*************step 7:find PSI - another integral******************
PSItag=Vz.*PHItag;%check the .*

%PSI_Z0=0; check if need to be 0 or Z0
PSI_Z0=Z0;

partialCalc_gpu=gpuArray(PSItag(endIndices)) + gpuArray(PSItag(startIndices));
integral_on_edges_gpu=partialCalc_gpu.*edgeVectors_gpu;
integral_on_edges=gather(integral_on_edges_gpu);

PSI_Z0=PSI_Z0+0.000000000000000001i;
PSI=treeCumSum(uint32(Z0index), PSI_Z0, integral_on_edges, startIndices, endIndices);

%*************step 8:find f******************
f=PHI+conj(PSI);