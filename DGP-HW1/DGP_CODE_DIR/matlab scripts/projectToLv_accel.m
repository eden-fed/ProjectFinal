% clc
a=size(C_sizeA,1);
n=size(C_sizeA,2);

%*************step 2:Evaluate gz and gz_gag******************
deltaS=cageVerteciesB4Map-circshift(cageVerteciesB4Map,1);
deltaD=cageVerteciesAfterMap-circshift(cageVerteciesAfterMap,1);

deltaS=circshift(deltaS,size(deltaS,1)-1);
deltaD=circshift(deltaD,size(deltaD,1)-1);

gz = 0.5*(abs(deltaD) + abs(deltaS)) .* deltaD ./ (abs(deltaD).*deltaS); %affine transformation with unit normal
gz_gag = 0.5*(abs(deltaD) - abs(deltaS)) .* deltaD ./ (abs(deltaD).*conj(deltaS)); %affine transformation with unit normal

gz_enc=repelem_ours(gz,NumOfVerticesInEdgesSizeA);
gz_gag_enc=repelem_ours(gz_gag,NumOfVerticesInEdgesSizeA);

%*************step 3,4:extract argument from gz, and evaluate log(gz), Vg on A******************

l_gz=logarithmExtraction(cageVerteciesB4Map_sizeA, gz_enc, cageVerteciesAfterMap, NumOfVerticesInEdgesSizeA);
Vg=(conj(gz_gag_enc))./gz_enc;

%*************step 5:solve 22 - obtain l(z), V(z)******************

for ii=1:max_iterations
    
    l = p_inv*l_gz;
    v=p_inv*Vg;
    l_gz=C_sizeA*l;
    
    cvx_begin
        variable R_l_gz(a);
        variable Vg(a) complex;
        minimize norm(real(l_gz)-R_l_gz,2)+norm(C_sizeA*v-Vg,2);
        subject to
            abs(Vg)<=k;
            abs(Vg)<=log(SIGMA)-R_l_gz;
            sigma*exp(-R_l_gz)+abs(Vg)<=1;
    cvx_end
    
    l_gz=complex(R_l_gz, imag(l_gz));

end

%check the distortion:
sum(abs(C_sizeA*v)>k)
sum(abs(C_sizeA*v)>log(SIGMA)-real(C_sizeA*l))
sum(sigma*exp(-real(C_sizeA*l))+abs(C_sizeA*v)>1)
%if all three are 0 then there is no distortion

Vz=C_sizeM*v;
lz=C_sizeM*l;

%*************step 6:find phi(z) - integral******************
PHItag=exp(lz);

if(exist('treeCumSum', 'file') ~= 3)
    fullPathName = which('treeCumSum.cpp');
    [folderName, fileName, ext] = fileparts(fullPathName);
    mex(fullPathName, '-outdir', folderName);
end

Cz0=C_sizeM(Z0index,:);
cageAfterMapSizeN=EmcCageVerteciesEdgeWise( cageVerteciesAfterMap, NumOfVerticesInEdgesSizeNlarge );
PHI_Z0=Cz0*cageAfterMapSizeN;

%calc the integral on the edges
partialCalc=PHItag(endIndices) + PHItag(startIndices);
integral_on_edges=partialCalc.*edgeVectors;
% integral_on_edges=gather(integral_on_edges_gpu);

% find the integral on all the spanning tree
PHI_Z0=PHI_Z0+0.000000000000000001i;
PHI = treeCumSum(uint32(Z0index), PHI_Z0, integral_on_edges, startIndices, endIndices);

%*************step 7:find PSI - another integral******************
PSItag=Vz.*PHItag;

PSI_Z0=0;

partialCalc=PSItag(endIndices) + PSItag(startIndices);
integral_on_edges=partialCalc.*edgeVectors;
% integral_on_edges=gather(integral_on_edges_gpu);

PSI_Z0=PSI_Z0+0.000000000000000001i;
PSI=treeCumSum(uint32(Z0index), PSI_Z0, integral_on_edges, startIndices, endIndices);

%*************step 8:find f******************
f=PHI+conj(PSI);