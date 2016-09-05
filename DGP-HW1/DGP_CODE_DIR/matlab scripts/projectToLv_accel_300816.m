% clc
a=size(C_sizeA,1);
n=size(C_sizeA,2);

%*************step 1:Evaluate gz ******************
sourceEdges=cageVerteciesB4Map_sizeA([2:end 1])-cageVerteciesB4Map_sizeA;

cageVerteciesAfterMap_sizeA=EmcCageVerteciesEdgeWise( cageVerteciesAfterMap, NumOfVerticesInEdgesSizeA );
destEdges=cageVerteciesAfterMap_sizeA([2:end 1])-cageVerteciesAfterMap_sizeA;

% gz=destEdges/sourceEdges;
% gz_gag=0;
%*************step 2:evaluate log(gz), Vg on A******************
sourceEdges_minus1=sourceEdges([end 1:end-1]);
destEdges_minus1=destEdges([end 1:end-1]);

d=log(destEdges./destEdges_minus1)-log(sourceEdges./sourceEdges_minus1);
d(1)=0;
l_gz=log(destEdges(1)./sourceEdges(1)) + cumsum(d); 
BBB=l_gz;
%*************step 5:obtain l(z)******************
gamma1=log(SIGMA);
gamma2=log(sigma);
x_axis=1:size(l_gz,1);

for ii=1:iterations
  
cvx_begin quiet
    variable  l(n) complex;
    minimize norm(C_sizeA*l-l_gz,2);
cvx_end

%disp(cvx_optval);

l_gz=C_sizeA*l;

disp(norm(BBB-l_gz,2));

plot(x_axis, exp(real(l_gz)));
title(num2str(ii));
ylim([0 5]);
drawnow

%l_gz(real(l_gz)>gamma1)=gamma1;
%l_gz(real(l_gz)<gamma2)=gamma2;
x = max(min(real(l_gz), gamma1_gpu), gamma2_gpu);
l_gz = complex(x, imag(l_gz));

end


lz=C_sizeM*l;

%*************step 3:find phi(z) - integral******************
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
% partialCalc_gpu=PHItag(endIndices_gpu) + PHItag(startIndices_gpu);
% integral_on_edges_gpu=partialCalc_gpu.*edgeVectors_gpu;
% integral_on_edges=gather(integral_on_edges_gpu);
% startIndices=gather(startIndices_gpu);
% endIndices=gather(endIndices_gpu);
partialCalc_gpu=gpuArray(PHItag(endIndices)) + gpuArray(PHItag(startIndices));
integral_on_edges_gpu=partialCalc_gpu.*edgeVectors_gpu;
integral_on_edges=gather(integral_on_edges_gpu);

% find the integral on all the spanning tree
PHI_Z0=PHI_Z0+0.000000000000000001i;
PHI = treeCumSum(uint32(Z0index), PHI_Z0, integral_on_edges, startIndices, endIndices);

%*************step 4:find f******************
f=PHI;