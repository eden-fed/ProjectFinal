% clc
a=size(C_sizeA,1);
n=size(C_sizeA,2);
% max_iterations=1000;
% epsilon=0.0001;
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
%*************step 5:obtain l(z)******************
gamma1=log(SIGMA);
gamma2=log(sigma);
%p_inv=pinv(C_sizeA);

for ii=1:max_iterations
    l = p_inv*l_gz;
    l_gz=C_sizeA*l;
    
    if(~any(real(l_gz)>(gamma1+epsilon))) && (~any(real(l_gz)<(gamma2-epsilon)))
        break;
    end
    
    x = max(min(real(l_gz), gamma1), gamma2);
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

partialCalc=PHItag(endIndices) + PHItag(startIndices);
integral_on_edges=partialCalc.*edgeVectors;

% find the integral on all the spanning tree
PHI_Z0=PHI_Z0+1e-10*1i;
PHI = treeCumSum(uint32(Z0index), PHI_Z0, integral_on_edges, startIndices, endIndices);

%*************step 4:find f******************
f=PHI;