% clc
a=size(C_sizeA,1);
n=size(C_sizeA,2);

%*************step 2:Evaluate gz and gz_gag******************
[gz_enc,gz_gag_enc]=EvalGzAndGzBar( cageVerteciesB4Map,cageVerteciesAfterMap,NumOfVerticesInEdgesSizeA );

%*************step 3,4:extract argument from gz, and evaluate log(gz), Vg on A******************

init_l_gz=logarithmExtraction(cageVerteciesB4Map_sizeA, gz_enc, cageVerteciesAfterMap, NumOfVerticesInEdgesSizeA);
init_Vg=(conj(gz_gag_enc))./gz_enc;
%        l_gz=init_l_gz;Vg=init_Vg;
l = p_inv*init_l_gz;%first make l_gz and Vg a map ***temporary***
v=p_inv*init_Vg;
l_gz=C_sizeA*l;
Vg=C_sizeA*v;
%*************step 5:solve 22 - obtain l(z), V(z)******************
cvx_begin quiet
    variable  l(n) complex;
    variable v(n) complex;
 %   minimize norm(C_sizeA*l-l_gz,2)+lambda*norm(C_sizeA*v-Vg,2);
    minimize sum_square_abs(C_sizeA*l-l_gz)+lambda*sum_square_abs(C_sizeA*v-Vg);
    subject to
        abs(C_sizeA*v)<=k;
        abs(C_sizeA*v)<=log(SIGMA)-real(C_sizeA*l);
        sigma*exp(-real(C_sizeA*l))+abs(C_sizeA*v)<=1;
cvx_end

l_gz=C_sizeA*l;
Vg=C_sizeA*v;
fprintf('max k = %d\n',max(abs(Vg)));
fprintf('max SIGMA = %d\n',max(exp(real(l_gz)).*(1+abs(Vg))));
fprintf('min sigma = %d\n',min(exp(real(l_gz)).*(1-abs(Vg))));
fprintf('energy = %d\n\n',norm(l_gz-init_l_gz,2)+lambda*norm(Vg-init_Vg,2));

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
if (isreal(PHI_Z0))
PHI_Z0=complex(PHI_Z0);
end
PHI = treeCumSum(uint32(Z0index), PHI_Z0, integral_on_edges, startIndices, endIndices);

%*************step 7:find PSI - another integral******************
PSItag=Vz.*PHItag;

PSI_Z0=0;

partialCalc=PSItag(endIndices) + PSItag(startIndices);
integral_on_edges=partialCalc.*edgeVectors;
% integral_on_edges=gather(integral_on_edges_gpu);

if (isreal(PSI_Z0))
PSI_Z0=complex(PSI_Z0);
end
if (isreal(integral_on_edges))
integral_on_edges=complex(integral_on_edges);
end
PSI=treeCumSum(uint32(Z0index), PSI_Z0, integral_on_edges, startIndices, endIndices);

%*************step 8:find f******************
f=PHI+conj(PSI);