%********************check if mex files exists********************
make_mex;
%***********************start algorithm****************************
a=size(C_sizeA,1);
n=size(C_sizeA,2);

%*************step 2:Evaluate gz and gz_gag******************
[gz_enc,gz_gag_enc]=EvalGzAndGzBar( cageVerteciesB4Map,cageVerteciesAfterMap,NumOfVerticesInEdgesSizeA );

%*************step 3,4:extract argument from gz, and evaluate log(gz), Vg on A******************

init_l_gz=logarithmExtraction(cageVerteciesB4Map_sizeA, gz_enc, cageVerteciesAfterMap, NumOfVerticesInEdgesSizeA);
init_Vg=(conj(gz_gag_enc))./gz_enc;

%*************step 5:solve 22 - obtain l(z), V(z)******************
% [m,intersectionX]=findLineApproxForCurve(sigma,SIGMA,k);  %temp - create line instead of third condition

l = p_inv*init_l_gz;
v=p_inv*init_Vg;
l_gz_first_step=C_sizeA*l;
Vg_first_step=C_sizeA*v;

l_gz=l_gz_first_step;
Vg=Vg_first_step;

l_gz_gpu=gpuArray(l_gz);
Vg_gpu=gpuArray(Vg);
SIGMA_gpu=gpuArray(SIGMA);
sigma_gpu=gpuArray(sigma);
k_gpu=gpuArray(k);
intersectionX_gpu=gpuArray(intersectionX);
p_inv_gpu=gpuArray(p_inv);
C_sizeA_gpu=gpuArray(C_sizeA);
epsilon_gpu=gpuArray(epsilon);
tic
for iter=1:max_iterations

%     if(all(abs(Vg)<=k+epsilon)) && (all(abs(Vg)<=log(SIGMA)-real(l_gz)+epsilon)) && (all(log(sigma)+m*abs(C_sizeA*v)<=real(C_sizeA*l)+epsilon))
%         break;
%     end
    
    %local
    
    [ abs_Vg,R_l_gz ] = projectToPoly( abs(Vg_gpu),real(l_gz_gpu),SIGMA_gpu,sigma_gpu,k_gpu,intersectionX_gpu);
    
    l_gz_local=complex(R_l_gz, imag(l_gz_gpu));
    Vg_local=abs_Vg.*exp(1i*angle(Vg_gpu));

    x_prevGlobal=[l_gz_gpu;Vg_gpu];
    x_local=[l_gz_local;Vg_local];
    n0=x_prevGlobal-x_local;  
    if(norm(n0)<epsilon_gpu)%stop condition from the article
        break;
    end

	
    %****
    %global
    l = p_inv_gpu*l_gz_local;
    v=p_inv_gpu*Vg_local;
    
    l_gz_gpu=C_sizeA_gpu*l;
    Vg_gpu=C_sizeA_gpu*v;
    
    
end
toc

fprintf('iterations = %d\n\n',iter);

v=gather(v);
l=gather(l);
Vz=C_sizeM*v;
lz=C_sizeM*l;

%*************step 6:find phi(z) - integral******************
PHItag=exp(lz);

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