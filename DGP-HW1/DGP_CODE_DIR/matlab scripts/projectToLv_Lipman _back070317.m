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
[m,intersectionX]=findLineForTest(sigma,SIGMA,k);%temp - create line instead of third condition

l = p_inv*init_l_gz;%for the case of 1 iteration
v=p_inv*init_Vg;
l_gz_first_step=C_sizeA*l;
Vg_first_step=C_sizeA*v;

l_gz=l_gz_first_step;
Vg=Vg_first_step;

%******preprocess******
T=blkdiag(C_sizeA,C_sizeA);
T_trans=T';
% M=sparse(T_trans*T);
% [M_L, M_U, M_p, M_q] = lu(M, 'vector');
%**********************

% figure('position', [1220, 0, 600, 1400])
% h=plot(0,0);
for iter=1:max_iterations
    
    %local
    [ abs_Vg,R_l_gz ] = projectToPoly( abs(Vg),real(l_gz),SIGMA,sigma,k,intersectionX);
    
    l_gz_local=complex(R_l_gz, imag(l_gz));
    Vg_local=abs_Vg.*exp(1i*angle(Vg));
    %****
    
    %global
%     energy=sum_square_abs(l_gz_first_step-l_gz)+sum_square_abs(Vg_first_step-Vg);
%     delete(h);
%     h=convex_graph_with_map( sigma, SIGMA, k, l_gz, Vg , energy, m, log(sigma));
%     
    n0=[l_gz;Vg]-[l_gz_local;Vg_local];
    if(norm(n0)<epsilon)%stop condition from the article
        break;
    end

    KKTmat=[T'*T , T'*n0 ; n0'*T , 0];
    KKTresult=KKTmat \ [T'*[l_gz;Vg] ; n0'*[l_gz_local;Vg_local]];
    %fprintf('cond(KKTmat) = %d\n',cond(KKTmat));

%     eta_0=T_trans*n0;
%     d_0=n0'*[l_gz_local;Vg_local];
%     c_0=T_trans*[l_gz;Vg];
% 
%     y_c = M_U\(M_L\(c_0(M_p,:)));
%     y_eta = M_U\(M_L\(eta_0(M_p,:)));
%     
%     KKTresult=y_c-(y_eta*(eta_0'*y_c-d_0)/(eta_0'*y_eta));

%     if(any(abs(KKTresult_old(1:end-1)-KKTresult)>1e-6))
%        warning('different result'); 
%     end
    l=KKTresult(1:n);
    v=KKTresult(n+1:2*n);%(n+1:end)
    %***
    l_gz=C_sizeA*l;
    Vg=C_sizeA*v;
    
end
if(all(abs(Vg)<=k+epsilon)) && (all(abs(Vg)<=log(SIGMA)-real(l_gz)+epsilon)) && (all(m*abs(C_sizeA*v)+log(sigma)<=real(C_sizeA*l)+epsilon))
    fprintf('the constraints are satisfied\n');
end

fprintf('max k = %d\n',max(abs(Vg)));
fprintf('max SIGMA = %d\n',max(exp(real(l_gz)).*(1+abs(Vg))));
fprintf('min sigma = %d\n',min(exp(real(l_gz)).*(1-abs(Vg))));
fprintf('energy = %d\n',sum_square_abs(l_gz_first_step-l_gz)+sum_square_abs(Vg_first_step-Vg));
fprintf('iterations = %d\n\n',iter);

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