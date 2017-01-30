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
[m,b]=findLineForTest(sigma,SIGMA,k);  %temp - create line instead of third condition

l = p_inv*init_l_gz;
v=p_inv*init_Vg;
l_gz_first_step=C_sizeA*l;
Vg_first_step=C_sizeA*v;

l_gz=l_gz_first_step;
Vg=Vg_first_step;
figure('position', [610, 0, 600, 1400])
h=plot(0,0);
for iter=1:max_iterations
    
    energy=sum_square_abs(l_gz_first_step-l_gz)+sum_square_abs(Vg_first_step-Vg);
    delete(h);
    h=convex_graph_with_map( sigma, SIGMA, k, l_gz, Vg , energy, m, b);
    if(all(abs(Vg)<=k+epsilon)) && (all(abs(Vg)<=log(SIGMA)-real(l_gz)+epsilon)) && (all(m*abs(C_sizeA*v)+b<=real(C_sizeA*l)+epsilon))
        break;
    end
    
    %local
%     cvx_begin quiet
%         variable R_l_gz(a);
%         variable abs_Vg(a);
%         minimize sum_square_abs(real(l_gz)-R_l_gz)+sum_square_abs(abs(Vg)-abs_Vg);
%         subject to
%             abs_Vg>=0;
%             abs_Vg<=k;
%             abs_Vg<=log(SIGMA)-R_l_gz;
%             m*abs_Vg+b<=R_l_gz;
%             %sigma*exp(-R_l_gz)+abs_Vg<=1;
%     cvx_end
%     
%     l_gz=complex(R_l_gz, imag(l_gz));
%     Vg=abs_Vg.*exp(1i*angle(Vg));

    cvx_begin quiet
        variable l_gz_local(a) complex;
        variable Vg_local(a) complex;
        E1=sum_square_abs(l_gz-l_gz_local)+sum_square_abs(Vg-Vg_local);
        minimize E1;
        subject to
            abs(Vg_local)<=k;
            abs(Vg_local)<=log(SIGMA)-real(l_gz_local);
            m*abs(Vg_local)+b<=real(l_gz_local);
    cvx_end
    
    l_gz=l_gz_local;
    Vg=Vg_local;
    %****
    
    %global
    l = p_inv*l_gz;
    v=p_inv*Vg;
    
    l_gz=C_sizeA*l;
    Vg=C_sizeA*v;
    
    
    %     if(energy>=407.61)
    %         break;
    %     end
    %     if(all(abs(Vg)<=k+epsilon)) && (all(abs(Vg)<=log(SIGMA)-real(l_gz)+epsilon)) && (all(sigma*exp(-real(l_gz))+abs(Vg)<=1+epsilon))
    %         break;
    %     end
    
    
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