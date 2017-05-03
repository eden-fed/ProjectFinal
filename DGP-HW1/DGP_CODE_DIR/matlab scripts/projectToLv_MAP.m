%********************check if mex files exists********************
% make_mex;
%***********************start algorithm****************************
a=size(C_sizeA,1);
n=size(C_sizeA,2);

%*************step 2:Evaluate gz and gz_gag******************
[gz_enc,gz_gag_enc]=EvalGzAndGzBar( cageVerteciesB4Map,cageVerteciesAfterMap,NumOfVerticesInEdgesSizeA );

%*************step 3,4:extract argument from gz, and evaluate log(gz), Vg on A******************

init_l_gz=logarithmExtraction(cageVerteciesB4Map_sizeA, gz_enc, cageVerteciesAfterMap, NumOfVerticesInEdgesSizeA);
init_Vg=(conj(gz_gag_enc))./gz_enc;

% *************step 5:solve 22 - obtain l(z), V(z)******************

l = p_inv*init_l_gz;
v=p_inv*init_Vg;
l_gz_first_step=C_sizeA*l;
Vg_first_step=C_sizeA*v;

l_gz=l_gz_first_step;
Vg=Vg_first_step;

% figure('position', [610, 0, 600, 1400])
% h=plot(0,0);

for iter=1:max_iterations
    
    %graph:
%     energy=sum_square_abs(l_gz_first_step-l_gz)+sum_square_abs(Vg_first_step-Vg);
%     delete(h);
%     h=convex_graph_with_map( sigma, SIGMA, k, l_gz, Vg , energy, m, log(sigma));
    
    %local step:
    prev_abs_Vg=abs(Vg);
    [ abs_Vg,R_l_gz ] = projectToPoly( prev_abs_Vg,real(l_gz),SIGMA,sigma,k,intersectionX);
    
    l_gz_local=complex(R_l_gz, imag(l_gz));
    Vg_local=abs_Vg.*Vg./prev_abs_Vg;
    
    %stop condition check:
    x_prevGlobal=[l_gz;Vg];    
    x_local=[l_gz_local;Vg_local];
        
    n0=x_prevGlobal-x_local;
    if(norm(n0)<epsilon)%stop condition
        break;
    end


    %global step:
    l = p_inv*l_gz_local;
    v=p_inv*Vg_local;
    
    l_gz=C_sizeA*l;
    Vg=C_sizeA*v;

end

% fprintf('max k = %d\n',max(abs(Vg)));
% fprintf('max SIGMA = %d\n',max(exp(real(l_gz)).*(1+abs(Vg))));
% fprintf('min sigma = %d\n',min(exp(real(l_gz)).*(1-abs(Vg))));
% energy=sum_square_abs(l_gz_first_step-l_gz)+sum_square_abs(Vg_first_step-Vg);
% fprintf('energy = %d\n',energy);
% fprintf('iterations = %d\n\n',iter);

% figure('position', [610, 0, 600, 1400])
% convex_graph_with_map( sigma, SIGMA, k, l_gz, Vg, energy , m, log(sigma));

Vz=C_sizeM*v;
lz=C_sizeM*l;

%*************step 6:find phi(z) - integral******************
PHItag=exp(lz);

Cz0=C_sizeM(Z0index,:);
cageAfterMapSizeN=EmcCageVerteciesEdgeWise( cageVerteciesAfterMap, NumOfVerticesInEdgesSizeNlarge, n );
PHI_Z0=Cz0*cageAfterMapSizeN;

%calc the integral on the edges
partialCalc=PHItag(endIndices) + PHItag(startIndices);
integral_on_edges=partialCalc.*edgeVectors;

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

if (isreal(PSI_Z0))
    PSI_Z0=complex(PSI_Z0);
end
if (isreal(integral_on_edges))
    integral_on_edges=complex(integral_on_edges);
end
PSI=treeCumSum(uint32(Z0index), PSI_Z0, integral_on_edges, startIndices, endIndices);

%*************step 8:find f******************
f=PHI+conj(PSI);