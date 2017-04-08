%********************check if mex files exists********************
% make_mex;
%***********************start algorithm****************************
% timerLG=tic;
a=size(C_sizeA,1);
n=size(C_sizeA,2);

%*************step 2:Evaluate gz and gz_gag******************
[gz_enc,gz_gag_enc]=EvalGzAndGzBar( cageVerteciesB4Map,cageVerteciesAfterMap,NumOfVerticesInEdgesSizeA );

%*************step 3,4:extract argument from gz, and evaluate log(gz), Vg on A******************

init_l_gz=logarithmExtraction(cageVerteciesB4Map_sizeA, gz_enc, cageVerteciesAfterMap, NumOfVerticesInEdgesSizeA);

init_Vg=(conj(gz_gag_enc))./gz_enc;

% *************step 5:solve 22 - obtain l(z), V(z)******************
% [m,intersectionX]=findLineApproxForCurve(sigma,SIGMA,k);  %temp - create line instead of third condition

l = p_inv*init_l_gz;
v=p_inv*init_Vg;
l_gz_first_step=C_sizeA*l;
Vg_first_step=C_sizeA*v;

l_gz=l_gz_first_step;
Vg=Vg_first_step;
% 
% energy_graph=zeros(max_iterations,1);%*
% temp_idx=1;
% 
% energy_graph2=zeros(2*max_iterations,1);%*
% temp_idx2=1;
% % 
% figure('position', [610, 0, 600, 1400])
% h=plot(0,0);

% converge=ones(max_iterations,1);
% energy_=ones(max_iterations,1);

for iter=1:max_iterations
    
%     energy=sum_square_abs(l_gz_first_step-l_gz)+sum_square_abs(Vg_first_step-Vg);
%     delete(h);
%     h=convex_graph_with_map( sigma, SIGMA, k, l_gz, Vg , energy, m, log(sigma));
% 	energy_graph(temp_idx)=energy;%*
%     temp_idx=temp_idx+1;
% 	x_prevGlobal=[l_gz;Vg];
%     if(all(abs(Vg)<=k+epsilon)) && (all(abs(Vg)<=log(SIGMA)-real(l_gz)+epsilon)) && (all(log(sigma)+m*abs(Vg)<=real(l_gz)+epsilon))
% %         fprintf('the constraints are satisfied\n');
%         break;
%     end
    
    %local

    prev_abs_Vg=abs(Vg);
    [ abs_Vg,R_l_gz ] = projectToPoly( prev_abs_Vg,real(l_gz),SIGMA,sigma,k,intersectionX);
    
    l_gz_local=complex(R_l_gz, imag(l_gz));
%     Vg_local=abs_Vg.*exp(1i*angle(Vg));  
    Vg_local=abs_Vg.*Vg./prev_abs_Vg;
    
    x_prevGlobal=[l_gz;Vg];    
    x_local=[l_gz_local;Vg_local];
        
    n0=x_prevGlobal-x_local;
    if(norm(n0)<epsilon)%stop condition from the article
        break;
    end
    
% 	E1=sum_square_abs(C_sizeA*l-l_gz)+sum_square_abs(C_sizeA*v-Vg);
%  	energy_graph2(temp_idx2)=E1;%*
%     temp_idx2=temp_idx2+1;

%****
    %global
    l = p_inv*l_gz_local;
    v=p_inv*Vg_local;
%     E1=sum_square_abs(C_sizeA*l-l_gz)+sum_square_abs(C_sizeA*v-Vg);
% 	energy_graph2(temp_idx2)=E1;%*
%     temp_idx2=temp_idx2+1;
    
    l_gz=C_sizeA*l;
    Vg=C_sizeA*v;
    %     if(all(abs(Vg)<=k+epsilon)) && (all(abs(Vg)<=log(SIGMA)-real(l_gz)+epsilon)) && (all(sigma*exp(-real(l_gz))+abs(Vg)<=1+epsilon))
    %         break;
    %     end
%     x=[l_gz;Vg];
%     converge(iter)=norm(x-x_infty)/(norm(x_prevGlobal-x_infty));
end

% fprintf('max k = %d\n',max(abs(Vg)));
% fprintf('max SIGMA = %d\n',max(exp(real(l_gz)).*(1+abs(Vg))));
% fprintf('min sigma = %d\n',min(exp(real(l_gz)).*(1-abs(Vg))));
% energy=sum_square_abs(l_gz_first_step-l_gz)+sum_square_abs(Vg_first_step-Vg);
% fprintf('energy = %d\n',energy);
% fprintf('iterations = %d\n\n',iter);

% figure('position', [610, 0, 600, 1400])
% convex_graph_with_map( sigma, SIGMA, k, l_gz, Vg, energy , m, log(sigma));

% x=1:iter;%*
% energy_graph=energy_graph(x);%*
% figure%*
% plot(x,energy_graph,'LineWidth',3);%*
% title('energy = sumSquareAbs(l_gzFirstStep-l_gz)+sumSquareSbs(V_gFirstStep-V_g)');
% xlabel('iterations');%*
% ylabel('energy');%*
% 
% x=1:2*iter;%*
% energy_graph2=energy_graph2(x);%*
% figure%*
% plot(x,energy_graph2,'LineWidth',3);%*
% title(['energy = ','E1+delta*E2']);
% ylabel('energy');%*

Vz=C_sizeM*v;
lz=C_sizeM*l;

%*************step 6:find phi(z) - integral******************
PHItag=exp(lz);

Cz0=C_sizeM(Z0index,:);
cageAfterMapSizeN=EmcCageVerteciesEdgeWise( cageVerteciesAfterMap, NumOfVerticesInEdgesSizeNlarge, n );
PHI_Z0=Cz0*cageAfterMapSizeN;
% PHI_Z0=Cz0*cageVerteciesAfterMapSizeNLarge;

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
% timerLGtime=toc(timerLG)