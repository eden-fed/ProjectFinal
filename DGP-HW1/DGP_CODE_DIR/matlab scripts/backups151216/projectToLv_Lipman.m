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
mag=zeros(max_iterations,1);%*
energy=zeros(max_iterations,1);%*
    l = p_inv*init_l_gz;%for the case of 1 iteration
    v=p_inv*init_Vg;
    l_gz=C_sizeA*l;
    Vg=C_sizeA*v;
for iter=1:max_iterations
    
    %local
    cvx_begin quiet
        variable R_l_gz(a);
        variable abs_Vg(a);
       % minimize norm(real(l_gz)-R_l_gz,2)+norm(abs(Vg)-abs_Vg,2);
        minimize sum_square_abs(real(l_gz)-R_l_gz)+sum_square_abs(abs(Vg)-abs_Vg);
        subject to
            abs_Vg>=0;
            abs_Vg<=k;
            abs_Vg<=log(SIGMA)-R_l_gz;
            sigma*exp(-R_l_gz)+abs_Vg<=1;
    cvx_end
    %***
    
    l_gz_local=complex(R_l_gz, imag(l_gz));
    Vg_local=abs_Vg.*exp(1i*angle(Vg));
    
    %global
    
    n0=[l_gz;Vg]-[l_gz_local;Vg_local]; 
    mag(iter)=norm(n0);%*
    energy(iter)=norm(l_gz-l_gz_local,2)+norm(Vg-Vg_local,2);%*
    if(norm(n0)<epsilon)%termination condition from the article
        break;
    end
    T=blkdiag(C_sizeA,C_sizeA);
    M=[T'*T , T'*n0 ; n0'*T , 0];
    KKTresult=M \ [T'*[l_gz;Vg] ; n0'*[l_gz_local;Vg_local]];
    %fprintf('cond(M) = %d\n',cond(M));
    
    l=KKTresult(1:n);
    v=KKTresult(n+1:2*n);
    %***
    l_gz=C_sizeA*l;
    Vg=C_sizeA*v;

end
fprintf('max k = %d\n',max(abs(Vg)));
fprintf('max SIGMA = %d\n',max(exp(real(l_gz)).*(1+abs(Vg))));
fprintf('min sigma = %d\n',min(exp(real(l_gz)).*(1-abs(Vg))));
fprintf('energy = %d\n',norm(l_gz-init_l_gz,2)+norm(Vg-init_Vg,2));
fprintf('iterations = %d\n\n',iter);


%graphs:
% x=1:iter;%*
% mag=mag(x);%*
% figure%*
% plot(x,mag,'-o');%*
% xlabel('iterations');%*
% ylabel('norm(n0)');%*
% 
% energy=energy(x);%*
% figure%*
% plot(x,energy,'-o');%*
% xlabel('iterations');%*
% ylabel('energy');%*
    
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