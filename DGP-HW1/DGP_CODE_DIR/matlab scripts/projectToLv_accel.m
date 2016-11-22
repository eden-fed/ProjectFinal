%********************check if mex files exists********************
make_mex;
%***********************start algorithm****************************
a=size(C_sizeA,1);
n=size(C_sizeA,2);

%*************step 2:Evaluate gz and gz_gag******************
[gz_enc,gz_gag_enc]=EvalGzAndGzBar( cageVerteciesB4Map,cageVerteciesAfterMap,NumOfVerticesInEdgesSizeA );

%*************step 3,4:extract argument from gz, and evaluate log(gz), Vg on A******************

l_gz=logarithmExtraction(cageVerteciesB4Map_sizeA, gz_enc, cageVerteciesAfterMap, NumOfVerticesInEdgesSizeA);
Vg=(conj(gz_gag_enc))./gz_enc;
first_l_gz=l_gz;first_Vg=Vg;

%*************step 5:solve 22 - obtain l(z), V(z)******************
mag=zeros(max_iterations,1);%*
energy=zeros(max_iterations,1);%*
for ii=1:max_iterations
    
    %local
    if sigma==1 && SIGMA==1
        R_l_gz=zeros(size(l_gz));
        abs_Vg=complex(zeros(size(Vg)));
    else
        [abs_Vg,R_l_gz]=localStep(A,B,abs(Vg),real(l_gz),k,log(SIGMA));
    end
    
    l_gz_local=complex(R_l_gz, imag(l_gz));
    Vg_local=abs_Vg.*exp(1i*angle(Vg));
    
     n0=[l_gz;Vg]-[l_gz_local;Vg_local]; %*
     
     l_gz=l_gz_local;Vg=Vg_local;%*
    %global
    l = p_inv*l_gz;
    v=p_inv*Vg;
    
    l_gz=C_sizeA*l;
    Vg=C_sizeA*v;
      
    if(all(abs(Vg)<k+epsilon)) && (all(abs(Vg)<log(SIGMA)-real(l_gz)+epsilon)) && (all(sigma*exp(-real(l_gz))+abs(Vg)<1+epsilon))
        break;
    end
    
    mag(ii)=norm(n0);%*
    energy(ii)=norm(C_sizeA*l-l_gz_local,2)+norm(C_sizeA*v-Vg_local,2);%*
end
x=1:ii;%*
mag=mag(x);%*
figure
plot(x,mag,'-o');%*
xlabel('iterations');%*
ylabel('norm(n0)');%*

energy=energy(x);
figure
plot(x,energy,'-o');%*
xlabel('iterations');%*
ylabel('energy');%*

fprintf('iterations = %d\n',ii);

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