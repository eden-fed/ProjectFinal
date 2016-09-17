% clc
tic
a=size(C_sizeA,1);
n=size(C_sizeA,2);

%*************step 2:Evaluate gz and gz_gag******************
deltaS=cageVerteciesB4Map([2:end 1])-cageVerteciesB4Map;
deltaD=cageVerteciesAfterMap([2:end 1])-cageVerteciesAfterMap;

gz = 0.5*(abs(deltaD) + abs(deltaS)) .* deltaD ./ (abs(deltaD).*deltaS); %affine transformation with unit normal
gz_gag = 0.5*(abs(deltaD) - abs(deltaS)) .* deltaD ./ (abs(deltaD).*conj(deltaS)); %affine transformation with unit normal

gz_enc=repelem_ours(gz,NumOfVerticesInEdgesSizeA);
gz_gag_enc=repelem_ours(gz_gag,NumOfVerticesInEdgesSizeA);

%*************step 3,4:extract argument from gz, and evaluate log(gz), Vg on A******************

l_gz=logarithmExtraction(cageVerteciesB4Map_sizeA, gz_enc, cageVerteciesAfterMap, NumOfVerticesInEdgesSizeA);
Vg=(conj(gz_gag_enc))./gz_enc;

%*************step 5:solve 22 - obtain l(z), V(z)******************

for ii=1:max_iterations
    
    %global
    l = p_inv*l_gz;
    v=p_inv*Vg;
    
    l_gz=C_sizeA*l;
    Vg=C_sizeA*v;
    
    if(~any(abs(Vg)>k+epsilon)) && (~any(abs(Vg)>log(SIGMA)-real(l_gz)+epsilon)) && (~any(sigma*exp(-real(l_gz))+abs(Vg)>1+epsilon))
        break;
    end
    
    %local
    %try using segments
    if sigma==1 && SIGMA==1
        R_l_gz=zeros(size(l_gz));
        abs_Vg=complex(zeros(size(Vg)));
    else
        crossPoint_vAxis1=-(lambertw((sigma/SIGMA)*exp(1))-1);
        if(k>=crossPoint_vAxis1)% the case that only the first and third equetions hold   
       %    crossPoint_vAxis2=-crossPoint_vAxis1;
            vAxis_segments=linspace(0,crossPoint_vAxis1,201); 
        else %the case that all equetions hold
            vAxis_segments=linspace(0,k,201);
        end
        lAxis_segments=-log((1-vAxis_segments)/sigma);

        A=(lAxis_segments(1:end-1)-lAxis_segments(2:end))./(vAxis_segments(1:end-1)-vAxis_segments(2:end));
        B=lAxis_segments(1:end-1)-A.*vAxis_segments(1:end-1);
        A_length=length(A);
        A=A';B=B';

   %     [abs_Vg,R_l_gz]=localStep(A,B,abs(Vg),real(l_gz),k,log(SIGMA));
        R_l_gz=zeros(size(l_gz));
        abs_Vg=zeros(size(Vg));
        yy1=zeros(size(l_gz));
        xx1=zeros(size(Vg));
        for jj=1:length(l_gz)
        cvx_begin
            variable y1;
            variable x1;
            minimize norm(real(l_gz(jj))-y1,2)+norm(abs(Vg(jj))-x1,2);
            subject to
                x1>=0;
                x1<=k;
                x1<=log(SIGMA)-y1;
                sigma*exp(-y1)+x1<=1;
        cvx_end
         xx1(jj)=x1;
         yy1(jj)=y1;
         [x1 y1]
         cvx_begin quiet
            variable y2;
            variable x2;
            minimize norm(real(l_gz(jj))-y2,2)+norm(abs(Vg(jj))-x2,2);
            subject to
                x2>=0;
                x2<=k;
                x2<=log(SIGMA)-y2;
                y2>=A*x2+B;
         cvx_end
         R_l_gz(jj)=y2;
         abs_Vg(jj)=x2;
         [x2 y2]
        end
    end
    l_gz=complex(R_l_gz, imag(l_gz));
    Vg=abs_Vg.*exp(1i*angle(Vg));

end

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
PHI_Z0=PHI_Z0+1e-10*1i;
PHI = treeCumSum(uint32(Z0index), PHI_Z0, integral_on_edges, startIndices, endIndices);

%*************step 7:find PSI - another integral******************
PSItag=Vz.*PHItag;

PSI_Z0=0;

partialCalc=PSItag(endIndices) + PSItag(startIndices);
integral_on_edges=partialCalc.*edgeVectors;
% integral_on_edges=gather(integral_on_edges_gpu);

PSI_Z0=complex(PSI_Z0);
PSI=treeCumSum(uint32(Z0index), PSI_Z0, integral_on_edges, startIndices, endIndices);

%*************step 8:find f******************
f=PHI+conj(PSI);
toc