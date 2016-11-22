% clc
a=size(Ctag,1);
n=size(Ctag,2);
%k=(SIGMA-sigma)/(SIGMA+sigma);
%lambda=1;
%*************step 2:Evaluate gz and gz_gag******************
[gz_enc,gz_gag_enc]=EvalGzAndGzBar( cageVerteciesB4Map,cageVerteciesAfterMap,NumOfVerticesInEdges );

%*************step 3:solve 13 - obtain r, psi(z)******************
Cz0=C_sizeM(Z0index,:);

%R=abs(Ctag_tempAN*cageVerteciesAfterMap);
abs_gz=abs(gz_enc);
conj_gs_gag=conj(gz_gag_enc);

cvx_begin
    variable  psai(n) complex;
    variable r(a);
    minimize norm(r-abs_gz, 2) + lambda*norm(Ctag*psai-conj_gs_gag, 2);
    %minimize norm(r-abs_gz, 2) + lambda*norm(Ctag*psai-conj(conj_gs_gag), 2);
    subject to
		Cz0*psai == zeros(size(Cz0*psai));
        abs(Ctag*psai)<=k*r;
        abs(Ctag*psai)<=SIGMA-r;
        abs(Ctag*psai)<=r-sigma;
cvx_end

%*************step 4:solve 15 - obtain l(z)******************

% normB4map=cageVerteciesB4Map/norm(cageVerteciesB4Map);
% normAfterMap=cageVerteciesAfterMap/norm(cageVerteciesAfterMap);
cvx_begin
    variable  l(n) complex;
    minimize norm(real(C_sizeA*l)-log(r),2);
    subject to
        %imag(Cz0*l)==mean(angle(normB4map-mean(normB4map))-angle(normAfterMap-mean(normAfterMap)));
        %imag(Cz0*l)==2*mean(angle(cageVerteciesB4Map-cageVerteciesAfterMap));
      %  imag(Cz0*l)==2*mean(angle(cageVerteciesB4Map-cageVerteciesAfterMap));
        imag(Cz0*l)==0;
cvx_end
%*************step 5:find derivative of phi(z)******************
PHItag=exp(C_sizeM*l);

%*************step 6:find phi(z) - integral******************
%we need to send :cageVerteciesAfterMap,cageVerteciesB4Map, startIndices, endIndices, integral_on_edges(or edgeVectors_gpu)
if(exist('treeCumSum', 'file') ~= 3)
    fullPathName = which('treeCumSum.cpp');
    [folderName, fileName, ext] = fileparts(fullPathName);
    mex(fullPathName, '-outdir', folderName);
end

%PHI_Z0=Z0+mean(cageVerteciesAfterMap)-mean(cageVerteciesB4Map);
PHI_Z0=Z0;
%calc the integral on the edges
partialCalc=PHItag(endIndices) + PHItag(startIndices);
integral_on_edges=partialCalc.*edgeVectors;
% integral_on_edges=gather(integral_on_edges_gpu);

% find the integral on all the spanning tree
PHI_Z0=PHI_Z0+0.000000000000000001i;
PHI = treeCumSum(uint32(Z0index), PHI_Z0, integral_on_edges, startIndices, endIndices);

%*************step 7:find f******************
PSI=C_sizeM*psai;
f=PHI+conj(PSI);

%######################################################################################################################
%find f for cage in size a
%*************step 5:find derivative of phi(z)******************
PHItag=exp(C_sizeA*l);

%*************step 6:find phi(z) - integral******************
PHI_rootPointOnCage=cageVerteciesB4Map_sizeA(1);

%calc the integral on the edges
startIndices_a=1:(size(cageVerteciesB4Map_sizeA,1)-1);
endIndices_a=2:size(cageVerteciesB4Map_sizeA,1);

integral_on_edges_gpu=(gpuArray(PHItag(endIndices_a)) + gpuArray(PHItag(startIndices_a))).*((gpuArray(cageVerteciesB4Map_sizeA(endIndices_a)) - gpuArray(cageVerteciesB4Map_sizeA(startIndices_a)))./2);
integral_on_edges=gather(integral_on_edges_gpu);

% find the integral on all the spanning tree
PHI_rootPointOnCage=PHI_rootPointOnCage+0.000000000000000001i;
PHI_on_cage = treeCumSum(uint32(1), PHI_rootPointOnCage, integral_on_edges, uint32(startIndices_a), uint32(endIndices_a));

%*************step 7:find f******************
PSI=C_sizeA*psai;
f_on_cage=PHI_on_cage+conj(PSI);
%######################################################################################################################

%************determine rotation and translation**************
%create p_internal - the map f on the internal points
p_internal=zeros(size(f,1),2);%convert the complex vector to matrix
p_internal(:,1)=real(f);
p_internal(:,2)=imag(f);

%create p_on_cage - the map f on the cage (size a)
p_on_cage=zeros(size(f_on_cage,1),2);%convert the complex vector to matrix
p_on_cage(:,1)=real(f_on_cage);
p_on_cage(:,2)=imag(f_on_cage);


sampledCageVerteciesAfterMap=EmcCageVerteciesEdgeWise( cageVerteciesAfterMap, NumOfVerticesInEdges );

%convert the complex vector to matrix
q=zeros(size(sampledCageVerteciesAfterMap,1),2);
q(:,1)=real(sampledCageVerteciesAfterMap);
q(:,2)=imag(sampledCageVerteciesAfterMap);

%calc rotation and translation
[ R,T ] = CalcTranslationAndRotation4H( q , p_on_cage);

%f_cartezian=p_internal*R+repmat(T',size(p_internal,1),1);
f_cartezian=(R*(p_internal'))'+repmat(T',size(p_internal,1),1);
f=f_cartezian(:,1)+1i*f_cartezian(:,2);

