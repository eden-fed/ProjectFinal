
m=size(verticesCrd,1);%number of input points

minX=min(verticesCrd(:,1)); minY=min(verticesCrd(:,2)); minZ=min(verticesCrd(:,3));
maxX=max(verticesCrd(:,1)); maxY=max(verticesCrd(:,2)); maxZ=max(verticesCrd(:,3));

%if epsilon is not set by user, we copute it to be the average distance between points. (divided by 10 for better results)
if(eps==0)    
    dist=pdist2(verticesCrd,verticesCrd);
    eps=(sum(sum(dist))/(m^2-m))/10;

end

%for models containing up to n vertices - we should have 2n equetions 
%i.e. RBF num=vertices num
if(m<n)
    n=m;
end

%if the actual num of input points (m) is larger then n, randomly pick n points
%out of the entire point set
randVerticesInd=sort(randperm(m,n));%the index of the selected vertices

d=[zeros(m,1);(ones(m,1))*eps];%create the vector d (n zeros and then n epsilons)

os=verticesCrd+eps*normalsCrd;%create the off surface points

K=zeros(2*m,2*n);%this will be the matrix K

for ii=1:n
  temp=verticesCrd-repmat(verticesCrd(randVerticesInd(ii),:),[length(verticesCrd) 1]);
  K(1:m,ii)=(sum((temp).^2,2)).^(1.5);%norm(pi-pj)^3
  temp=os-repmat(verticesCrd(randVerticesInd(ii),:),[length(os) 1]);
  K(m+1:2*m,ii)=(sum(temp.^2,2)).^(1.5);%norm((pj+eps*nj)-pi)^3
end

for ii=(n+1):(2*n)
    temp=verticesCrd-repmat(os(randVerticesInd(ii-n),:),[length(verticesCrd) 1]);
    K(1:m,ii)=(sum((temp).^2,2)).^(1.5);%norm(pj-(pi+eps*ni))^3
    temp=os-repmat(os(randVerticesInd(ii-n),:),[length(os) 1]);
    K(m+1:2*m,ii)=(sum(temp.^2,2)).^(1.5);%norm((pj+eps*nj)-pi)^3
end

%compute the least squares solution to the liner system K*w=d
w=K\d;

verticesCrd_n=verticesCrd(randVerticesInd,:);
os_n=os(randVerticesInd,:);
%calculate the SDF function at any point in 3D space
calc_sdf_gpu
%calc_sdf_not_gpu

%extract the zero level set
[faces,vertices] = isosurface(x,y,z,sdf,0);%create faces and vetices from the isosurface function
