%calculate the SDF function at any point in 3D space

sdf_calc=tic;

xRes = gpuArray.linspace(minX,maxX,grid); yRes = gpuArray.linspace(minY,maxY,grid); zRes = gpuArray.linspace(minZ,maxZ,grid);
[x,y,z] = meshgrid(xRes,yRes,zRes); % x y and z are gridXgridXgrid
sdf = zeros(size(x),'gpuArray');%sdf is gridXgridXgrid

verticesCrd_gpu=gpuArray(verticesCrd_n);
os_gpu=gpuArray(os_n);
w_gpu=gpuArray(w);
n_gpu=gpuArray(n);

for ii=1:n				
    sdf=sdf + w_gpu(ii)*((x-verticesCrd_gpu(ii,1)).^2 +(y-verticesCrd_gpu(ii,2)).^2 +(z-verticesCrd_gpu(ii,3) ).^2).^1.5 + ...
     + w_gpu(ii+n_gpu)*((x-os_gpu(ii,1)).^2+(y-os_gpu(ii,2)).^2+(z-os_gpu(ii,3)).^2).^1.5;								
end

sdf=gather(sdf);
x=gather(x);
y=gather(y);
z=gather(z);

gpu_time=toc(sdf_calc)
