%calculate the SDF function at any point in 3D space

sdf_calc=tic;

xRes = linspace(minX,maxX,grid); yRes = linspace(minY,maxY,grid); zRes = linspace(minZ,maxZ,grid);
[x,y,z] = meshgrid(xRes,yRes,zRes); % x y and z are gridXgridXgrid
sdf = zeros(size(x));%sdf is gridXgridXgrid

for ii=1:n				
    sdf=sdf + w(ii)*((x-verticesCrd_n(ii,1)).^2 +(y-verticesCrd_n(ii,2)).^2 +(z-verticesCrd_n(ii,3) ).^2).^1.5 + ...
     + w(ii+n)*((x-os_n(ii,1)).^2+(y-os_n(ii,2)).^2+(z-os_n(ii,3)).^2).^1.5;								
end

not_gpu_time=toc(sdf_calc)
