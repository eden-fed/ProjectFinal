function [ A,B ] = createLineSegments( sigma,SIGMA,k,numOfSegments)
 crossPoint_vAxis1=-(lambertw((sigma/SIGMA)*exp(1))-1);
if(k>=crossPoint_vAxis1)% the case that only the first and third equetions hold   
    vAxis_segments=linspace(0,crossPoint_vAxis1,numOfSegments); 
else %the case that all equetions hold
    vAxis_segments=linspace(0,k,numOfSegments);
end
lAxis_segments=-log((1-vAxis_segments)/sigma);

A=(lAxis_segments(1:end-1)-lAxis_segments(2:end))./(vAxis_segments(1:end-1)-vAxis_segments(2:end));
B=lAxis_segments(1:end-1)-A.*vAxis_segments(1:end-1);
A=A';B=B';
end

