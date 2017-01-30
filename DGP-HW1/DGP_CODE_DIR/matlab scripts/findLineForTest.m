function [ a,b ] = findLineForTest( sigma,SIGMA,k)
    %use the return value as y>=a*x+b.
    point1_vAxis=0;
    point1_lAxis=-log(1/sigma);
    
	point2_vAxis=-(lambertw((sigma/SIGMA)*exp(1))-1);
    if(k<=point2_vAxis)
        point2_vAxis=k;
    end
    point2_lAxis=-log((1-point2_vAxis)/sigma);
    
    a=(point2_lAxis-point1_lAxis)/(point2_vAxis-point1_vAxis);
    b=point1_lAxis;

end


