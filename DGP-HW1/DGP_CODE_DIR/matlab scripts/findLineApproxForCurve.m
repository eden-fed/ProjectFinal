function [ m,intersectionX ] = findLineApproxForCurve( sigma,SIGMA,k)
    %use the return value as y>=a*x+b.
    point1_vAxis=0;
    point1_lAxis=-log(1/sigma);
    
	point2_vAxis=-(lambertw((sigma/SIGMA)*exp(1))-1);
    if(k<=point2_vAxis)
        point2_vAxis=k;
    end
    point2_lAxis=-log((1-point2_vAxis)/sigma);
    
    if(point2_vAxis-point1_vAxis==0)
        m=1;
    else
        m=(point2_lAxis-point1_lAxis)/(point2_vAxis-point1_vAxis);
    end
    intersectionX=point2_vAxis;

end


