function [ x_poly,y_poly ] = projectToPoly( x,y,SIGMA,sigma,k,intersectionX )
%x = vector of x coordinates of the given points
%y = vector of y coordinate of the given points
%intersectionX = x coordinate of right point for the line approximation
if(k==intersectionX)
    [ x_poly,y_poly ] = projectToPolywithK( x,y,SIGMA,sigma,k);
else
    [ x_poly,y_poly ] = projectToPolyNok(x,y,SIGMA,sigma,intersectionX );
end

end

