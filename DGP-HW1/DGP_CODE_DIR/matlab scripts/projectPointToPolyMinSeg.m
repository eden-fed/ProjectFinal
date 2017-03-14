function [x_poly,y_poly] = projectPointToPolyMinSeg(x, y, xv, yv) 
%x = x coordinate of the given point
%y = y coordinate of the given point
%xv = vector of x coordinates of the vertices of the polygon
%yv = vector of y coordinates of the vertices of the polygon

%use cross product to determine if inside or outside polygon
I=(xv(2:end)-xv(1:end-1)).*(y-yv(1:end-1))-(yv(2:end)-yv(1:end-1)).*(x-xv(1:end-1));
if(all(I<0))
    x_poly=x;
    y_poly=y;
else
    %find the closest point on each line segment
    dot=(x-xv(1:end-1)).*(xv(2:end)-xv(1:end-1))+(y-yv(1:end-1)).*(yv(2:end)-yv(1:end-1));
    projectionOnLine=dot./(((xv(2:end)-xv(1:end-1)).^2+(yv(2:end)-yv(1:end-1)).^2));
    t=max(0,min(1,projectionOnLine));

    closestPointOnEachSegX=xv(1:end-1)+t.*(xv(2:end)-xv(1:end-1));
    closestPointOnEachSegY=yv(1:end-1)+t.*(yv(2:end)-yv(1:end-1));

    %find the minimun distace from all the line segments
    [~,minIdx]=min((x-closestPointOnEachSegX).^2+(y-closestPointOnEachSegY).^2);
    x_poly=closestPointOnEachSegX(minIdx);
    y_poly=closestPointOnEachSegY(minIdx);
end
end

