%*******************************************************************************
% function:	p_poly_dist
% Description:	distance from point to polygon whose vertices are specified by the
%              vectors xv and yv
% Input:  
%    x - point's x coordinate
%    y - point's y coordinate
%    xv - vector of polygon vertices x coordinates
%    yv - vector of polygon vertices x coordinates
% Output: 
%    d - distance from point to polygon (defined as a minimal distance from 
%        point to any of polygon's ribs, positive if the point is outside the
%        polygon and negative otherwise)
%    x_poly: x coordinate of the point in the polygon closest to x,y
%    y_poly: y coordinate of the point in the polygon closest to x,y
%
% Routines: p_poly_dist.m
% Revision history:
%    03/31/2008 - return the point of the polygon closest to x,y
%               - added the test for the case where a polygon rib is 
%                 either horizontal or vertical. From Eric Schmitz.
%               - Changes by Alejandro Weinstein
%    7/9/2006  - case when all projections are outside of polygon ribs
%    23/5/2004 - created by Michael Yoshpe 
% Remarks:
%*******************************************************************************
function [x_poly,y_poly] = p_poly_dist_vec(x, y, xv, yv) 

% If (xv,yv) is not closed, close it.
xv = xv(:);
yv = yv(:);
Nv = length(xv);
xv = [xv ; xv(1)];
yv = [yv ; yv(1)];


% linear parameters of segments that connect the vertices
% Ax + By + C = 0
A = -diff(yv);
B =  diff(xv);
C = yv(2:end).*xv(1:end-1) - xv(2:end).*yv(1:end-1);

% find the projection of point (x,y) on each rib
Np=length(x);
AB = 1./(A.^2 + B.^2);
vv = (A*x+B*y+repmat(C,1,Np));
rep_x=repmat(x,Nv,1);
rep_y=repmat(y,Nv,1);
xp = rep_x - repmat((A.*AB),1,Np).*vv;
yp = rep_y - repmat((B.*AB),1,Np).*vv;

% Test for the case where a polygon rib is 
% either horizontal or vertical. From Eric Schmitz
id = find(diff(xv)==0);
xp(id,:)=repmat(xv(id),1,Np);
clear id
id = find(diff(yv)==0);
yp(id,:)=repmat(yv(id),1,Np);

% find all cases where projected point is inside the segment
rep_xv_1=repmat(xv(1:end-1),1,Np);
rep_xv_2=repmat(xv(2:end),1,Np);
rep_yv_1=repmat(yv(1:end-1),1,Np);
rep_yv_2=repmat(yv(2:end),1,Np);
idx_x = (((xp>=rep_xv_1) & (xp<=rep_xv_2)) | ((xp>=rep_xv_2) & (xp<=rep_xv_1)));
idx_y = (((yp>=rep_yv_1) & (yp<=rep_yv_2)) | ((yp>=rep_yv_2) & (yp<=rep_yv_1)));
idx = idx_x & idx_y;

% distance from point (x,y) to the vertices
dv = sqrt((rep_xv_1-rep_x).^2 + (rep_yv_1-rep_y).^2);

x_poly=zeros(size(x));
y_poly=zeros(size(x));

if(~any(idx)) % all projections are outside of polygon ribs
   [~,I] = min(dv);
   x_poly = xv(I);
   y_poly = yv(I);
else
   % distance from point (x,y) to the projection on ribs
   dp = sqrt((xp.*idx-rep_x.*idx).^2 + (yp.*idx-rep_y.*idx).^2);
   dp(dp==0)=Nan;
   [min_dv,I1] = min(dv);
   [min_dp,I2] = min(dp);
   [~,I] = min([min_dv; min_dp]);
   
    xx=1:length(x);
    
    x_poly(I==1)=rep_xv_1(sub2ind(size(dv), I1(I==1), xx(I==1)));
    x_poly(I==2)=xp(sub2ind(size(dp), I2(I==2), xx(I==2)));
    y_poly(I==1)=rep_yv_1(sub2ind(size(dv), I1(I==1), xx(I==1)));
    y_poly(I==2)=yp(sub2ind(size(dp), I2(I==2), xx(I==2)));

end

end