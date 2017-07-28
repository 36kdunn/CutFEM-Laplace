function I = intOverPoly(f,verts,Tol)
%
% Computes the integral of a general simple polygon by breaking the polygon
% into triangles using an ear-clipping method.  When three vertices are
% chosen, the remaining vertices are checked to make sure they are not
% inside the triangle formed by the chosen vertices.  We also ensure that
% no errors are caused by three colinear vertices or repeated vertices.
% The value of the integral over each triangle is computed using a
% three-point Gaussian quadrature formula
%
% Input:
%   f: function to be integrated over the polygon
%   verts: vertices of polygon - MUST BE ORDERED COUNTERCLOCKWISE!!!
%   Tol: tolerance for two vertices to be "equal" ( abs(x1-x2)<Tol )
%
% Output:
%   I: value of integral over polygon

% verts = [0.15,0;0,0;0.5,0;0.5,0.5;0.5,0.5;0.25,0.05]

% Initialize value of integral and vertex index
I = 0;
vertInd = 1;
% Loop over all vertices
while numel(verts)/2 > 2
    
    % Avoid index out of bounds erros by checking indices
    if vertInd > numel(verts)/2
        vertInd = 1;
    end
    
    indplus1 = vertInd+1;
    indplus2 = vertInd+2;
    
    if indplus1 == numel(verts)/2+1
        indplus1 = 1;
        indplus2 = 2;
    end
    if indplus2 == numel(verts)/2+1
        indplus2 = 1;
    end
    
    % Check three indices of interest for colinearity and repeated vertices
    [verts,vertInd,indplus1,indplus2] = checkForProblems(verts,indplus1,vertInd,indplus2,Tol);
    
    % Assign values of vertices to (x1,y1), (x2,y2), and (x3,y3) for simplicity
    x1 = verts(vertInd,1);
    y1 = verts(vertInd,2);
    
    x2 = verts(indplus1,1);
    y2 = verts(indplus1,2);
    
    x3 = verts(indplus2,1);
    y3 = verts(indplus2,2);
    
    % Compute INWARD normal between each line segment
    n12 = -[y2-y1, -(x2-x1)];
    n23 = -[y3-y2, -(x3-x2)];
    
    % Compute midpoint of first segment
    xbar1 = (x1+x2)/2;
    ybar1 = (y1+y2)/2;
    % Compute midpoint of second segment
    xbar2 = (x2+x3)/2;
    ybar2 = (y2+y3)/2;
    
    % See if the rays eminating from (xbar1,y1) and (x3,y3) in the direction
    % of n12, n23 respectively intersect.
    dx = xbar2 - xbar1;
    dy = ybar2 - ybar1;
    D = n23(1)*n12(2)-n23(2)*n12(1);
    sgn1 = sign(dy * n23(1) - dx * n23(2))*sign(D);
    sgn2 = sign(dy * n12(1) - dx * n12(2))*sign(D);
    
    % The area inside triangle is part of the area of interest if either
    % sgn1 or sgn2 is nonnegative. This means that
    % (x,y) = t1*n12 + (xbar1,ybar1) = t2*n23 + (xbar2,ybar2)
    % for t1 or t2 nonnegative
    % If D == 0, then n12, n23 are parallel and do not intersect
    if (sgn1 > 0 || sgn2 > 0)% && D ~= 0
        % If the interior of the triangle formed by the three vertices is
        % in our area of interest
        flag = 0;
        ind = indplus2+1;
        if ind > numel(verts)/2
            ind = 1;
        end
        while ind ~= vertInd
            % Makes sure that every remaining point is outside the triangle
            % v0 and v1 form a basis of R^2. Find u,v s.t. v2 = u*v0 + v*v1
            v0 = [x2-x1,y2-y1];
            v1 = [x3-x1,y3-y1];
            v2 = [verts(ind,1)-x1,verts(ind,2)-y1];
            
            d = dot(v0,v0)*dot(v1,v1)-dot(v0,v1)^2;
            u = (dot(v2,v0)*dot(v1,v1)-dot(v0,v1)*dot(v2,v1))/d;
            v = (dot(v0,v0)*dot(v2,v1)-dot(v0,v1)*dot(v2,v0))/d;
            
            if ((u>=0 && u<=1) && (v>=0 && v<=1)) && ((u+v)<=1 && (u+v)>=0)
                % If, for u,v described above,
                %   (1) 0<=u<=1 and 0<=v<=1, and
                %   (2) 0<=u+v<=1
                % then the point lies within or on edges of the triangle. 
                % Change flag to 1 and we no longer need to check other vertices.
                flag = 1;
                break
            end
            
            % Else, increase loop index and ensure it does not exceed bounds
            ind = ind+1;
            if ind > numel(verts)/2
                ind = 1;
            end
        end
        
        % If no vertices were found, integrate over triangle and eliminate
        % middle vertex
        if flag == 0
            I = I + intOverTri(f,verts([vertInd;indplus1;indplus2],:));
            verts = verts([1:indplus1-1,indplus1+1:numel(verts)/2],:);
        else
            % Otherwise, increase vertex index by 1 and try new triangle
            vertInd = vertInd + 1;
            if vertInd > numel(verts)/2
                vertInd = 1;
            end
        end
    else
        % If this is not our area of interest, increase vertex index by 1
        % and try new triangle
        vertInd = vertInd + 1;
        if vertInd > numel(verts)/2
            vertInd = 1;
        end
    end
    
end

end

function [verts,prevInd,ind,nextInd] = checkForProblems(verts,ind,prevInd,nextInd,Tol)
%
% Checks for colinearity and repeated vertex coordinates.  Both cases may
% cause problems in integration function


% Compute outward normals
n12 = [verts(ind,2)-verts(prevInd,2), -(verts(ind,1)-verts(prevInd,1))];
n12 = n12/norm(n12); % normalize
n23 = [verts(nextInd,2)-verts(ind,2), -(verts(nextInd,1)-verts(ind,1))];
n23 = n23/norm(n23); % normalize

% Check to make sure we don't have three colinear vertices
if norm(n12-n23) < Tol || norm(n12+n23) < Tol
    % If normals point in the same direction or in opposite directions
    % then three vertices are colinear; either 1---2---3 or 2---1---3
    verts = verts([1:ind-1,ind+1:numel(verts)/2],:);
end

% Update indices in case the length of verts changed in last conditional
if ind > numel(verts)/2
    ind = 1;
end

prevInd = ind - 1;
if prevInd==0
    prevInd = numel(verts)/2;
end

nextInd = ind + 1;
if nextInd > numel(verts)/2
    nextInd = 1;
end

% Check to make sure we don't have a repeated vertex
if verts(ind,1) == verts(nextInd,1) && verts(ind,2) == verts(nextInd,2)
    verts = verts([1:ind-1,ind+1:numel(verts)/2],:);
end

end

function I = intOverTri(f,verts)
%
% Integrates general triangle with vertices given in verts.
%
% Inputs:
%   f: function to be integrated over triangular region
%   verts: vertices of triangle
%
% Outputs:
%   I: value of integral over triangular region

% Define Gaussian quadrature points in reference triangle
% gaussPts = [1/6, 1/6; 2/3, 1/6; 1/6, 2/3]';
% w = ones(3,1)/3;

gaussPts = [1/3, 1/3; 1/5, 3/5; 1/5, 1/5; 3/5, 1/5]';
w = [-27/48; 25/48; 25/48; 25/48];

%{
gaussPts = [0.33333333333333 0.33333333333333;
    0.47014206410511 0.47014206410511;
    0.47014206410511 0.05971587178977;
    0.05971587178977 0.47014206410511;
    0.10128650732346 0.10128650732346;
    0.10128650732346 0.79742698535309;
    0.79742698535309 0.10128650732346]';

 w=[0.22500000000000;
 0.13239415278851;
 0.13239415278851;
 0.13239415278851;
 0.12593918054483;
 0.12593918054483;
 0.12593918054483];
%}
%{
gaussPts = [0.24928674517091 0.24928674517091;
    0.24928674517091 0.50142650965818;
    0.50142650965818 0.24928674517091;
    0.06308901449150 0.06308901449150;
    0.06308901449150 0.87382197101700;
    0.87382197101700 0.06308901449150;
    0.31035245103378 0.63650249912140;
    0.63650249912140 0.05314504984482;
    0.05314504984482 0.31035245103378;
    0.63650249912140 0.31035245103378;
    0.31035245103378 0.05314504984482;
    0.05314504984482 0.63650249912140]';

w=[ 0.11678627572638;
 0.11678627572638;
 0.11678627572638;
 0.05084490637021;
 0.05084490637021;
 0.05084490637021;
 0.08285107561837;
 0.08285107561837;
 0.08285107561837;
 0.08285107561837;
 0.08285107561837;
 0.08285107561837];
%}

% Create matrix/vector for affine transformation of ref. coords [xhat,yhat]
% [x;y] = T*[xhat;yhat] + b where (xhat,yhat)
T = [verts(2,1)-verts(1,1), verts(3,1)-verts(1,1);
    verts(2,2)-verts(1,2), verts(3,2)-verts(1,2)];
b = [verts(1,1);verts(1,2)];

I = 0;
transGaussPts = zeros(2,size(gaussPts,2));
for i = 1:size(gaussPts,2)
    % Translate each pair of quadrature points to triangle given by verts
    transGaussPts(:,i) = T*gaussPts(:,i) + b;
    % Integral is approximately equal to the sum of function evaluations
    % multiplied by the area of triangle = abs(det(T))/2 divided by 3
    I = I + f(transGaussPts(1,i),transGaussPts(2,i))*(abs(det(T))/2)*w(i,1);
end

end
