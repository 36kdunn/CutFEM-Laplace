function I = gaussEdge(f,x1,x2,y1,y2)
%
% Gauss quadrature over edge of square with limits [x1,x2], [y1,y2]

% Two point quadrature
% gaussPts = [-1/sqrt(3), 1/sqrt(3)];
% w = [1,1];

% Three point quadrature
% gaussPts = [-sqrt(3/5), 0, sqrt(3/5)];
% w = [5/9, 8/9, 5/9];

% Four point quadrature
gaussPts = [-sqrt(3/7 + (2/7)*sqrt(6/5)), -sqrt(3/7 - (2/7)*sqrt(6/5)), sqrt(3/7 - (2/7)*sqrt(6/5)), sqrt(3/7 + (2/7)*sqrt(6/5))];
w = [(18-sqrt(30))/36, (18+sqrt(30))/36, (18+sqrt(30))/36, (18-sqrt(30))/36];

dx = x2-x1;
dy = y2-y1;
ds = sqrt(dx^2+dy^2);

% Compute translated Gaussian quadrature points
gaussX = 0.5*(gaussPts*dx + 2*x1 + dx);
gaussY = 0.5*(gaussPts*dy + 2*y1 + dy);

I = 0;
for i = 1:length(gaussPts)

    I = I + w(1,i)*f(gaussX(1,i),gaussY(1,i));
    
end

I = ds*I/2; % area (length in this case) = 1/N

