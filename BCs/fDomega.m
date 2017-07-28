function f = fDomega(x,y)

% EX 1: CONSTANT
%f = ((1+k2)*(x.^2+3*x+y.^2)-4)/k1;
% f=0;

% EX 2: LINEAR
% f = 3+x+y-4;

% EX 3: EXPONENTIALS/TRIG
f = x.*exp(-x.*y);

% f = 10*((x-1/2).^2 + (y-1/2).^2);

% EX 4: RBM Test Example
% f = 2-2*y;

% EX 4: THREE COUPLED PROBLEMS, FUNDAMENTAL SLN
% c = -log(sqrt(8))/(2*pi);
% scale = pi/log(sqrt(2));
% f = scale*(-c-log(sqrt((x-2).^2+(y-2).^2))/(2*pi));

% f = 0;

end