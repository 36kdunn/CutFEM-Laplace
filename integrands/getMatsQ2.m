function mats = getMatsQ2(N)
%
% Stored pre-computed matrices called in assembly of mass matrix and load
% vector.
% Intputs:
%   N = number of intervals on one side of the square domain Omega.
% Outputs:
%   mats = structure containing all pre-computed matrices.

% Values of basis functions evaluated at gauss quadrature points on [0,1]
x = 0.5*[-1/sqrt(3), -1/sqrt(3), 1/sqrt(3), 1/sqrt(3)] + 0.5;
y = 0.5*[-1/sqrt(3), 1/sqrt(3), -1/sqrt(3), 1/sqrt(3)] + 0.5;

phi0x = 2*(1-x).*(0.5-x);
phi5x = 4*(1-x).*x;
phi1x = 2*(x-0.5).*x;
phi0y = 2*(1-y).*(0.5-y);
phi5y = 4*(1-y).*y;
phi1y = 2*(y-0.5).*y;

mats.phiVals = [phi0x.*phi0y;
    phi0x.*phi5y;
    phi0x.*phi1y;
    phi5x.*phi0y;
    phi5x.*phi5y;
    phi5x.*phi1y;
    phi1x.*phi0y;
    phi1x.*phi5y;
    phi1x.*phi1y];

mats.uvOmega = [2/15, 1/15, -1/30;
    1/15, 8/15, 1/15;
    -1/30, 1/15, 2/15];
    
mats.graddotnvB = [2/5, -8/15, 2/15, 1/5, -4/15, 1/15, -1/10, 2/15, -1/30;
0, 0, 0, 0, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, 0, 0, 0, 0;
1/5, -4/15, 1/15, 8/5, -32/15, 8/15, 1/5, -4/15, 1/15;
0, 0, 0, 0, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, 0, 0, 0, 0;
-1/10, 2/15, -1/30, 1/5, -4/15, 1/15, 2/5, -8/15, 2/15;
0, 0, 0, 0, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, 0, 0, 0, 0];

mats.graddotnvT = [0, 0, 0, 0, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, 0, 0, 0, 0;
2/15, -8/15, 2/5, 1/15, -4/15, 1/5, -1/30, 2/15, -1/10;
0, 0, 0, 0, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, 0, 0, 0, 0;
1/15, -4/15, 1/5, 8/15, -32/15, 8/5, 1/15, -4/15, 1/5;
0, 0, 0, 0, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, 0, 0, 0, 0;
-1/30, 2/15, -1/10, 1/15, -4/15, 1/5, 2/15, -8/15, 2/5];

mats.graddotnvL = [2/5, 1/5, -1/10, -8/15, -4/15, 2/15, 2/15, 1/15, -1/30;
1/5, 8/5, 1/5, -4/15, -32/15, -4/15, 1/15, 8/15, 1/15;
-1/10, 1/5, 2/5, 2/15, -4/15, -8/15, -1/30, 1/15, 2/15;
0, 0, 0, 0, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, 0, 0, 0, 0];

mats.graddotnvR = [0, 0, 0, 0, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, 0, 0, 0, 0;
2/15, 1/15, -1/30, -8/15, -4/15, 2/15, 2/5, 1/5, -1/10;
1/15, 8/15, 1/15, -4/15, -32/15, -4/15, 1/5, 8/5, 1/5;
-1/30, 1/15, 2/15, 2/15, -4/15, -8/15, -1/10, 1/5, 2/5];

% Local stiffness matrix for Laplacian Q2
mats.gradugradv = [28/45, -1/5, -1/30, -1/5, -16/45, 1/9, -1/30, 1/9, -1/45;
-1/5, 88/45, -1/5, -16/45, -16/15, -16/45, 1/9, 0, 1/9;
-1/30, -1/5, 28/45, 1/9, -16/45, -1/5, -1/45, 1/9, -1/30;
-1/5, -16/45, 1/9, 88/45, -16/15, 0, -1/5, -16/45, 1/9;
-16/45, -16/15, -16/45, -16/15, 256/45, -16/15, -16/45, -16/15, -16/45;
1/9, -16/45, -1/5, 0, -16/15, 88/45, 1/9, -16/45, -1/5;
-1/30, 1/9, -1/45, -1/5, -16/45, 1/9, 28/45, -1/5, -1/30;
1/9, 0, 1/9, -16/45, -16/15, -16/45, -1/5, 88/45, -1/5;
-1/45, 1/9, -1/30, 1/9, -16/45, -1/5, -1/30, -1/5, 28/45];%*(1/N^2)

mats.uv = [4, 2, -1, 2, 1, -1/2, -1, -1/2, 1/4;
    2, 16, 2, 1, 8, 1, -1/2, -4, -1/2;
    -1, 2, 4, -1/2, 1, 2, 1/4, -1/2, -1;
    2, 1, -1/2, 16, 8, -4, 2, 1, -1/2;
    1, 8, 1, 8, 64, 8, 1, 8, 1;
    -1/2, 1, 2, -4, 8, 16, -1/2, 1, 2;
    -1, -1/2, 1/4, 2, 1, -1/2, 4, 2, -1;
    -1/2, -4, -1/2, 1, 8, 1, 2, 16, 2;
    1/4, -1/2, -1, -1/2, 1, 2, -1, 2, 4]/(225*N^2);

% Not multiplied by N because the ghost penalty term is multiplied by 1/N
mats.ghostvert = [2/15, 1/15, -1/30, -8/15, -4/15, 2/15, 4/5, 2/5, -1/5, -8/15, -4/15, 2/15, 2/15, 1/15, -1/30;
1/15, 8/15, 1/15, -4/15, -32/15, -4/15, 2/5, 16/5, 2/5, -4/15, -32/15, -4/15, 1/15, 8/15, 1/15;
-1/30, 1/15, 2/15, 2/15, -4/15, -8/15, -1/5, 2/5, 4/5, 2/15, -4/15, -8/15, -1/30, 1/15, 2/15;
-8/15, -4/15, 2/15, 32/15, 16/15, -8/15, -16/5, -8/5, 4/5, 32/15, 16/15, -8/15, -8/15, -4/15, 2/15;
-4/15, -32/15, -4/15, 16/15, 128/15, 16/15, -8/5, -64/5, -8/5, 16/15, 128/15, 16/15, -4/15, -32/15, -4/15;
2/15, -4/15, -8/15, -8/15, 16/15, 32/15, 4/5, -8/5, -16/5, -8/15, 16/15, 32/15, 2/15, -4/15, -8/15;
4/5, 2/5, -1/5, -16/5, -8/5, 4/5, 24/5, 12/5, -6/5, -16/5, -8/5, 4/5, 4/5, 2/5, -1/5;
2/5, 16/5, 2/5, -8/5, -64/5, -8/5, 12/5, 96/5, 12/5, -8/5, -64/5, -8/5, 2/5, 16/5, 2/5;
-1/5, 2/5, 4/5, 4/5, -8/5, -16/5, -6/5, 12/5, 24/5, 4/5, -8/5, -16/5, -1/5, 2/5, 4/5;
-8/15, -4/15, 2/15, 32/15, 16/15, -8/15, -16/5, -8/5, 4/5, 32/15, 16/15, -8/15, -8/15, -4/15, 2/15;
-4/15, -32/15, -4/15, 16/15, 128/15, 16/15, -8/5, -64/5, -8/5, 16/15, 128/15, 16/15, -4/15, -32/15, -4/15;
2/15, -4/15, -8/15, -8/15, 16/15, 32/15, 4/5, -8/5, -16/5, -8/15, 16/15, 32/15, 2/15, -4/15, -8/15;
2/15, 1/15, -1/30, -8/15, -4/15, 2/15, 4/5, 2/5, -1/5, -8/15, -4/15, 2/15, 2/15, 1/15, -1/30;
1/15, 8/15, 1/15, -4/15, -32/15, -4/15, 2/5, 16/5, 2/5, -4/15, -32/15, -4/15, 1/15, 8/15, 1/15;
-1/30, 1/15, 2/15, 2/15, -4/15, -8/15, -1/5, 2/5, 4/5, 2/15, -4/15, -8/15, -1/30, 1/15, 2/15];

% Not multiplied by N because the ghost penalty term is multiplied by 1/N
mats.ghosthoriz = [2/15, -8/15, 4/5, -8/15, 2/15, 1/15, -4/15, 2/5, -4/15, 1/15, -1/30, 2/15, -1/5, 2/15, -1/30;
-8/15, 32/15, -16/5, 32/15, -8/15, -4/15, 16/15, -8/5, 16/15, -4/15, 2/15, -8/15, 4/5, -8/15, 2/15;
4/5, -16/5, 24/5, -16/5, 4/5, 2/5, -8/5, 12/5, -8/5, 2/5, -1/5, 4/5, -6/5, 4/5, -1/5;
-8/15, 32/15, -16/5, 32/15, -8/15, -4/15, 16/15, -8/5, 16/15, -4/15, 2/15, -8/15, 4/5, -8/15, 2/15;
2/15, -8/15, 4/5, -8/15, 2/15, 1/15, -4/15, 2/5, -4/15, 1/15, -1/30, 2/15, -1/5, 2/15, -1/30;
1/15, -4/15, 2/5, -4/15, 1/15, 8/15, -32/15, 16/5, -32/15, 8/15, 1/15, -4/15, 2/5, -4/15, 1/15;
-4/15, 16/15, -8/5, 16/15, -4/15, -32/15, 128/15, -64/5, 128/15, -32/15, -4/15, 16/15, -8/5, 16/15, -4/15;
2/5, -8/5, 12/5, -8/5, 2/5, 16/5, -64/5, 96/5, -64/5, 16/5, 2/5, -8/5, 12/5, -8/5, 2/5;
-4/15, 16/15, -8/5, 16/15, -4/15, -32/15, 128/15, -64/5, 128/15, -32/15, -4/15, 16/15, -8/5, 16/15, -4/15;
1/15, -4/15, 2/5, -4/15, 1/15, 8/15, -32/15, 16/5, -32/15, 8/15, 1/15, -4/15, 2/5, -4/15, 1/15;
-1/30, 2/15, -1/5, 2/15, -1/30, 1/15, -4/15, 2/5, -4/15, 1/15, 2/15, -8/15, 4/5, -8/15, 2/15;
2/15, -8/15, 4/5, -8/15, 2/15, -4/15, 16/15, -8/5, 16/15, -4/15, -8/15, 32/15, -16/5, 32/15, -8/15;
-1/5, 4/5, -6/5, 4/5, -1/5, 2/5, -8/5, 12/5, -8/5, 2/5, 4/5, -16/5, 24/5, -16/5, 4/5;
2/15, -8/15, 4/5, -8/15, 2/15, -4/15, 16/15, -8/5, 16/15, -4/15, -8/15, 32/15, -16/5, 32/15, -8/15;
-1/30, 2/15, -1/5, 2/15, -1/30, 1/15, -4/15, 2/5, -4/15, 1/15, 2/15, -8/15, 4/5, -8/15, 2/15];


mats.ghostvert2ndOrder = [32/15, 16/15, -8/15, -64/15, -32/15, 16/15, 0, 0, 0, 64/15, 32/15, -16/15, -32/15, -16/15, 8/15;
16/15, 128/15, 16/15, -32/15, -256/15, -32/15, 0, 0, 0, 32/15, 256/15, 32/15, -16/15, -128/15, -16/15;
-8/15, 16/15, 32/15, 16/15, -32/15, -64/15, 0, 0, 0, -16/15, 32/15, 64/15, 8/15, -16/15, -32/15;
-64/15, -32/15, 16/15, 128/15, 64/15, -32/15, 0, 0, 0, -128/15, -64/15, 32/15, 64/15, 32/15, -16/15;
-32/15, -256/15, -32/15, 64/15, 512/15, 64/15, 0, 0, 0, -64/15, -512/15, -64/15, 32/15, 256/15, 32/15;
16/15, -32/15, -64/15, -32/15, 64/15, 128/15, 0, 0, 0, 32/15, -64/15, -128/15, -16/15, 32/15, 64/15;
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
64/15, 32/15, -16/15, -128/15, -64/15, 32/15, 0, 0, 0, 128/15, 64/15, -32/15, -64/15, -32/15, 16/15;
32/15, 256/15, 32/15, -64/15, -512/15, -64/15, 0, 0, 0, 64/15, 512/15, 64/15, -32/15, -256/15, -32/15;
-16/15, 32/15, 64/15, 32/15, -64/15, -128/15, 0, 0, 0, -32/15, 64/15, 128/15, 16/15, -32/15, -64/15;
-32/15, -16/15, 8/15, 64/15, 32/15, -16/15, 0, 0, 0, -64/15, -32/15, 16/15, 32/15, 16/15, -8/15;
-16/15, -128/15, -16/15, 32/15, 256/15, 32/15, 0, 0, 0, -32/15, -256/15, -32/15, 16/15, 128/15, 16/15;
8/15, -16/15, -32/15, -16/15, 32/15, 64/15, 0, 0, 0, 16/15, -32/15, -64/15, -8/15, 16/15, 32/15];

mats.ghosthoriz2ndOrder = [32/15, -64/15, 0, 64/15, -32/15, 16/15, -32/15, 0, 32/15, -16/15, -8/15, 16/15, 0, -16/15, 8/15;
-64/15, 128/15, 0, -128/15, 64/15, -32/15, 64/15, 0, -64/15, 32/15, 16/15, -32/15, 0, 32/15, -16/15;
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
64/15, -128/15, 0, 128/15, -64/15, 32/15, -64/15, 0, 64/15, -32/15, -16/15, 32/15, 0, -32/15, 16/15;
-32/15, 64/15, 0, -64/15, 32/15, -16/15, 32/15, 0, -32/15, 16/15, 8/15, -16/15, 0, 16/15, -8/15;
16/15, -32/15, 0, 32/15, -16/15, 128/15, -256/15, 0, 256/15, -128/15, 16/15, -32/15, 0, 32/15, -16/15;
-32/15, 64/15, 0, -64/15, 32/15, -256/15, 512/15, 0, -512/15, 256/15, -32/15, 64/15, 0, -64/15, 32/15;
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
32/15, -64/15, 0, 64/15, -32/15, 256/15, -512/15, 0, 512/15, -256/15, 32/15, -64/15, 0, 64/15, -32/15;
-16/15, 32/15, 0, -32/15, 16/15, -128/15, 256/15, 0, -256/15, 128/15, -16/15, 32/15, 0, -32/15, 16/15;
-8/15, 16/15, 0, -16/15, 8/15, 16/15, -32/15, 0, 32/15, -16/15, 32/15, -64/15, 0, 64/15, -32/15;
16/15, -32/15, 0, 32/15, -16/15, -32/15, 64/15, 0, -64/15, 32/15, -64/15, 128/15, 0, -128/15, 64/15;
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
-16/15, 32/15, 0, -32/15, 16/15, 32/15, -64/15, 0, 64/15, -32/15, 64/15, -128/15, 0, 128/15, -64/15;
8/15, -16/15, 0, 16/15, -8/15, -16/15, 32/15, 0, -32/15, 16/15, -32/15, 64/15, 0, -64/15, 32/15];
