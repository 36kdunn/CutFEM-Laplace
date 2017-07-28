function f = fNint(x,y,n)

% EX 1
strainTensor = [3*x.^2 - 3, x+y; x+y, 2*y];

% EX 2: LINEAR
% strainTensor = [1, 1;1, 2];

strainTensor = [3, 0; 0, 3];

f = strainTensor*n(:);

end