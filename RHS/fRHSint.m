function f = fRHSint(x,y)

% EX 1: EXPONENTIAL/TRIG
f = -(-2*y + x.*y.^2 + x.^3).*exp(-x.*y);

end