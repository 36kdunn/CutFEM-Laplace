function f = RHSBulkExt(x,y,x0,y0,N)

f = fRHSext(x,y).*basisFuncsVec(x,y,x0,y0,N);