function f = RHSBulkInt(x,y,x0,y0,N)

f = fRHSint(x,y).*basisFuncsVec(x,y,x0,y0,N);