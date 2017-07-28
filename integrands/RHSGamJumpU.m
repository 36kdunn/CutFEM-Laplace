function f = RHSGamJumpU(x,y,x0,y0,N,n)

phi0x = 2*(1-N*(x-x0)).*(0.5-N*(x-x0));
phi5x = 4*(1-N*(x-x0)).*(N*(x-x0));
phi1x = 2*(N*(x-x0)-0.5).*(N*(x-x0));

phi0y = 2*(1-N*(y-y0)).*(0.5-N*(y-y0));
phi5y = 4*(1-N*(y-y0)).*(N*(y-y0));
phi1y = 2*(N*(y-y0)-0.5).*(N*(y-y0));

phi0xprime = 2*(-N)*(0.5-N*(x-x0)) + 2*(1-N*(x-x0)).*(-N);
phi5xprime = 4*(-N)*(N*(x-x0)) + 4*(N)*(1-N*(x-x0));
phi1xprime = 2*(N)*(N*(x-x0)) + 2*(N)*(N*(x-x0)-0.5);

phi0yprime = 2*(-N)*(0.5-N*(y-y0)) + 2*(1-N*(y-y0)).*(-N);
phi5yprime = 4*(-N)*(N*(y-y0)) + 4*(N)*(1-N*(y-y0));
phi1yprime = 2*(N)*(N*(y-y0)) + 2*(N)*(N*(y-y0)-0.5);

gradVecX = [phi0xprime.*phi0y;
    phi0xprime.*phi5y;
    phi0xprime.*phi1y;
    phi5xprime.*phi0y;
    phi5xprime.*phi5y;
    phi5xprime.*phi1y;
    phi1xprime.*phi0y;
    phi1xprime.*phi5y;
    phi1xprime.*phi1y];

gradVecY = [phi0yprime.*phi0x;
    phi5yprime.*phi0x;
    phi1yprime.*phi0x;
    phi0yprime.*phi5x;
    phi5yprime.*phi5x;
    phi1yprime.*phi5x;
    phi0yprime.*phi1x;
    phi5yprime.*phi1x;
    phi1yprime.*phi1x];

f = zeros(9,1);

% Result of the inner product (g,grad(eps(v),n))
f =  fjumpu(x,y)*gradVecX*n(1,1) + fjumpu(x,y)*n(1,2)*gradVecY;

