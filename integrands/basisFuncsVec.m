function f = basisFuncsVec(x,y,x0,y0,N)

phi0x = 2*(1-N*(x-x0)).*(0.5-N*(x-x0));
phi5x = 4*(1-N*(x-x0)).*(N*(x-x0));
phi1x = 2*(N*(x-x0)-0.5).*(N*(x-x0));

phi0y = 2*(1-N*(y-y0)).*(0.5-N*(y-y0));
phi5y = 4*(1-N*(y-y0)).*(N*(y-y0));
phi1y = 2*(N*(y-y0)-0.5).*(N*(y-y0));

phi00 = phi0x.*phi0y;
phi05 = phi0x.*phi5y;
phi01 = phi0x.*phi1y;

phi50 = phi5x.*phi0y;
phi55 = phi5x.*phi5y;
phi51 = phi5x.*phi1y;

phi10 = phi1x.*phi0y;
phi15 = phi1x.*phi5y;
phi11 = phi1x.*phi1y;

f = [phi00;phi05;phi01;phi50;phi55;phi51;phi10;phi15;phi11];