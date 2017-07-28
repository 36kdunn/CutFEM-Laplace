% test integrateInterfaceElements
close all

% Add subfolders to path
addpath('BCs')
addpath('integrands')
addpath('RHS')
addpath('utilities')
% addpath('utils')

%load('ic16.mat')

% Test the ability to find (xp,yp) and (xn,yn)
Tol = 1e-15;
N = 32; % h = 1/N
% G = [flipud([0.6;0.65;0.35;0.4;0.45;0.5;0.55;0.6]),flipud(0.7*[1;1;1;1;1;1;1;1])];
% G = flipud(G);
theta = (0:pi/80:2*pi)';
G = [0.23*cos(theta) + 0.5, 0.23*sin(theta) + 0.5];
h_s = max(sqrt(sum((G(1:end-1,:)-G(2:end,:)).^2,2)));
h = [1/N,h_s];
% G = [0.41*cos(theta) + 0.5, 0.31*sin(theta) + 0.5];
% Gstep = sqrt((G(2,1)-G(1,1))^2+(G(2,2)-G(1,2))^2);

% Cardioid plus
%{
G = [0.35*cos(theta).*(1-sin(theta)) + 0.5, 0.15*sin(theta).*(1-sin(theta)) + 0.7];
G1 = [(0.15*cos(theta).*(1-cos(theta)) + 0.7), 0.35*sin(theta).*(1-cos(theta)) + 0.5];
G = (G+G1)/2;
Gstep = sqrt((G(2,1)-G(1,1))^2+(G(2,2)-G(1,2))^2);
%}

% stepRatio = Gstep/N;

%{
plot(G(:,1),G(:,2),'*-')
ax = gca;
grid on
%}

intMeanZero = 0;
extMeanZero = 0;
intNeum = 1;
extNeumGamma = 1;
extNeumOmega = 0;
jump = 0;
derivJump = 0;
extGamCouple = 1;
intGamCouple = 1;

flags = [intMeanZero, extMeanZero, intNeum, extNeumGamma, extNeumOmega,...
    extGamCouple, intGamCouple];

[extsln,intsln,x,y] = intExtFEM(N,G,flags);

%{
% f = @(s) [0.3*cos(s) + 0.505, 0.23*sin(s) + 0.55];
f = @(s) [0.28*cos(s) + 0.5, 0.23*sin(s) + 0.5];

h = 0.1;
searchstep = h/100;

intCoords = findInterfaceIntersectionsFunc(f,h,searchstep);

t = 0:pi/100:2*pi;
fvals = f(t(:));

fintvals = f(intCoords);

figure(1)
plot(fvals(:,1),fvals(:,2),fintvals(:,1),fintvals(:,2),'ro')
grid on;

% intOverCut(0.7,0.5,0.1,f,intCoords([19:20,1:2]),Tol)
intOverCut(0.7,0.5,0.1,f,intCoords,Tol)

%}

