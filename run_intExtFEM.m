% test integrateInterfaceElements
close all

% Add subfolders to path
addpath('BCs')
addpath('integrands')
addpath('RHS')
addpath('utils')


% Test the ability to find (xp,yp) and (xn,yn)
Tol = 1e-15;
N = 32; % h = 1/N

% Define interface as circle
theta = (0:pi/80:2*pi)';
G = [0.23*cos(theta) + 0.5, 0.23*sin(theta) + 0.5];
h_s = max(sqrt(sum((G(1:end-1,:)-G(2:end,:)).^2,2)));
h = [1/N,h_s];
% G = [0.41*cos(theta) + 0.5, 0.31*sin(theta) + 0.5];
% Gstep = sqrt((G(2,1)-G(1,1))^2+(G(2,2)-G(1,2))^2);

% Define interface as sum of two cardioids
%{
G = [0.35*cos(theta).*(1-sin(theta)) + 0.5, 0.15*sin(theta).*(1-sin(theta)) + 0.7];
G1 = [(0.15*cos(theta).*(1-cos(theta)) + 0.7), 0.35*sin(theta).*(1-cos(theta)) + 0.5];
G = (G+G1)/2;
Gstep = sqrt((G(2,1)-G(1,1))^2+(G(2,2)-G(1,2))^2);
%}

%{
plot(G(:,1),G(:,2),'*-')
ax = gca;
grid on
%}

% Specify whether we choose interior/exterior solutions to have 0 mean
intMeanZero = 0;
extMeanZero = 0;

flags = [intMeanZero, extMeanZero];

[extsln,intsln,x,y] = intExtFEM(N,G,flags);
