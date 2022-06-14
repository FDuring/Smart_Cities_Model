% Cube with 2 walls and feed-back
clc, clear all

%Physical values
%****************************
Kp = 1e5;              %P-controller gain: large for precision
Sc = 5*3*3; Si = Sc; Sg = 3*3; %surface [m2]: concrete, insulation, glass
Va = 3*3*3;             %air volume[m3]
rhoa = 1.2; ca = 1000;  %indoor air density; heat capacity
Vpa = 1*Va/3600;        %infiltration and ventilation air: volume/hour

% c: concrete; i: insulation;  g: glass
lamc = 2;       lami = 0.04;    lamg = 1.2;   %[W/m K]
rhoccc = 2.5e6; rhoici = 0.02e6; rhogcg = 2.0e6; %[J/m3 K]
wc = 0.2;       wi = 0.08;      wg = 0.01;    %[m]
epswLW = 0.9;   %long wave wall emmisivity
epswSW = 0.5;   %short wave wall emmisivity

epsgLW = 0.6;   %long wave glass emmisivity
taugSW = 0.7;   %short wave glass transmitance
alphagSW = 0.1; %short wave glass absortivity

sigma = 5.67e-8;%[W/m2 K4]
Fwg = 1/5;      %view factor wall - glass
Tm = 20 + 273;  %mean temp for radiative exchange

% convection coefficients
ho = 10;    hi = 4;   %[W/m2 K]
%***************************

% Conductances and capacities
Gc = lamc/wc*Sc; Cc = Sc*wc*rhoccc; %concrete
Gi = lami/wi*Si; Ci = Si*wi*rhoici; %insulation
Gg = lamg/wg*Sg; Cg = Sg*wg*rhogcg; %glass
Ca = Va*rhoa*ca;
% Convection
Gwo = ho*Sc; Gwi = hi*Si;           %convection wall out; wall in
Ggo = ho*Sg; Ggi = hi*Sg;           %convection glass out; glass in
% Long wave radiative exchange
GLW1 = epswLW/(1-epswLW)*Si*4*sigma*Tm^3;
GLW2 = Fwg*Si*4*sigma*Tm^3;
GLW3 = epsgLW/(1-epsgLW)*Sg*4*sigma*Tm^3;
GLW = 1/(1/GLW1 + 1/GLW2 +1/GLW3);  %long-wave exg. wall-glass
% ventilation & advection
Gv = Vpa*rhoa*ca;                   %air ventilation
% glass: convection outdoor & conduction
Ggs = 1/(1/Ggo + 1/(2*Gg));         %cv+cd glass

% Thermal network
% *****************************************************************
A(1,1) = 1;
A(2,1) = -1; A(2,2) = 1;
A(3,2) = -1; A(3,3) = 1;
A(4,3) = -1; A(4,4) = 1;
A(5,4) = -1; A(5,5) = 1;
A(6,5) = -1; A(6,6) = 1;
A(7,5) = -1; A(7,7) = 1;
A(8,6) = -1; A(8,7) = 1;
A(9,8) = 1;
A(10,6) = 1; A(10,8) = -1;
A(11,7) = 1;
A(12,7) = 1;

G = diag([Gwo 2*Gc 2*Gc 2*Gi 2*Gi GLW Gwi Ggi Ggs 2*Gg Gv Kp]');

C = diag([0 Cc 0 Ci 0 0 Ca Cg]);
C = diag([0 Cc 0 Ci 0 0 0   0]);

% Inputs
b = zeros(12,1); b([1 7 11 12]) = [10 70 110 120];
f = zeros(8,1); f([1 2 5 7 8]) = [1000 2000 5000 7000 8000];
y = ones(8,1);
u = nonzeros([b; f]);

% Thermal circuit -> state-space
[As,Bs,Cs,Ds] = fTC2SS(A,G,b,C,f,y);

yss = (-Cs*inv(As)*Bs + Ds)*u;
yade = inv(A'*G*A)*(A'*G*b + f);
[yss yade abs(yss-yade)<1e-10]
