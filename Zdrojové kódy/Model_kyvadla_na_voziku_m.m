clear all;
clc;
close all;

f = 0; %sila pusobici na vozik
M = 50; %hmostnost voziku
m = 5; %hmotnost tělesa na lane
dx = 0; %pocatecni rychlost voziku
x = 0; %pocatecni poloh voziku
dphi = 0; %pocatecni rychlost kyvadla
phi = 45*pi/180; %pocatecni poloha kyvadla ve °
l = 1; %delka zavesu
g = 9.81; %gravitacni sila

%%
syms ddphi_s ddx_s dx dphi x phi M m l f g

ddx_v = (-m*l*ddphi_s+m*l*(dphi^2)*phi+f)/(m+M)
ddphi_v = (ddx_s+g*phi)/(l)

ddx = subs(ddx_v, ddphi_s, ddphi_v)
ddphi = subs(ddphi_v, ddx_s, ddx_v)

ddx = simplifyFraction(ddx)
ddphi = simplifyFraction(ddphi)
% ddx_v = subs(ddx_v, ddphi, ddphi_v)
% ddphi_v = subs(ddphi_v, ddx, ddx_v)
% 
% ddx_v = (l*m*phi*dphi^2 + f - m*(ddx + g*phi))/(M + m)
%  
% ddphi_v = (g*phi + (l*m*phi*dphi^2 + f - m*(ddx + g*phi))/(M + m))/l


