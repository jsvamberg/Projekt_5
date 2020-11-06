clear all;
clc;
close all;

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

%%
syms x1 x2 x3 x4 u m g M l 

f1 = x3
f2 = x4
f3 = (-m*g*sin(x2)*cos(x2)+m*x4^2*sin(x2)+u)/(M+m+m*(cos(x2))^2)
f4 = (m*l*x4^2*sin(x2)*cos(x2)+u*cos(x2)+m*g*sin(x2)+M*g*sin(x2))/(l*(m+M+m*(cos(x2))^2))

df1_dx1 = diff(f1, x1);
df1_dx2 = diff(f1, x2);
df1_dx3 = diff(f1, x3);
df1_dx4 = diff(f1, x4);

df2_dx1 = diff(f2, x1);
df2_dx2 = diff(f2, x2);
df2_dx3 = diff(f2, x3);
df2_dx4 = diff(f2, x4);

df3_dx1 = diff(f3, x1);
df3_dx2 = diff(f3, x2);
df3_dx3 = diff(f3, x3);
df3_dx4 = diff(f3, x4);

df4_dx1 = diff(f4, x1);
df4_dx2 = diff(f4, x2);
df4_dx3 = diff(f4, x3);
df4_dx4 = diff(f4, x4);

df1_du = diff(f1, u);
df2_du = diff(f2, u);
df3_du = diff(f3, u);
df4_du = diff(f4, u);

A = [df1_dx1 df1_dx2 df1_dx3 df1_dx4;
    df2_dx1 df2_dx2 df2_dx3 df2_dx4;
    df3_dx1 df3_dx2 df3_dx3 df3_dx4;
    df4_dx1 df4_dx2 df4_dx3 df4_dx4]

B = [df1_du;
    df2_du;
    df3_du;
    df4_du]

%%
df1_dx1 = subs(df1_dx1, [x1 x2 x3 x4], [0 0 0 0]);
df1_dx2 = subs(df1_dx2, [x1 x2 x3 x4], [0 0 0 0]);
df1_dx3 = subs(df1_dx3, [x1 x2 x3 x4], [0 0 0 0]);
df1_dx4 = subs(df1_dx4, [x1 x2 x3 x4], [0 0 0 0]);

df2_dx1 = subs(df2_dx1, [x1 x2 x3 x4], [0 0 0 0]);
df2_dx2 = subs(df2_dx2, [x1 x2 x3 x4], [0 0 0 0]);
df2_dx3 = subs(df2_dx3, [x1 x2 x3 x4], [0 0 0 0]);
df2_dx4 = subs(df2_dx4, [x1 x2 x3 x4], [0 0 0 0]);

df3_dx1 = subs(df3_dx1, [x1 x2 x3 x4], [0 0 0 0]);
df3_dx2 = subs(df3_dx2, [x1 x2 x3 x4], [0 0 0 0]);
df3_dx3 = subs(df3_dx3, [x1 x2 x3 x4], [0 0 0 0]);
df3_dx4 = subs(df3_dx4, [x1 x2 x3 x4], [0 0 0 0]);

df4_dx1 = subs(df4_dx1, [x1 x2 x3 x4], [0 0 0 0]);
df4_dx2 = subs(df4_dx2, [x1 x2 x3 x4], [0 0 0 0]);
df4_dx3 = subs(df4_dx3, [x1 x2 x3 x4], [0 0 0 0]);
df4_dx4 = subs(df4_dx4, [x1 x2 x3 x4], [0 0 0 0]);

df1_du = subs(df1_du, [x1 x2 x3 x4], [0 0 0 0]);
df2_du = subs(df2_du, [x1 x2 x3 x4], [0 0 0 0]);
df3_du = subs(df3_du, [x1 x2 x3 x4], [0 0 0 0]);
df4_du = subs(df4_du, [x1 x2 x3 x4], [0 0 0 0]);

A = [df1_dx1 df1_dx2 df1_dx3 df1_dx4;
    df2_dx1 df2_dx2 df2_dx3 df2_dx4;
    df3_dx1 df3_dx2 df3_dx3 df3_dx4;
    df4_dx1 df4_dx2 df4_dx3 df4_dx4]

B = [df1_du;
    df2_du;
    df3_du;
    df4_du]
%%
f = 0; %sila pusobici na vozik
M = 50; %hmostnost voziku
m = 5; %hmotnost tělesa na lane
dx = 0; %pocatecni rychlost voziku
x = 0; %pocatecni poloh voziku
dphi = 0; %pocatecni rychlost kyvadla
phi = 0; %pocatecni poloha kyvadla ve °
l = 1; %delka zavesu
g = 9.81; %gravitacni sila

A = [0 0 1 0;
    0 0 0 1;
    0 -(g*m)/(M + 2*m) 0 0;
    0 (M*g + g*m)/(l*(M + 2*m)) 0 0]
B = [0;
    0;
    1/(M + 2*m);
    1/(l*(M + 2*m))]
C = [1 0 0 0;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1]
D = 0

sys = ss(A, B, C, D)

tf(ss(A, B, C, D))
%step(sys)