clear all;
clc;
close all;


syms ddx ddphi dx dphi x phi M m g l f
u = l*sin(phi)+x
du = l*dphi*cos(phi)+dx

v = -l*cos(phi)
dv = l*dphi*sin(phi)

T = (1/2)*M*dx^2+(1/2)*m*dx^2+(1/2)*m*du^2+(1/2)*m*dv^2
V = -m*g*l*cos(phi)
L = T - V


dL_ddx = diff(L, dx)
dL_dx = diff(L, x)
dL_ddphi = diff(L, dphi)
dL_dphi = diff(L, phi)

D_dL_ddx = M*ddx+2*m*ddx+m*l*ddphi*cos(phi)-m*l*dphi*sin(phi)
D_dL_ddphi = ddphi*l^2*m*sin(phi)^2+dphi*l^2*m*cos(phi)*2*sin(phi)+ddx*l*m*cos(phi)-dx*l*m*sin(phi)+m*l^2*ddphi*cos(phi)^2-m*l^2*dphi*sin(phi)*2*cos(phi)

eq1 = simplify(D_dL_ddx-dL_dx) == f
eq2 = simplify(D_dL_ddphi-dL_dphi) == 0
%% Dosazení jedné rovnice za druhou
clear;
syms ddphi_s ddx_s dx dphi x phi M m l f g

ddx_v = (-ddphi_s*l*m*cos(phi)-dphi*l*m*sin(phi)+f)/(M+2*m)
ddphi_v = (-ddx_s*l*m*cos(phi)+dx*l*m*sin(phi)-g*l*m*sin(phi)-l*m*dphi*dx*sin(phi))/(l^2*m)


ddx = subs(ddx_v, ddphi_s, ddphi_v)
ddphi = subs(ddphi_v, ddx_s, ddx_v)

ddx_simp = simplifyFraction(ddx)
ddphi_simp = simplifyFraction(ddphi)

ddx_simp - ddx
ddphi_simp - ddphi
% ddx_v = subs(ddx_v, ddphi, ddphi_v)
% ddphi_v = subs(ddphi_v, ddx, ddx_v)
% 
% ddx_v = (l*m*phi*dphi^2 + f - m*(ddx + g*phi))/(M + m)
%  
% ddphi_v = (g*phi + (l*m*phi*dphi^2 + f - m*(ddx + g*phi))/(M + m))/l

%% Parciální derivování
syms y1 y2 y3 y4 u m g M l 

f1 = y3
f2 = y4
f3 = (u-y4*l*m*sin(y2)-y3*m*cos(y2)*sin(y2)+g*m*cos(y2)*sin(y2)+y4*y3*m*cos(y2)*sin(y2))/(M+2*m-m*(cos(y2))^2)
f4 = -(u*cos(y2)-2*y3*m*sin(y2)+2*g*m*sin(y2)-M*y3*sin(y2)+M*g*sin(y2)+M*y4*y3*sin(y2)+2*y3*y4*m*sin(y2)-y4*l*m*cos(y2)*sin(y2))/(l*(M+2*m+m*(cos(y2))^2))

df1_dy1 = diff(f1, y1);
df1_dy2 = diff(f1, y2);
df1_dy3 = diff(f1, y3);
df1_dy4 = diff(f1, y4);

df2_dy1 = diff(f2, y1);
df2_dy2 = diff(f2, y2);
df2_dy3 = diff(f2, y3);
df2_dy4 = diff(f2, y4);

df3_dy1 = diff(f3, y1);
df3_dy2 = diff(f3, y2);
df3_dy3 = diff(f3, y3);
df3_dy4 = diff(f3, y4);

df4_dy1 = diff(f4, y1);
df4_dy2 = diff(f4, y2);
df4_dy3 = diff(f4, y3);
df4_dy4 = diff(f4, y4);

df1_du = diff(f1, u);
df2_du = diff(f2, u);
df3_du = diff(f3, u);
df4_du = diff(f4, u);

A = [df1_dy1 df1_dy2 df1_dy3 df1_dy4;
    df2_dy1 df2_dy2 df2_dy3 df2_dy4;
    df3_dy1 df3_dy2 df3_dy3 df3_dy4;
    df4_dy1 df4_dy2 df4_dy3 df4_dy4]

B = [df1_du;
    df2_du;
    df3_du;
    df4_du]

%% linearizace v bode [0, 0, 0, 0]
df1_dy1 = subs(df1_dy1, [y1 y2 y3 y4], [0 0 0 0]);
df1_dy2 = subs(df1_dy2, [y1 y2 y3 y4], [0 0 0 0]);
df1_dy3 = subs(df1_dy3, [y1 y2 y3 y4], [0 0 0 0]);
df1_dy4 = subs(df1_dy4, [y1 y2 y3 y4], [0 0 0 0]);

df2_dy1 = subs(df2_dy1, [y1 y2 y3 y4], [0 0 0 0]);
df2_dy2 = subs(df2_dy2, [y1 y2 y3 y4], [0 0 0 0]);
df2_dy3 = subs(df2_dy3, [y1 y2 y3 y4], [0 0 0 0]);
df2_dy4 = subs(df2_dy4, [y1 y2 y3 y4], [0 0 0 0]);

df3_dy1 = subs(df3_dy1, [y1 y2 y3 y4], [0 0 0 0]);
df3_dy2 = subs(df3_dy2, [y1 y2 y3 y4], [0 0 0 0]);
df3_dy3 = subs(df3_dy3, [y1 y2 y3 y4], [0 0 0 0]);
df3_dy4 = subs(df3_dy4, [y1 y2 y3 y4], [0 0 0 0]);

df4_dy1 = subs(df4_dy1, [y1 y2 y3 y4], [0 0 0 0]);
df4_dy2 = subs(df4_dy2, [y1 y2 y3 y4], [0 0 0 0]);
df4_dy3 = subs(df4_dy3, [y1 y2 y3 y4], [0 0 0 0]);
df4_dy4 = subs(df4_dy4, [y1 y2 y3 y4], [0 0 0 0]);

df1_du = subs(df1_du, [y1 y2 y3 y4], [0 0 0 0]);
df2_du = subs(df2_du, [y1 y2 y3 y4], [0 0 0 0]);
df3_du = subs(df3_du, [y1 y2 y3 y4], [0 0 0 0]);
df4_du = subs(df4_du, [y1 y2 y3 y4], [0 0 0 0]);

A = [df1_dy1 df1_dy2 df1_dy3 df1_dy4;
    df2_dy1 df2_dy2 df2_dy3 df2_dy4;
    df3_dy1 df3_dy2 df3_dy3 df3_dy4;
    df4_dy1 df4_dy2 df4_dy3 df4_dy4]

B = [df1_du;
    df2_du;
    df3_du;
    df4_du]

%% Konkretizace hodnot a zjisteni konkretniho modelu a prenosovych funkci
clear;
close;
f = 0; %sila pusobici na vozik
M = 15; %hmostnost voziku
m = 5; %hmotnost tělesa na lane
dx = 0; %pocatecni rychlost voziku
x = 0; %pocatecni poloh voziku
dphi = 0; %pocatecni rychlost kyvadla
phi = 0*pi/180; %pocatecni poloha kyvadla ve °
l = 1; %delka zavesu
g = 9.81; %gravitacni sila


A = [0 0 1 0;
    0 0 0 1;
    0 (g*m)/(M + m) 0 0;
    0 -(M*g + 2*g*m)/(l*(M + 3*m)) 0 0]
B = [0;
    0;
    1/(M + m);
    -1/(l*(M + 3*m))]
C = [1 0 0 0;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1]
D = 0

sys = ss(A, B, C, D)

tf = tf(ss(A, B, C, D))


out=sim('Model_kyvadla_na_voziku_m_slx.slx',20)

figure;
subplot(4,1,1)
plot(out.lin_x);
legend('x')
title('Linearizovany model - Poloha vozíku')
xlabel('t [s]')
ylabel('s [m]')
subplot(4,1,2)
plot(out.lin_phi);
legend('\phi')
title('Linearizovany model - Uhel kyvadla \phi')
xlabel('t [s]')
ylabel('\phi [°]')
subplot(4,1,3)
plot(out.lin_dx);
legend('dx')
title('Linearizovany model - Rychlost voziku')
xlabel('t [s]')
ylabel('v [m/s]')
subplot(4,1,4)
plot(out.lin_dphi);
legend('d\phi')
title('Linearizovany model - Uhlova rychlost kyvadla')
xlabel('t [s]')
ylabel('\omega [rad/s]')

figure;
subplot(4,1,1)
plot(out.nelin_x);
legend('x')
title('Nelinearizovany model - Poloha vozíku')
xlabel('t [s]')
ylabel('s [m]')
subplot(4,1,2)
plot(out.nelin_phi);
legend('\phi')
title('Nelinearizovany model - Uhel kyvadla \phi')
xlabel('t [s]')
ylabel('\phi [°]')
subplot(4,1,3)
plot(out.nelin_dx);
legend('dx')
title('Nelinearizovany model - Rychlost voziku')
xlabel('t [s]')
ylabel('v [m/s]')
subplot(4,1,4)
plot(out.nelin_dphi);
legend('d\phi')
title('Nelinearizovany model - Uhlova rychlost kyvadla')
xlabel('t [s]')
ylabel('\omega [rad/s]')