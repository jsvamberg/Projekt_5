clear all;
clc;
close all;


% syms ddx ddphi dx dphi x phi M m g l f
% u = l*sin(phi)+x
% du = l*dphi*cos(phi)+dx
% 
% v = -l*cos(phi)
% dv = l*dphi*sin(phi)
% 
% T = (1/2)*(M*dx^2)+(1/2)*m*(du)^2+(1/2)*m*(dv)^2
% V = -m*g*l*cos(phi)
% L = T - V
% 
% 
% dL_ddx = diff(L, dx)
% dL_dx = diff(L, x)
% dL_ddphi = diff(L, dphi)
% dL_dphi = diff(L, phi)
% 
% %%
% D_dL_ddx = M*ddx+m*ddx+m*ddphi*l*cos(phi)-m*dphi*l*sin(phi)
% D_dL_ddphi = ddphi*l^2*m*(sin(phi))^2+dphi*l^2*m*cos(phi)*2*sin(phi)
% 
% eq1 = D_dL_ddx-dL_dx == f
% eq2 = D_dL_ddphi-dL_dphi == 0

%% Dosazení jedné rovnice za druhou
clear;
syms dx dphi x phi M m l f g ddphi ddx
 
eq1 = ddx == (-m*l*ddphi*cos(phi)+m*l*dphi^2*sin(phi)+f)/(M+m)
eq2 = ddphi == (sin(phi)*(dx*dphi-dx+g)+ddx*cos(phi))/(l)

[ddx, ddphi] = solve([eq1, eq2], [ddx, ddphi])

%% Parciální derivování
clear;
clc;
syms y1 y2 y3 y4 u m g M l 

f1 = y3
f2 = y4
f3 = (l*m*sin(y2)*y4^2-y3*m*cos(y2)*sin(y2)*y4+u+y3*m*cos(y2)*sin(y2)-g*m*cos(y2)*sin(y2))/(m*cos(y2)^2+m+M)
f4 = (u*cos(y2)-y3*m*sin(y2)+g*m*sin(y2)-M*y3*sin(y2)+M*g*sin(y2)+M*y4*y3*sin(y2)+y3*y4*m*sin(y2)+y4^2*l*m*cos(y2)*sin(y2))/(l*(m*cos(y2)^2+m+M))

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

% %%
% clear;
% close;
% u = 0; %sila pusobici na vozik
% M = 15; %hmostnost voziku
% m = 5; %hmotnost tělesa na lane
% y3 = 0; %pocatecni rychlost voziku
% y1 = 0; %pocatecni poloh voziku
% y4 = 0; %pocatecni rychlost kyvadla
% y2 = 0*pi/180; %pocatecni poloha kyvadla ve °
% l = 1; %delka zavesu
% g = 9.81; %gravitacni sila
% 
% A = [0, 0, 1, 0;
%     0, 0, 0, 1;
%     0, (g*m*sin(y2)^2 - g*m*cos(y2)^2 + m*y3*cos(y2)^2 - m*y3*sin(y2)^2 + m*y3*y4*sin(y2)^2 + l*m*y4^2*cos(y2) - m*y3*y4*cos(y2)^2)/(m*cos(y2)^2 + M + m) + (2*m*cos(y2)*sin(y2)*(l*m*sin(y2)*y4^2 - m*y3*cos(y2)*sin(y2)*y4 + u + m*y3*cos(y2)*sin(y2) - g*m*cos(y2)*sin(y2)))/(m*cos(y2)^2 + M + m)^2, (m*cos(y2)*sin(y2) - m*y4*cos(y2)*sin(y2))/(m*cos(y2)^2 + M + m), -(m*y3*cos(y2)*sin(y2) - 2*l*m*y4*sin(y2))/(m*cos(y2)^2 + M + m);
%     0, (g*m*cos(y2) - M*y3*cos(y2) - u*sin(y2) - m*y3*cos(y2) + M*g*cos(y2) + l*m*y4^2*cos(y2)^2 - l*m*y4^2*sin(y2)^2 + M*y3*y4*cos(y2) + m*y3*y4*cos(y2))/(l*(m*cos(y2)^2 + M + m)) + (2*m*cos(y2)*sin(y2)*(u*cos(y2) - M*y3*sin(y2) + g*m*sin(y2) - m*y3*sin(y2) + M*g*sin(y2) + M*y3*y4*sin(y2) + m*y3*y4*sin(y2) + l*m*y4^2*cos(y2)*sin(y2)))/(l*(m*cos(y2)^2 + M + m)^2), -(M*sin(y2) + m*sin(y2) - M*y4*sin(y2) - m*y4*sin(y2))/(l*(m*cos(y2)^2 + M + m)), (M*y3*sin(y2) + m*y3*sin(y2) + 2*l*m*y4*cos(y2)*sin(y2))/(l*(m*cos(y2)^2 + M + m))]
% 
% B = [0;
%     0;
%     1/(m*cos(y2)^2 + M + m);
%     cos(y2)/(l*(m*cos(y2)^2 + M + m))]
% 
% C = [1 0 0 0;
%     0 1 0 0;
%     0 0 1 0;
%     0 0 0 1]
% 
% D = 0
% 
% sys = ss(A, B, C, D)
%% linearizace v bode [0, 0, 0, 0]
df1_dy1 = subs(df1_dy1, [y1 y2 y3 y4], [0 pi 0 0]);
df1_dy2 = subs(df1_dy2, [y1 y2 y3 y4], [0 pi 0 0]);
df1_dy3 = subs(df1_dy3, [y1 y2 y3 y4], [0 pi 0 0]);
df1_dy4 = subs(df1_dy4, [y1 y2 y3 y4], [0 pi 0 0]);

df2_dy1 = subs(df2_dy1, [y1 y2 y3 y4], [0 pi 0 0]);
df2_dy2 = subs(df2_dy2, [y1 y2 y3 y4], [0 pi 0 0]);
df2_dy3 = subs(df2_dy3, [y1 y2 y3 y4], [0 pi 0 0]);
df2_dy4 = subs(df2_dy4, [y1 y2 y3 y4], [0 pi 0 0]);

df3_dy1 = subs(df3_dy1, [y1 y2 y3 y4], [0 pi 0 0]);
df3_dy2 = subs(df3_dy2, [y1 y2 y3 y4], [0 pi 0 0]);
df3_dy3 = subs(df3_dy3, [y1 y2 y3 y4], [0 pi 0 0]);
df3_dy4 = subs(df3_dy4, [y1 y2 y3 y4], [0 pi 0 0]);

df4_dy1 = subs(df4_dy1, [y1 y2 y3 y4], [0 pi 0 0]);
df4_dy2 = subs(df4_dy2, [y1 y2 y3 y4], [0 pi 0 0]);
df4_dy3 = subs(df4_dy3, [y1 y2 y3 y4], [0 pi 0 0]);
df4_dy4 = subs(df4_dy4, [y1 y2 y3 y4], [0 pi 0 0]);

df1_du = subs(df1_du, [y1 y2 y3 y4], [0 pi 0 0]);
df2_du = subs(df2_du, [y1 y2 y3 y4], [0 pi 0 0]);
df3_du = subs(df3_du, [y1 y2 y3 y4], [0 pi 0 0]);
df4_du = subs(df4_du, [y1 y2 y3 y4], [0 pi 0 0]);

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
phi = 5*pi/180; %pocatecni poloha kyvadla ve °
l = 1; %delka zavesu
g = 9.81; %gravitacni sila


A = [ 0, 0, 1, 0;
      0, 0, 0, 1;
      0, -(g*m)/(M + 2*m), 0, 0;
      0, -(M*g + g*m)/(l*(M + 2*m)), 0, 0]
B = [0;
     0;
     1/(M + 2*m);
     -1/(l*(M + 2*m))]
C = [1 0 0 0;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1]
D = [0;
    0;
    0;
    0]

sys = ss(A, B, C, D)

% Riditelnost
rank(ctrb(sys.A, sys.B)) %system je riditelny, protoze ma plnou radkovou hodnost


tf = tf(ss(A, B, C, D))

%% Navrh stavoveho regulatoru

p1 = -10 + 0i
p2 = -10 + 0i
p3 = -10 + 2.8014i
p4 = -10 - 2.8014i

K = acker(sys.A, sys.B, [p1, p2, p3, p4])
sys_cl = ss(sys.A-sys.B*K,sys.B,sys.C,sys.D)
%%
out=sim('Model_kyvadla_na_voziku_m_slx.slx',50)

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

figure;
subplot(4,1,1)
plot(out.stav_x);
legend('x')
title('Linearni model se stav. reg. - Poloha vozíku')
xlabel('t [s]')
ylabel('s [m]')
subplot(4,1,2)
plot(out.stav_phi);
legend('\phi')
title('Linearni model se stav. reg. - Uhel kyvadla \phi')
xlabel('t [s]')
ylabel('\phi [°]')
subplot(4,1,3)
plot(out.stav_dx);
legend('dx')
title('Linearni model se stav. reg. - Rychlost voziku')
xlabel('t [s]')
ylabel('v [m/s]')
subplot(4,1,4)
plot(out.stav_dphi);
legend('d\phi')
title('Linearni model se stav. reg. - Uhlova rychlost kyvadla')
xlabel('t [s]')
ylabel('\omega [rad/s]')