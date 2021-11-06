%autor: Zdenek Boucek

%nastaveni konstant
m=15; %hmotnost voziku
M=5; %hmotnost telesa na kyvadle
l1=1; %delka kyvadla
k=0; %tuhost pruzin
g=9.81; %gravitacni konstanta
r0=0; %pocatecni poloha voziku
theta0=deg2rad(5); %pocatecni odchyleni kyvadla

sim('sm_test_wout_simmech',100)
figure
subplot(2,1,1); %podgraf s polohou voziku
plot(simout.time, simout.signals.values(:,1));
title('vozík');
xlabel('t[s]');
ylabel('x[m]');
subplot(2,1,2); %podgraf s odchylkou kyvadla
plot(simout1.time, simout1.signals.values(:,1));
title('kyvadlo');
xlabel('t[s]');
ylabel('\Theta[rad]');