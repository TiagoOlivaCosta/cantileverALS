%% TESTE ESTATICO REALIZADO COM STRAIN GAUGES NA VIGA
%% O codigo le a medicao para ambos os strain gauges
%% e calcula a covariancia associada atraves do metodo ALS

close all
dt = 2.4414e-04; % [S]
n = 16384;
T = (n-1)*dt; % [s] Per?odo da janela de tempo
% nptf = 4097;
% f_Ny = delta_f*(nptf-1);
% f_hz = 0:delta_f:f_Ny;
% omega = 2*pi*f_hz;
%F_s = 2*f_Ny; % Sampling frequency
t = 0:dt:T;
path(path,'.')
path(path,'.\dados')
load strain_S1_X
S1 = Signal.y_values.values;
load strain_S2_X
S2 = Signal.y_values.values;

Q = diag([0.1, 0.1, 0.1, 0.1])*10^-12;
R = diag([223, 6])*10^-20;
P = eye(4)*0.1;

u = [0, 0]';
x1 = [0, 0, 0, 0]';
x2 = x1;
N = 10; %number of delays
run syst.m

for i=1:T % percorre vetor de medicao
  z1 = S1(i);                 % medicoes
  z2 = S2(i);
  [r1, y1] = inov(sys,x1,z1); % residuos de inovacao
  [r2, y2] = inov(sys,x2,z2);
  x1 = sysUpdate(sys,x1,u);   % propagacao de estados
  x2 = sysUpdate(sys,x2,u);



end

figure
subplot(2,1,1)
plot(t,S1,'k')
hold on
subplot(2,1,2)
plot(t,S2,'r')
hold off
xlabel('Time [s]')
ylabel('Strain 1')

figure
subplot(2,1,1)
plot(t(100:end),S1(100:end),'k')
hold on
subplot(2,1,2)
plot(t(100:end),S2(100:end),'k')
hold off
xlabel('Time [s]')
ylabel('Strain 2')
