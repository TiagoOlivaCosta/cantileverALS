f1 = 38.457;
f2 = 258.488;
omega_1 = 2*pi*f1; % Natural frequency 1
omega_2 = 2*pi*f2; % Natural frequency 2

% FORCEI OS MODOS REAIS (Leopoldo) parte imagin�ria despres�vel
psi_1 = [9.812e-1; 2.673]; % Modo de deslocamento REAL #1
psi_2 = [-2.617;   1.682]; % Modo de deslocamento REAL #2

zeta_1 = 0.26/100;
zeta_2 = 1.32/100;
omega_r2 = [(omega_1^2),0;0,(omega_2^2)];
PHI = [psi_1,psi_2]; % Displacement Modal Matrix

% Matriz massa
M = inv(PHI.')*inv(PHI); % Maia & Silva 1998 pg. 83

% Matriz rigidez
K = inv(PHI.')*omega_r2*inv(PHI); % Maia & Silva 1998 pg. 83

% Damping matrix as a function of modal damping ratio

DA = [(2*zeta_1*omega_1),0;0,(2*zeta_2*omega_2)];

C = inv(PHI.')*DA*inv(PHI);

%% Strain Modal Matrix PHI_ep

C_i1_j2_m1 = 6.0130e4; % Residues
C_i1_j2_m2 = 1.2322e5;

C_i2_j2_m1 = 4.6020e5;
C_i2_j2_m2 = 6.1853e6;

psi_sup_eps_1 = [C_i1_j2_m1; C_i2_j2_m1]./PHI(2,1);
psi_sup_eps_2 = [-C_i1_j2_m2; C_i2_j2_m2]./PHI(2,2);

PSI = 1e-6*[psi_sup_eps_1, psi_sup_eps_2]; % Strain Modal Matrix

% Mass, damping and stiffness MODAL MATRICES for strain EMA

H1 = PSI*inv(PHI);

% Mass modal matrix strain
M_ep = inv(H1')*M*inv(H1);

% Damping modal matrix strain
C_ep = inv(H1')*C*inv(H1);

% Stiffness modal matrix strain
K_ep = inv(H1')*K*inv(H1);

% Mass modal matrix strain
M_m = PSI'*M_ep*PSI;

% Damping modal matrix strain
C_m = PSI'*C_ep*PSI;

% Stiffness modal matrix strain
K_m = PSI'*K_ep*PSI;

Z = [(zeta_1),0;0,(zeta_2)]; % Damping factors matrix
OMEGA = [omega_1, 0; 0, omega_2]; % Natural frequencies matrix
b = [[0,0];[0,1]]; % matriz booleana de aplica��o da forca

A = [[zeros(2), eye(2)]; [-OMEGA.^2, -2*Z*OMEGA]];
B = [zeros(2); PHI'*b];

%H_a = [[-PHI, zeros(2,2)];[zeros(2,2), -PHI]];
%H_a = [[-inv(M)*K, zeros(2,2)];[zeros(2,2), -inv(M)*C]];
C_ac = PHI*[-OMEGA.^2, -2*Z*OMEGA];
C_ss = [PSI, zeros(2,2)];

H = [C_ac; C_ss];
D = [PHI*PHI'*b];
C = C_ss;

sys = ss(A,B,C,D);
