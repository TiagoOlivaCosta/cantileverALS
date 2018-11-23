% Chemical reactor example
% Odelson, Lutz and Rawlings
% IEEE Transactions on Control Systems Technology, May 2006

clear all
%% ************* reactor system  **************************

% variables
syms caf cas ccs... % [mol] concentration
     E  ... % [J/mol]
     R  ... % [J/mol/K]
     k0 ... % [l/s]
     V  ... % [l]
     T  ... % [K]
     Fw Fa ... % [l/s] Flow rates

% numerical values
    caf = 10.58; cas = 0.034236; ccs = 0.04;
    E = 4.544*10^4;
    R = 8.314;
    k0 = 2.6*10^5;
    V = 0.45;
    T = 296;
    Fw = 1.7611*10^-3;
    Fa = 9.0745*10^-6;
    %comment
    
a = (Fw+Fa)/V;
b = k0*exp(-E/(R*T));
A = [[-(a + b), 0];...
     [2*b, -a]];
B = [(caf-cas)/V; -ccs/V];
C = [0 1];
size_d = 1; % number of disturbance inputs;

%%% Augmented State-Space model to account for disturbances
   Bd = B;   Cd = zeros(size_d); % pure input disturbance assumption
%  Bd = zeros(size(A,1),1); Cd = eye(size_d);                 % pure output disturbance assumption
Gd = eye(size_d);

 Aa = [[A, Bd]; [zeros(1, size(A,2)), eye(size_d)]];
% Aa = [[A, Bd]; [zeros(size(A)), eye(size_d)]];
Ba = [B; zeros(1, size(B,2))];
Ca = [C, Cd];
% Ca = eye(3);
G = B; % hypothesis: the process disturbances actuate over
Ga = [G; Gd];
L = [-98.7857, -2.36762, 5.75749]';
% L = [-0.57, -1.79, 1.86]';

L = [-0.57857, -1.79762, 1.86749]';
A = Aa;
B = Ba;
C = Ca;
G = Ga;
G = eye(3);

N = 5;
[p,n] = size(C);
% **********************************************************

%% ALS matrices setup
% Initial matrices setup
% Abar = A-A*L*C;
Abar = A-A*L*C;
Abar2 = A - L*C;
Gbar = [G, -L];
O = OB(Abar,C,N);
Iota = I(A,C,N);
ds = dirSum(-A*L,N);
PS = Iota*ds;
In2 = eye(n^2);

% ALS procedure

D = kron(O,O)*(In2-kron(Abar,Abar))^-1 + kron(Iota, Iota)*iPn(n,N);
Als = [D*kron(G,G), D*kron(A*L, A*L)+(ks(PS,PS)+eye(p^2*N^2))*iPn(p,N)];

Aps = (A'*A)^-1*A';

eigs(A-Abar)
eigs(A-Abar2)
