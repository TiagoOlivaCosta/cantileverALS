% function ALS(
%B.J. Odelson et al
N = 3;

run syst.m
L = [38.5726    0.7684;
   -1.8059    0.2360;
  -48.4661   -0.9655;
   77.4313  -10.1172];
 G = B;
 Abar = A-A*L*C;
 Gbar = [G, -A*L];
 
 O = OB(Abar,C,N);
 Iota = I(A,C,N);
 ds = dirSum(-A*L,N);
 PS = Iota*ds;
 %As = [D*kron(G,G), D*kron(A*L,A*L)+PP
D = [kron(O,O)*(In2-kron(Abar,Abar))^-1 + kron(Iota, Iota)*perm];
    
 ACM = O*P*O' + Iota*(dirSum(Gbar*Qbar*Gbar',N))*Iota' + ...
      PS*dirSum(Rv,N)+dirSum(Rv,N)*PS'+dirSum(Rv,N);
  

    In2 = eye(length(Abar)^2);
    Rv = eye(6);
    
    D = [kron(O,O)*(In2-kron(Abar,Abar))^-1 + kron(Iota, Iota)*perm];
    As = [D*kron(G,G), D*kron(A*L, A*L)+(ks(PS,PS)+eye(p^2*N^2)]
    