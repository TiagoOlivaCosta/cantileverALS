function [R1, R2, Q1, Q2, r] = covEstfun(bias_a, n_a, Z, n_c, t)
%acelerometro
v0 = 0;

aa = 0 + bias_a + n_a;
v_a = cumtrapz(t,aa)+v0;
x_a = cumtrapz(t, v_a);
%camera
f = 2;
u = f*0./Z + n_c;

du = [0;diff(u)];
yc = [u'; du'];

Q01 =  eye(2)*0.1;
R01 =  cov(x_a,v_a);
Q02 = eye(2)*0.1;
R02 =  cov(u,du);
x1 = zeros(2,1);     % Initial condition on the state
x2 = zeros(2,1);     % Initial condition on the state
P = Q01;         % Initial error covariance

x = [0;0];

R1 = R01;
R2 = R02;
Q1 = Q01;
Q2 = Q02;
S1 = R1;
S2 = R2;

va = 0;
t_step = 0.01;
for i = 1:length(t)-1
    va = va+aa(i)*t_step;
    xa = x1(1)+va*t_step;
    z1 = [xa, va]';
    z2 = yc(:,i);
     
    %estimacao das matrizes de covariancia
    window = 15;
    uw = ones(1,window)/window;
    r1 = x1-z1; %residuos
    r2 = x2-z2;
    if i<=window
        v1(:, i) = r1;
        v2(:, i) = r2;
    else
        v1 = circshift(v1,-1,2);
        v1(:,end) = r1;
        filtered1 = [conv(v1(1,:),uw,'same'); conv(v1(2,:),uw,'same')] ;
        
        v2 = circshift(v2,-1,2);
        v2(:,end) = r1;
        filtered2 = [conv(v2(1,:),uw,'same'); conv(v2(2,:),uw,'same')] ;
        
        r1 = filtered1(:,end);
        r2 = filtered2(:,end);
    end
    
%     alfa = 0.99;
%     r1 = alfa*r1old+(1-alfa)*r1;
%     r2 = alfa*r1old+(1-alfa)*r2;
%     r1old = r1;
%     r2old = r2;
%     r1 = filter(Bb, Ab, r1);
%     r2 = filter(Bb, Ab, r2);
     Ck1 = r1*r1';
     Ck2 = r2*r2';
     r(i) = r1(1);
     

     Cp1 = (diag(Ck1)-diag(S1));
     Cp2 = (diag(Ck2)-diag(S2));
    [a1(i), R1, Q1, M(:, i)] = covFuzzy(r1, R1, Q1, Cp1);
    [a2(i), R2, Q2, M2(:, i)] = covFuzzy(r2, R2, Q2, Cp2);
         
    [x1, P1] = KF(x, z1, P, Q1, R1);
    [x2, P2] = KF(x, z2, P, Q2, R2);
    
   
   
% CovInt
Paa = P1;
Pbb = P2;
I1 = eye(2)/Paa;
I2 = eye(2)/Pbb;

%[w, x, P] = CI(x1, x2, I1, I2, 'trace');
[w, x, P] = CI2(x1, x2, I1, I2);
 
end
end