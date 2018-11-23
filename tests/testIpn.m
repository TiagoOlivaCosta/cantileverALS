clear all
p = 2;
n = 3;
A = magic(p);
I = eye(p);
% for i = 1:p
%     I = eye(p);
% end
T = eye(n);

for i = 1:n
    t(:,:,i) = T(:,i);   
end
for i =1:n
    
    II(:,:,i) = I;
end
B = krp(II,t);
Ipn = kron(B,eye(p));
