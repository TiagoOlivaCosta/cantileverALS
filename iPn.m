function U = iPn(p, n)

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
U = kron(B,eye(p));

end