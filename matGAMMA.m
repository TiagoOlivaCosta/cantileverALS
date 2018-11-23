function G = matGAMMA(A, C, N)
G = zeros(size(C,1)*N, size(C,2)*N);
O = OB(A,C,N);


d1 = size(C*A,1); %captura as dimensoes dos componentes da matriz de observabilidade
d2 = size(C*A,2);
o1 = size(O,1); %captura as dimensoes da matriz de observabilidade
o2 = size(O,2);

G = zeros(o1, (N-1)*o2+d2);

for i=1:N-1
    G(d1*(i+1):end, 1+(i-1)*d2:i*d2) = O(1:end-d1*i, :)

end
end
