function matI = I(A, C, N)
matI = zeros(size(C,1)*N, size(C,2)*N);

for i=1:N-1
    for i=1:N-1
        P = C*A^(i-1);
      matI(2+size(P,1)*(i-1):1+size(P,1)*(i-1)+size(P,1), 1+size(P,2)*(i-1):i*size(P,2)) = P;
    end
end
end
