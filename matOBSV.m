function O = matOBSV(A,C, N)
    O = zeros(N*size(C,1), size(C,2));
for i = 1:N
    O(1+size(C,1)*(i-1):i*size(C,1), 1:size(C,2)) = C*A^(i-1);
end

end
