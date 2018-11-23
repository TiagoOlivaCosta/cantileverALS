function S = ks(A, B)

S = kron(A, eye(length(B))) + kron(eye(length(A)), B );
end