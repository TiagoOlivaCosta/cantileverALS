function [AB] = krp(A,B)
% ***** KRB Khatri-Rao product *********
% A and B are four-dimensional matrices where the 3rd and 4th dimensions 
% define a matrix partition. i.e. each A(:, :, k,z) of the tensor 
% corresponds to a matrix.

[a1,a2,a3,a4] = size(A); % get matrices dimensions
[b1,b2,b3,b4] = size(B);
if (a3~=b3) ||  (a4~=b4)
   error(' Error in krp.m - Matrices dimensions are incompatible for Kahtri-Rao product')
end
S = zeros(a1*b1, a2*b2, a3*b3);
%AB = zeros(a1*b1*(a3*b3), a2*b2*(a3*b3));

for i = 1:a3
    for j = 1:a4
        P(:,:) = kron(A(:,:,i,j),B(:,:,i,j));
        %S(:,:, i+j-1) = P;
        AB(1+(i-1)*size(P,1):i*size(P,1), 1+(j-1)*size(P,2):j*size(P,2)) = P;
    end
end

% C = permute(S,[1 3 2]);
% AB = reshape(C,[],size(A,2),1);

end