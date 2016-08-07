function B  = PolyBasis(Grid)
% B  = PolyBasis(K,Z)
% Polynomial basis functions.  Using 2nd order polynomial
%
% inputs
% Grid.KZ has two columsn
% K    n x 1   points for K
% Z    n x 1   points for Z
%
% outputs
% B    n x 6   array of basis functions: 1, K, Z, K^2, K*Z, Z^2

K = Grid.KZ(:,1);
Z = Grid.KZ(:,2);
B = [ones(size(K)) K Z K.^2 K.*Z Z.^2];


end

