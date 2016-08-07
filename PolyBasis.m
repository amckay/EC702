function B  = PolyBasis(K,Z)
% B  = PolyBasis(K,Z)
% Polynomial basis functions.  Using 2nd order polynomial
%
% inputs
% K    n x 1   points for K
% Z    n x 1   points for Z
%
% outputs
% B    n x 6   array of basis functions: 1, K, Z, K^2, K*Z, Z^2

B = [ones(size(K)) K Z K.^2 K.*Z Z.^2 K.^3 K.^2.*Z K.*Z.^2 Z.^3];


end

