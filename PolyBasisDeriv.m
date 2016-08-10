function B  = PolyBasis(K,Z)
% B  = PolyBasis(K,Z)
% Derivative of polynomial basis functions.  Using 2nd order polynomial
%
% inputs
% K    n x 1   points for K
% Z    n x 1   points for Z
%     or scalar for Z
%
% outputs
% B    n x 6   array of basis functions: 1, K, Z, K^2, K*Z, Z^2

%B = [ones(size(K)) K Z K.^2 K.*Z Z.^2 K.^3 K.^2.*Z K.*Z.^2 Z.^3];
if numel(Z) ==1
    Zb = Z*ones(size(K));
else
    Zb = Z;
end
%B = [ones(size(K)) K Zb K.^2 K.*Zb Zb.^2];
B = [zeros(size(K)) ones(size(K)) zeros(size(K)) 2*K Zb zeros(size(K))];

end

