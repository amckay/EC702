function MU = MargUtil(Par,K,Z,b)
% MU = MargUtil(Par,K,Z,b)
%   Given capital K, TFP Z, and policy rule coefficients b, find u'(c)
%   where c = f(K,Z) - Kp

MU =  PolyBasis(K,Z) * b ;


end

