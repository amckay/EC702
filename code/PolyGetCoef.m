function b = PolyGetCoef(K,Z,Y)
% b = PolyGetCoef(Grid,Y)
%   Fits the polynomial from PolyBasis to the function(s) in column(s) of
%   Y.
%
% inputs
% K    n x 1   points for K
% Z    n x 1   points for Z
% Y    n x 1   valies for function at (K,Z)
%
% outputs
% b    6 x 1   basis coefficients

b = PolyBasis(K,Z) \ Y;


end

