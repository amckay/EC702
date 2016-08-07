function b = PolyGetCoef(K,Z,Y)
%b = PolyGetCoef(Grid,Y)
%   Fits the polynomial from PolyBasis to the function(s) in column(s) of
%   Y.

X = PolyBasis(K,Z);

b = (X'* X) \ (X' * Y);


end

