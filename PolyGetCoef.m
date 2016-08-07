function b = PolyGetCoef(Grid,Y,reg)
%b = PolyGetCoef(Grid,Y)
%   Fits the polynomial from PolyBasis to the function(s) in column(s) of
%   Y.

if ~exist('reg','var')
    reg = 0;
end

b = (Grid.BasisFuncsInv + reg*eye(size(Grid.BasisFuncsInv)))\(Grid.BasisFuncs'*Y);


end

