function y  = fprime(Par, K,Z )
% y  = f( K,Z )
%   Derivative of production function gross of undepreciated capital

y =  Par.alpha * exp(Z) .* K.^(Par.alpha-1) + (1-Par.delta);
end

