function y  = f(Par, K,Z )
% y  = f( K,Z )
%   Production function gross of undepreciated capital

y =  exp(Z) .* K.^Par.alpha + (1-Par.delta)*K;
end

