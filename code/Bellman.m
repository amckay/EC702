function V  = Bellman( Par, b, K, Z, Kp )
% V  = Bellman( Par, b, K, Z, Kp )
%   Evaluate the RHS of the Bellman equation
%
% Inputs
% Par   Parameter structure
% b     6 x 1 coefficients in polynomial for E[ V(K',Z') | Z ]
% K     n x 1 array of current capital
% Z     n x 1 array of current TFP
% Kp    n x 1 array of savings
%
% Output
% V     n x 1 array of value function
%

C = f(Par,K,Z) - Kp;
u = C.^(1-Par.gamma) / (1-Par.gamma);
V = u + Par.beta * PolyBasis(Kp,Z) * b;



end

