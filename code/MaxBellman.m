function [V, Kp] = MaxBellman(Par,b,Grid)
% [V, Kp] = MaxBellman(Par,b,Grid)
%   Maximizes the RHS of the Bellman equation using golden section search
%
% Inputs
% Par       Parameter structure
% b     6 x 1 coefficients in polynomial for E[ V(K',Z') | Z ]
% Grid      Grid structure


p = (sqrt(5)-1)/2;

A = Grid.K(1) * ones(size(Grid.KK));
D = min(f(Par,Grid.KK,Grid.ZZ) - 1e-3, Grid.K(end)); % -1e-3 so we always have positve consumption.

B = p*A+(1-p)*D;
C = (1-p)*A + p * D;


fB = Bellman(Par,b,Grid.KK,Grid.ZZ,B);
fC = Bellman(Par,b,Grid.KK,Grid.ZZ,C);


MAXIT = 1000;
for it_inner = 1:MAXIT

    if all(D-A < 1e-6)
        break
    end
        
    
    I = fB > fC;
    
    D(I) = C(I);
    C(I) = B(I);
    fC(I) = fB(I);
    B(I) = p*C(I) + (1-p)*A(I);
    fB(I) = Bellman(Par,b,Grid.KK(I),Grid.ZZ(I),B(I));
    
    A(~I) = B(~I);
    B(~I) = C(~I);
    fB(~I) = fC(~I);
    C(~I) = p*B(~I) + (1-p)*D(~I);
    fC(~I) = Bellman(Par,b,Grid.KK(~I),Grid.ZZ(~I),C(~I));
    
end

% At this stage, A, B, C, and D are all within a small epsilon of one
% another.  We will use the average of B and C as the optimal level of
% savings.
Kp = (B+C)/2;

% evaluate the Bellman equation at the optimal policy to find the new
% value function.
V = Bellman(Par,b,Grid.KK,Grid.ZZ,Kp);

end