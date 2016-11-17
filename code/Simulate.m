function Sim = Simulate(Par,bKp,Mode,T)
% Sim = Simulate(Par,bKp,Mode)
% 
% Simulates the model.
%
% Inputs:
% Par       Parameter structure
% bKp       Polynomial coefficients for polynomial for Kp policy rule
% Mode      Mode = 'random' -> draw shocks
%           Mode = 'irf'    -> impulse response function



Sim.K = zeros(T,1);
Sim.Z = zeros(T,1);

Sim.K(1) = Par.Kstar;


if strcmp(Mode,'irf')
    Sim.Z(1) = Par.sigma;
    Eps = zeros(T,1);
elseif strcmp(Mode,'random')
    Sim.Z(1) = 0;
    Eps = Par.sigma * randn(T,1);
else
    error('Unrecognized Mode in Simulate.m');
end



for t = 2:T
    Sim.K(t) = PolyBasis(Sim.K(t-1),Sim.Z(t-1)) * bKp;
    Sim.Z(t) = Par.rho * Sim.Z(t-1) + Eps(t);
end


% Compute quantities from state variables
Ti = 2:T-1;
Kp = Sim.K(Ti+1);
Sim.K = Sim.K(Ti);
Sim.Z = Sim.Z(Ti);
Sim.Y = f(Par,Sim.K,Sim.Z) - (1-Par.delta) * Sim.K;
Sim.C = f(Par,Sim.K,Sim.Z) - Kp;

end
