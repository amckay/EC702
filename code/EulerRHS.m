function C = EulerRHS(Par,Grid,Kp,b)
% C = EulerRHS(Par,Kp,b,nK)
% RHS of Euler equation
% 
% inputs
% Par    Parameter structure
% Grid   Grid structure
% Kp     nK*Grid.nZ x 1  array of savings
% b      consumption function polynomial coefficients
%
% outputs
% C      consumption implied by the Euler equation



MpOfZ = zeros(Grid.nK*Grid.nZ,Grid.nZ);
for iZp = 1:Grid.nZ
    MpOfZ(:,iZp) = Par.beta* (PolyBasis(Kp,Grid.Z(iZp)) * b).^(-Par.gamma) .* fprime(Par,Kp,Grid.Z(iZp));
end

MU = sum( kron(Grid.PZ',ones(Grid.nK,1)) .* MpOfZ ,2);
C = (MU).^(-1/Par.gamma);