% Solve the model using iteration on the Euler equation
%
% Description of the model.
% V(K,Z) = max_{C,K'} u(C) + beta E_{Z'} V(K',Z')
% C + K' = e^Z K^alpha + (1-delta) K
% Z follows a stochastic process
% Z' = rho  Z + epsilon' with epsilon ~  N(0,sigma^2)
% u(C) = C^{1-\gamma}/(1-gamma)
%
% Parameter values
% beta = 0.99;
% gamma = 2;
% alpha = 0.36;
% delta = 0.03;
% rho = 0.95;
% sigma = 0.007;
%
% Alisdair McKay, 8/5/16

setup;

%% Initial guess of policy rule for consumption

C = f(Par, Grid.KK,Grid.ZZ ) - Par.delta*Grid.KK;
b = PolyGetCoef(Grid.KK,Grid.ZZ,C);

%% Iteration on the Euler equation




MpOfZ = zeros(Grid.nK*Grid.nZ,Grid.nZ);

Kp0 = zeros(size(Grid.KK));

for it = 1:1000

    Kp = f(Par,Grid.KK,Grid.ZZ) - (PolyBasis(Grid.KK,Grid.ZZ)*b);
    
    % RHS of Euler equation
    for iZp = 1:Grid.nZ
        MpOfZ(:,iZp) = Par.beta* (PolyBasis(Kp,Grid.Z(iZp)) * b).^(-Par.gamma) .* fprime(Par,Kp,Grid.Z(iZp));
    end
    
    MU = sum( kron(Grid.PZ',ones(Grid.nK,1)) .* MpOfZ ,2);
    C = (MU).^(-1/Par.gamma);
    
    
    %of course, RHS = LHS
    b = PolyGetCoef(Grid.KK,Grid.ZZ,C);
    
    % see how much our policy rule has changed
    test = max(abs(Kp - Kp0));
    Kp0 = Kp;
    
    
    disp(['iteration ' num2str(it) ', test = ' num2str(test)])
    
    if test < 1e-5
        break
    end
    
    
end

%% Make plot of policy rule

PolicyRulePlot;

%%  Consumption and savings at K = 29, Z = 0.03
C2903 = PolyBasis(29,0.03) * b
Kp2903 = f(Par,29,0.03) - C2903