% Solve the model using time iteration.
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

%Note we will approximate u'(C) and call it M
M = C.^(-Par.gamma);


b0 = PolyGetCoef(Grid.KK,Grid.ZZ,M);


%% Iteration on the Euler equation


X = PolyBasis(Grid.KK,Grid.ZZ);


MpOfZ = zeros(Grid.nK*Grid.nZ,Grid.nZ);
M = zeros(Grid.nK,Grid.nZ);

for it = 1:1000

    Kp = f(Par,Grid.KK,Grid.ZZ) - (X*b0).^(-1/Par.gamma);
    
    % RHS of Euler equation
    for iZp = 1:Grid.nZ
        MpOfZ(:,iZp) = Par.beta* MargUtil(Par,Kp,Grid.Z(iZp),b0) .* fprime(Par,Kp,Grid.Z(iZp));
    end
    
    
    for iZ = 1:Grid.nZ
        M(:,iZ) = MpOfZ( (1:Grid.nK)+ (iZ-1)*Grid.nK, : ) * Grid.PZ(:,iZ);
    end
    
    
    %of course, RHS = LHS
    b1 = PolyGetCoef(Grid.KK,Grid.ZZ,M(:));
    
    % see how much our coefficients have changed
    test = max(abs(b1 - b0));
    
    disp(['iteration ' num2str(it) ', test = ' num2str(test)])
    if test < 1e-5
        break
    end
    b0 = 0.75*b0 + 0.25*b1;
    
end


%%  Consumption and savings at K = 29, Z = 0.03
C = (PolyBasis(29,0.03) * b0)^(-1/Par.gamma)
Kp = f(Par,29,0.03) - C