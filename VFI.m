% Solve the model using value function iteration.
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


% create a structure that holds our constant parameters
Par.beta = 0.99;
Par.gamma = 2;
Par.alpha = 0.36;
Par.delta = 0.03;
Par.rho = 0.95;
Par.sigma = 0.007;


%% Solve for the steady state
Kstar = SteadyState(Par,0);

%% Create a grid for Z
meanZ = 0;
stdZ = Par.sigma / sqrt(1-Par.rho^2);
Grid.nZ = 7;
numStdZ = 2;
[Grid.Z, Grid.PZ]  = tauchen(N, mu, rho, sigma, m);  
    
Grid.PZ = Grid.PZ'; % this is a 7 x 7 transition matrix for which the columns sum to 1
% the (i,j) element is the probability of moving from j to i.

%% Create a grid for K
Grid.nK = 20;
Grid.K = linspace(SteadyState(Par,Grid.Z(1)), SteadyState(Par,Grid.Z(end)),Grid.nK)';  % this is a 20 x 1 array

%% Create a product of the two grids, the basis functions, and invert them
Grid.KZ = [kron(ones(Grid.nZ,1),Grid.K)  kron(Grid.Z,ones(Grid.nK,1))];
Grid.BasisFuncs = PolyBasis(Grid);
Grid.BasisFuncsInv = inv(Grid.BasisFuncs' * Grid.BasisFuncs);


%% Initial guess of value function -> all zeros

%Note: we will approximate E[V(K',Z') | Z] so the unknown function is a
%function of K' and Z.
b0 = zeros(6,1);


%% Bellman iteration
Kp = zeros(Grid.nK,Grid.nZ);  % Array to hold the savings decisions.

u = @(C) C.^(1-Par.gamma)/(1-Par.gamma);
f = @(K,Z) exp(Z) .* K.^Par.alpha + (1-Par.delta)*K;
Bellman = @(K,Z,Kp,b) u(f(K,Z) - Kp) + Par.beta * PolyBasis(Kp,Z) * b;

% for iK = 1:Grid.nK
%     for iZ = 1:Grid.nZ
%         Kp(iK,iZ) = fminbnd(