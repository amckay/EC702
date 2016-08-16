% Setup the model parameters and grids
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
Par.Kstar = Kstar;

%% Create a grid for Z
meanZ = 0;
Grid.nZ = 7;
numStdZ = 2;
[Grid.Z, Grid.PZ]  = tauchen(Grid.nZ, meanZ, Par.rho, Par.sigma, numStdZ);  
    
Grid.PZ = Grid.PZ'; % this is a 7 x 7 transition matrix for which the columns sum to 1
% the (i,j) element is the probability of moving from j to i.

%% Create a grid for K
Grid.nK = 20;
Grid.K = linspace(0.75*Kstar, 1.25*Kstar,Grid.nK)';  % this is a 20 x 1 array

%% Create a product of the two grids
[ZZ,KK] =meshgrid(Grid.Z,Grid.K);
Grid.KK = KK(:);
Grid.ZZ = ZZ(:);