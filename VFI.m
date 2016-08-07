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
Grid.K = linspace(0.9*Kstar, 1.1*Kstar,Grid.nK)';  % this is a 20 x 1 array

%% Create a product of the two grids, the basis functions, and invert them
Grid.KZ = [kron(ones(Grid.nZ,1),Grid.K)  kron(Grid.Z,ones(Grid.nK,1))];


%% Initial guess of value function -> all zeros

%Note: we will approximate E[V(K',Z') | Z] so the unknown function is a
%function of K' and Z.
EV = zeros(Grid.nK,Grid.nZ);


%% Bellman iteration
Kp = zeros(Grid.nK,Grid.nZ);  % Array to hold the savings decisions.
V = zeros(Grid.nK,Grid.nZ);   % temporary Array to hold the values

u = @(C) C.^(1-Par.gamma)/(1-Par.gamma);
f = @(K,Z) exp(Z) .* K.^Par.alpha + (1-Par.delta)*K;
Bellman = @(Kp,K,Z,EV) u( f(K,Z) - Kp) ...
                     + Par.beta * interp1(Grid.K,EV,Kp);
NegBellman = @(Kp,K,Z,EV) -Bellman(Kp,K,Z,EV);

for it = 1:1000
    for iZ = 1:Grid.nZ    
        for iK = 1:Grid.nK
            Kp(iK,iZ) = fminbnd(NegBellman,Grid.K(1),Grid.K(2),optimset('TolX',1e-12),...
                Grid.K(iK),Grid.Z(iZ),EV(:,iZ));
        end
        V(:,iZ) = Bellman(Kp(:,iZ),Grid.K,repmat(Grid.Z(iZ),Grid.nK,1),EV(:,iZ));
    end

    % take the expectation of the value function from the perspective of
    % the previous Z
    V = V * Grid.PZ; 

    % see how much our coefficients have changed
    test = max(abs(V(:) - EV(:)));
    
    disp(['iteration ' num2str(it) ', test = ' num2str(test)])
    if test < 1e-5
        break
    end
    EV = V;
end
