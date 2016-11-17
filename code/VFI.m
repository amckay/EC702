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
clear all;
setup;

%% Initial guess of value function -> all zeros

%Note: we will approximate E[V(K',Z') | Z] so the unknown function is a
%function of K' and Z.
b = zeros(6,1);


%% Bellman iteration

Kp0 = zeros(size(Grid.KK));
MAXIT = 2000;
for it = 1:MAXIT

    [V, Kp] = MaxBellman(Par,b,Grid);

    
    % take the expectation of the value function from the perspective of
    % the previous Z
    EV = reshape(V,Grid.nK,Grid.nZ) * Grid.PZ; 
    
    % update our polynomial coefficients
    b = PolyGetCoef(Grid.KK,Grid.ZZ,EV(:));
    
    % see how much our policy rule has changed
    test = max(abs(Kp0 - Kp));
    Kp0 = Kp;
    disp(['iteration ' num2str(it) ', test = ' num2str(test)])
    if test < 1e-5
        break
    end
end

%% Make plot of policy rule


PolicyRulePlot;

%%  Consumption and savings at K = 29, Z = 0.03

bKp =  PolyGetCoef(Grid.KK,Grid.ZZ,Kp);
Kp2903 = PolyBasis(29,0.03) * bKp
C2903 = f(Par,29,0.03) - Kp2903


%% Accuracy plot
figure
bC = PolyGetCoef(Grid.KK,Grid.ZZ,f(Par,Grid.KK,Grid.ZZ)-Kp);
Accuracy(Par,Grid,bC);

%% Simulate
SimulateScript;

%% Impulse response functions
IRFScript;