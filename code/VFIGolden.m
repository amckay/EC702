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
b0 = zeros(6,1);


%% Bellman iteration

MAXIT = 2000;

for it = 1:MAXIT

    [V, Kp] = MaxBellman(Par,b0,Grid);

    
    % take the expectation of the value function from the perspective of
    % the previous Z
    EV = reshape(V,Grid.nK,Grid.nZ) * Grid.PZ; 
    
    % update our polynomial coefficients
    b1 = PolyGetCoef(Grid.KK,Grid.ZZ,EV(:));
    
    % see how much our coefficients have changed
    test = max(abs(b1 - b0));
    
    disp(['iteration ' num2str(it) ', test = ' num2str(test)])
    if test < 1e-5
        break
    end
    b0 = b1;
end


%%  Consumption and savings at K = 29, Z = 0.03
% we need to get C (or Kp)
% we know V
% from the Envelope condition we have V_K = u'(C) f_K(K,Z)
% from this we can do C = (V_K / f_K)^(-1/Par.gamma)
% how do we get V_K?
% V = sum_i b_i Fi(K,Z) where Fi are the basis functions
% so 
% V_K = sum_i b_i Fi_K(K,Z)

bKp =  PolyGetCoef(Grid.KK,Grid.ZZ,Kp);
Kp2903 = PolyBasis(29,0.03) * bKp
C2903 = f(Par,29,0.03) - Kp2903