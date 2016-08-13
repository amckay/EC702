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

u = @(C) C.^(1-Par.gamma)/(1-Par.gamma);

Bellman = @(Kp,K,Z,b) u( f(Par,K,Z) - Kp) ...
                     + Par.beta * PolyBasis(Kp,Z) * b;


MinKp = Grid.K(1) * ones(size(Grid.KK));                 
MaxKp = min(f(Par,Grid.KK,Grid.ZZ) - 1e-3, Grid.K(end)); % -1e-3 so we always have positve consumption.

p = (sqrt(5)-1)/2;

for it = 1:1000
    
    A = MinKp;
    D = MaxKp;
    
    MAXIT_INNER = 1000;
    for it_inner = 1:MAXIT_INNER
        B = p*A+(1-p)*D;
        C = (1-p)*A + p * D;

        fB = Bellman(B,Grid.KK,Grid.ZZ,b0);
        fC = Bellman(C,Grid.KK,Grid.ZZ,b0);

        I = fB > fC;

        D(I) = C(I);
        C(I) = B(I);
        fC(I) = fB(I);
        B(I) = p*C(I) + (1-p)*A(I);
        fB(I) = Bellman(B(I),Grid.KK(I),Grid.ZZ(I),b0);

        A(~I) = B(~I);
        B(~I) = C(~I);
        fB(~I) = fC(~I);
        C(~I) = p*B(~I) + (1-p)*D(~I);
        fC(~I) = Bellman(C(~I),Grid.KK(~I),Grid.ZZ(~I),b0);

        if all(D-A) < 1e-6
            break
        end
    
    end

    % At this stage, A, B, C, and D are all within a small epsilon of one
    % another.  We will use the average of B and C as the optimal level of
    % savings.
    Kp = (B+C)/2;
            
    % evaluate the Bellman equation at the optimal policy to find the new
    % value function.
    V = Bellman(Kp,Grid.KK,Grid.ZZ,b0);
    
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

bV =  PolyGetCoef(Grid.KK,Grid.ZZ,V);
V_K = PolyBasisDeriv(29,0.03) * bV;
C = (V_K / fprime(Par,29,0.03))^(-1/Par.gamma)
Kp = f(Par,29,0.03) - C