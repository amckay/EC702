% Solve the model using policy iteration (Howard improvement).
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

%% Initial guess of value function -> all zeros

%Note: we will approximate E[V(K',Z') | Z] so the unknown function is a
%function of K' and Z.
b0 = zeros(6,1);


%% Bellman iteration
Kp = zeros(Grid.nK,Grid.nZ);  % Array to hold the savings decisions.

u = @(C) C.^(1-Par.gamma)/(1-Par.gamma);

Bellman = @(Kp,K,Z,b) u( f(Par,K,Z) - Kp) ...
    + Par.beta * PolyBasis(Kp,Z) * b;
NegBellman = @(Kp,K,Z,b) -Bellman(Kp,K,Z,b);

InnerIsDone = false;
for itOuter = 1:1000
    for itInner = 1:1000
        if itInner == 1
            for iK = 1:Grid.nK
                for iZ = 1:Grid.nZ
                    % first find the point at which consumption is negative
                    % Kp = f(K);
                    MaxKp = f(Par,Grid.K(iK),Grid.Z(iZ)) - 1e-3; % -1e-3 so we always have positve consumption.

                    Kp(iK,iZ) = fminbnd(NegBellman,Grid.K(1),MaxKp,optimset('TolX',1e-12),...
                        Grid.K(iK),Grid.Z(iZ),b0);
                end
            end
        end
        
        % evaluate the Bellman equation at the optimal policy to find the new
        % value function.
        V = Bellman(Kp(:),Grid.KK,Grid.ZZ,b0);
        
        % take the expectation of the value function from the perspective of
        % the previous Z
        EV = reshape(V,Grid.nK,Grid.nZ) * Grid.PZ;
        
        % update our polynomial coefficients
        b1 = PolyGetCoef(Grid.KK,Grid.ZZ,EV(:));
        
        if itInner == 1
            % see how much our coefficients have changed
            test = max(abs(b1 - b0));

            disp(['iteration ' num2str(itOuter) ', test = ' num2str(test)])
            if test < 1e-5
                InnerIsDone = true;
                break
            end
        end
        
        b0 = b1;
    end
    if InnerIsDone
        break
    end
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