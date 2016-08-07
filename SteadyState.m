function [Kstar] = SteadyState(Par,Zstar)
  % Solve for the steady state
  % inputs
  % Par -- parameter structure
  % Zstar -- steady state log productivity
  % outputs
  % Kstar -- steady state capital stock
  

  % f'(Kstar) beta = 1
  % alpha exp(Zstar) Kstar^(alpha-1) + 1 - delta = 1/beta
  % alpha exp(Zstar) Kstar^(alpha-1)  = 1/beta - 1 + delta
  % Kstar^(alpha-1)  = (1/beta - 1 + delta)/ ( alpha exp(Zstar))
  Kstar = ((1/Par.beta - 1 + Par.delta)./(Par.alpha * exp(Zstar))).^(1/(Par.alpha-1));

end
