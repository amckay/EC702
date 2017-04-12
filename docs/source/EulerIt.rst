

Iterating on the Euler equation
==================================

.. highlight:: matlab

We will now discsuss another method for solving the model.  There are two important reasons for considering this alternative. First, it is often more accurate to approximate the policy rules rather than the value function.  To understand this, suppose we had a model with a constant productivity so capital is the only state.  If we are interested in the optimal policy rules, we don't need to know the level of the value function as this doesn't change the behavior in the maximization.  So if we approximate the value function, we are wasting one of our polynomial coefficients on the level of the value function, which we are not interested in.  Instead we would do better to approximate the derivative of the value function.  But from the envelope condition, the derivative of the value function is a one-to-one function of the consumption function.  So if we approximate the consumption function we focus on what is important.  You will see that this method gives a slightly different answer from the previous method and then in the next chapter we will see that this method is more accurate.

The second reason to consider this alternative method is that it can easily be adapted to solve models with many frictions so that a competitive equilibrium is not equal to a planner's problem.  One can solve for competitive equilibrium of distorted economies using value function iteration, but typically one needs to impose the equilibrium conditions in an outer loop.  So one solves the decision problems of the agents using value function iteration then checks if their actions are compatitible with one another in the sense of equilirbium and if not one adjusts the prices and tries again.  Using the method presented here, we can compute an equilibrium and solve for the decision rules at the same time.   `This paper <http://onlinelibrary.wiley.com/doi/10.3982/QE364/abstract>`_ is a good example of this method applied to a distorted equilibrium.

Relative to value function iteration, the drawback of this method is that it will not necessarily converge so it is important here to have a decent initial guess of the solution.

.. It would be nice to come up with a simple example application that doesn't fit directly into VFI.  Maybe a labor income tax.

Here we are going to work with a set of functional equations that need to be satisfied by the policy rule that solves the model without including a maximization anywhere.  So instead of the Bellman equation we use the Euler equation.  I have repeated the equations here for your reference
 .. math::

   C(K,Z)^{-\gamma} &= \beta \mathbb E [ f_K(K',Z') {C(K',Z')}^{-\gamma}] \\
   K' &= f(K,Z) - C(K,Z) \\
   f(K,Z) &= e^Z K^\alpha + (1-\delta)K \\
   Z' &= \rho Z + \varepsilon'.

We would usually refer to such a system of equations as the "equilibrium conditions" of the model although in this context that term isn't entirely appropriate because we are looking a planner's problem and not a competitive equilibrium.


We start the script ``EulerIteration.m`` by defining the same parameters and grids as we did in ``VFI.m`` and I don't show those steps here.  We next create an initial guess for the consumption function.
::

  C = f(Par, Grid.KK,Grid.ZZ ) - Par.delta*Grid.KK;
  b0 = PolyGetCoef(Grid.KK,Grid.ZZ,C);

We guess that we consume output net of depreciation and then fit our second-order polynomial to that function.


Now let's create a function that takes a value for :math:`K_{t+1}` and computes consumption at date :math:`t` from the right-hand side of the Euler equation.
::

  function C = EulerRHS(Par,Grid,Kp,b)
  % C = EulerRHS(Par,Kp,b,nK)
  % RHS of Euler equation
  %
  % inputs
  % Par    Parameter structure
  % Grid   Grid structure
  % Kp     nK*Grid.nZ x 1  array of savings
  % b      consumption function polynomial coefficients
  %
  % outputs
  % C      consumption implied by the Euler equation

  MpOfZ = zeros(Grid.nK*Grid.nZ,Grid.nZ);
  for iZp = 1:Grid.nZ
      MpOfZ(:,iZp) = Par.beta* (PolyBasis(Kp,Grid.Z(iZp)) * b).^(-Par.gamma) .* fprime(Par,Kp,Grid.Z(iZp));
  end

  MU = sum( kron(Grid.PZ',ones(Grid.nK,1)) .* MpOfZ ,2);
  C = (MU).^(-1/Par.gamma);

The array ``MpOfZ`` will hold the value for :math:`\beta f_K(K',Z') {C(K',Z')}^{-\gamma}` that will be acheived conditional on a given outcome for :math:`Z'` for each grid point for :math:`K'`.  The rows of ``MpOfZ`` correspond to grid points and the columns to realizations of :math:`Z'`.

 Next, we loop over all future productivity levels
::

  for iZp = 1:Grid.nZ
      MpOfZ(:,iZp) = Par.beta* (PolyBasis(Kp,Grid.Z(iZp)) * b).^(-Par.gamma) .* fprime(Par,Kp,Grid.Z(iZp));
  end

and for each we calculate :math:`\beta f_K(K',Z') {C(K',Z')}^{-\gamma}` and store it in ``MpOfZ``.  The function ``f_prime`` is analogous to our function ``f`` but it evaluates the derivative of the gross production function :math:`\alpha e^Z K^{\alpha-1} + 1 -\delta`.  To compute marginal utility, we again use our approximated consumption function but now evaluated at the future state ``(Kp,Grid.Z(iZp))``.

The next step is to take the expectation over the outcomes for :math:`Z'`
::

   MU = sum( kron(Grid.PZ',ones(Grid.nK,1)) .* MpOfZ ,2);

Let's understand this step from the inside out.  ``kron`` is the `Kroenecker product <https://en.wikipedia.org/wiki/Kronecker_product>`_ and we use it here as a way of repeating the entries of our Markov chain transition matrix ``Grid.PZ``.  In particular, we transpose the transition matrix so now the columns correspond to future outcomes, and then we repeat each row ``Grid.nK`` or 20 times.  We do this because the first 20 rows of our grid are different levels of capital with the first productivity level.  The next 20 rows of our grid are different capital level with the second productivity level.  We take this result and multiply it element-by-element against ``MpOfZ`` by using ``.*``.  Notice that the columns ``MpOfZ`` correspond to different realizations of :math:`Z'` just like the columns of our transposed transition matrix. Finally we take the sum across each row to take the expectation.

Why couldn't we take an expectation with matrix multiplication like we did before?  In value function iteration, we had :math:`V(K',Z')` evaluated at the grid points of capital and productivity and then we used matrix multiplication to compute :math:`\mathbb E \left[ V(K',Z') | Z \right]`.  In this application we have a function of :math:`(K',Z')` but it is evaluated at a different :math:`K'` for each :math:`Z` because we are using the :math:`K'` that is implied by the consumption choice at :math:`(K,Z)` (see below). If we computed ``MpOfZ * Grid.PZ`` we would have an array of size 140 x 7, with the conditional expectations for each value of :math:`Z` in the columns, but we only really want the conditional expecation that is appropriate for the value of :math:`Z` that corresponds to that row of the grid.  We could start with ``MpOfZ * Grid.PZ`` and then look at the right column of the result, but doing that would be just as much work for us and more work for the computer because we would compute expectations that we do not need.

The rest is pretty straightforward.  ``MU`` is now the right-hand side of the Euler equation, which is equal to the left-hand side, which is marginal utility of consumption.  Inverting that we find consumption with ``C = (MU).^(-1/Par.gamma);``


The function ``EulerRHS`` takes a set of polynomial coefficients that approximate the consumption function in the next period and computes the level of consumption in the current period.  The algorithm for solving the model iterates on finding ``C`` from ``EulerRHS`` and then using that result to update the polynomial coefficients ``b``.

.. The main part of the algorithm works as follows: given a pair :math:`(K,Z)` we use our consumption function and the aggregate resource constraint to find :math:`K'`.  For that value of :math:`K'` and averaging over the possible realizations of :math:`Z'` we compute the right-hand side of the Euler equation using our consumption function again.  From the left-hand side of the Euler equation we can then determine a new value for :math:`C(K,Z)`.  We do this for all pairs  :math:`(K,Z)` in our grid and then fit a new polynomial to the values for :math:`C(K,Z)` that come from the left-hand side of the Euler equation.  We iterate on this procedure until the consumption function converges.

The following code implements the algorithm
::

  Kp0 = zeros(size(Grid.KK));

  for it = 1:1000

      Kp = f(Par,Grid.KK,Grid.ZZ) - (PolyBasis(Grid.KK,Grid.ZZ)*b);

      C = EulerRHS(Par,Grid,Kp,b);

      b = PolyGetCoef(Grid.KK,Grid.ZZ,C);

      % see how much our policy rule has changed
      test = max(abs(Kp - Kp0));
      Kp0 = Kp;


      disp(['iteration ' num2str(it) ', test = ' num2str(test)])

      if test < 1e-5
          break
      end

  end

There is just one part of this algorithm that is not self explanatory. When we call ``EulerRHS``, we need to know ``Kp``, but the optimal level of savings depends on the consumption function that we are trying to solve for so we don't know it yet.  Instead, we simply use our current approximate consumption function to calculate ``Kp``.

Run ``EulerIteration.m`` and see how it works.
