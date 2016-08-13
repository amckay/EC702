.. NumericalAnalysis documentation master file, created by
   sphinx-quickstart on Thu Aug 11 20:18:56 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Numerical analysis notes for EC 702
=============================================

.. toctree::
   :maxdepth: 1

   Overview <index>
   Function approximation <FuncApprox>
   



Overview
=================================================

We started this part of the course asking what does it mean to solve a model?  Our answer was that solving a model means finding a function or stochastic process that satisfies some conditions.  But how do you do that?  As we have seen, in very special cases we can do this analytically but in most cases there is no closed-form solution.  Most macroeconomic models are analyzed numerically.  These notes present two approaches to computing numerical solutions of models. The first approach searches for a value function that satisfies the Bellman equation.  The second approach searches for a policy rules that satisfy the equilibrium conditions of the model (e.g. we solve the functional Euler equation).  The advantages of value function methods are that they are applicable to complex models and they have good convergence properties.  The advantage of working directly with the equilibrium conditions is that they are often faster. In the course of these notes we will cover some fundamental techniques of numerical analysis.  We will discuss how we can represent  function on the computer and we will discuss how we can take expectations on the computer.


In these notes, I will use the programming language Matlab to demonstrate the methods.  For the applications we will consider we could use any number of different languages, but graduate students typically start out with Matlab and so that is why I have chosen it.  I think Matlab's popularity stems in part from the fact that it has good documentation.  There are some disadvantages of Matlab that might lead you to consider alternatives.  First, you need to buy a license in order to install it on your computer while some of the alternatives are free open-source software.  Second, for some applications Matlab is slower than some of the alternatives.  You can find a thorough comparison of programming languages from the perspective of economists `here <http://www.nber.org/papers/w20263>`_.



First application: stochastic neoclassical growth model
========================================================

Let's consider an application to fix ideas.  Suppose we want to solve the following version of the stochastic neoclassical growth model
 .. math::

 		V(K,Z) &=  \max_{K',C} \left\{ \frac{C^{1-\gamma}}{1-\gamma} + \beta \mathbb E V(K',Z') \right\} \\
        K' &= f(K,Z) - C \\
        f(K,Z) &= e^Z K^\alpha + (1-\delta)K \\
        Z' &= \rho Z + \varepsilon'.

Our first approach to solving the model is to find a function :math:`V(K,Z)` that satisfies the Bellman equation.  We could alternatively express the model in terms of equilibrium conditions using the Euler equation
 .. math::

 		C(K,Z)^{-\gamma} &= \beta \mathbb E [ f_K(K',Z') {C(K',Z')}^{-\gamma}] \\
        K' &= f(K,Z) - C \\
        f(K,Z) &= e^Z K^\alpha + (1-\delta)K \\
        Z' &= \rho Z + \varepsilon'.

Our second approach to solving the model will be to find a policy rule for consumption that satisfies the Euler equation.

In order to solve the model numerically, we need values for the parameters.  We will discuss how we use data to choose parameters, but for now we will simply assume :math:`\gamma = 2`, :math:`\beta = 0.99`, :math:`\beta = 0.99`, :math:`\delta = 0.03`, and :math:`\rho = 0.95`.  We also need a distribution for the innovations and we will assume :math:`\varepsilon \sim N(0,\sigma^2)` where the standard deviation :math:`\sigma = 0.007`.


Before we get to the discussion of solving the model, we will need to discuss some fundamental techniques of numerical analysis.

Function approximation
========================================================

Suppose there is a function :math:`f(x)` defined on some domain :math:`\mathcal{X} \subset \mathbb R`.   People sometimes say that a function is "an infinite dimensional object" meaning that one could represent the function as an infinite table of values associated with each :math:`x \in X`.  Representing a function this way would require an infinite amount of memory on a computer so we use a different approach.  Instead we will represnt :math:`f(x)` as a weighted some of a finite number of pre-defined functions, which we call "basis functions."  For example we know that any continuous function can be well approximated by a polynomial so we could use polynomials as our basis functions and then look for the coefficients on the polynomials to approximate :math:`f(x)`.  Concretely, let :math:`B_i(x) = x^i` so we have :math:`B_0(x) = 1`, :math:`B_1(x) = x`, :math:`B_2(x) = x^2`, and so on.  Let's now suppose we will use a cubic polynomial to approximate :math:`f(x)`.  We then have
 .. math::

 		f(x) \approx \sum_{i=0}^3 b_i B_i(x)

for some coefficients :math:`\left\{b_i \right\}_{i=0}^3`.  Using this approach the "infinite dimensional" function can now be represented by the four couefficients :math:`\left\{b_i \right\}_{i=0}^3`.

We still need to find the coefficients :math:`b\equiv \left\{b_i \right\}_{i=0}^3` that approximate our function. In our application we are looking for an unknown function that is implied by the Bellman equation (or Euler equation), but for now we will suppose that we know how to evaluate the function :math:`f(x)` at any value of :math:`x`.  We can then choose a few points :math:`X \subset \mathcal{X}` and evaluate :math:`f(x)` to define :math:`Y = \left\{f(x) : x \in X \right\}.`  Suppose there are :math:`J` points in :math:`X` and :math:`Y`.  We can then set up the following system of :math:`J` equations in our four unknown coefficients :math:`b`
 .. math::

 		\underbrace{
 		\left( \begin{array}{cccc}
			1 & x_1 & x_1^2 & x_1^3 \\
			1 & x_2 & x_2^2 & x_2^3 \\
			\vdots & \vdots & \vdots & \vdots \\
			1 & x_J & x_J^2 & x_J^3 \\
			\end{array} \right)
			}_{\equiv B(X)}
		\left( \begin{array}{c}
			b_0 \\
			b_1 \\
			b_2 \\
			b_3
			\end{array} \right)
			=
		\left( \begin{array}{c}
			y_1 \\
			y_2 \\
			\vdots \\
			y_J
			\end{array} \right).

In the matrix :math:`B(X)` is the *basis matrix* in which each row corresponds to one of our points in :math:`X` and each column corresponds to a basis function evaluated at that value of :math:`x`.  If we have four distinct points in :math:`X` then the basis matrix will be invertable and we can solve immediately for :math:`b`.  If  :math:`J>4` the system is overdetermined and there may not be any :math:`b` for which the equations hold exactly.  In this case we can find the :math:`b` that minimizes the sum of squared errors in the equations.  Some matrix calculus leads to the usual least-squares formula
 .. math::

 		b = \left[ B(X)^T B(X)\right]^{-1} B(X)^T Y

where the :math:`T` superscript indicates a matrix transpose.


Interpolation: using the approximation
----------------------------------------

Suppose we only know the value of the function at the points :math:`X` and don't know the value of the function at other points in :math:`\mathcal{X}`.  Suppose there is some set :math:`\tilde X` of points for which we would like to know the value of :math:`f(x)`.  We can use our approximation of the function to *interpolate* the values at :math:`X` to the values  at :math:`\tilde X`.  Doing so is simply a matter of constructing our basis functions at :math:`\tilde X` and weighting them with our coefficients :math:`b`.  We can phrase this in terms of matrix multiplication
 .. math::

 			\tilde Y = B(\tilde X) b

where the column vector :math:`\tilde Y` now contains the function values associated with :math:`\tilde X`.

Assessing accuracy
----------------------------------------

The easiest way to check the accuracy of our approximation is to select a set of points :math:`\tilde X \neq X` and then evaluate the function at those points both directly and through interpolation and check the difference
 .. math::

 			R(\tilde X) = B(\tilde X) b - f(\tilde X)

where :math:`R(\tilde X)` now contains the residuals at :math:`\tilde X`.  You might summarize this vector of residuals with the maximum absolute value, which is called the supremum norm and often written :math:`\| R(\tilde X) \|_\infty`.


Choosing the basis functions and  grid
----------------------------------------

In this class we will keep things simple.  We will work with the simple polynomial basis functions as shown above and we will choose evenly spaced grid points :math:`X \subset \mathcal{X}.` In a real applicaiton you might make different choices.  Using simple polynomials can lead to a basis matrix that is hard to badly behaved meaning that small deviations in :math:`Y` will lead to very different outcomes for :math:`b`.  This is a problem because even if you knew :math:`Y` exactly, a computer is forced to round it to some finite number of decimal places.  This is not likely to be a problem with a cubic polynomial but if you want a more accurate approximation you will need to include more basis functions and as the degree of the polynomial increases this "ill-conditioning" becomes a problem.  In that case you may want to learn about Chebyshev polynomials.

Similarly, the way we select the grid points will affect the accuracy of our approximation.  First, more grid points provide more information about the function.  However, the more gridpoints we have the longer it will take your computer to evaluate the function at all of the grid points and solve the least-squares fitting problem.  So there is a trade off between speed and accuracy.  Second, grid points near the boundaries of :math:`X` provide more information than those near the center (for a cubic polynomial we can learn a lot about :math:`b_3` by observing the function value at large positive and negative values of :math:`x` while observing :math:`f(0)` provides no information.)  While we will use evenly spaced grid points, you might want to learn about the Chebyshev nodes, which can increase accuracy for a given number of grid points.


Application to the stochastic growth model
--------------------------------------------

We will now apply these ideas to the stochastic growth model.  In that application, the functions we want to approximate are functions of two variables: the capital stock and the level of productivity.  We will use the components of a complete second order polynomial as our basis functions so our basis matrix will have the form
 .. math::

 		B(K,Z) =  \left( \begin{array}{cccccc}
			1 & K_1 & Z_1 & K_1^2 & K_1 Z_1 & Z_1^2 \\
			1 & K_2 & Z_2 & K_2^2 & K_2 Z_2 & Z_2^2 \\
			\vdots & \vdots & \vdots & \vdots & \vdots & \vdots \\
			1 & K_J & Z_J & K_J^2 & K_J Z_J & Z_J^2
			\end{array} \right)

and our basis coefficient vector :math:`b` will have size elements.

Turning to the computer, it will be convenient to have a function that takes the grid points :math:`K` and :math:`Z` and generates the basis matrix.  In Matlab we can do this as follows

.. highlight:: matlab

::

	function B  = PolyBasis(K,Z)
	% B  = PolyBasis(K,Z)
	% Polynomial basis functions.  Using 2nd order polynomial
	%
	% inputs
	% K    n x 1   points for K
	% Z    n x 1   points for Z
	%     or scalar for Z
	%
	% outputs
	% B    n x 6   array of basis functions: 1, K, Z, K^2, K*Z, Z^2

	Zb = Z.*ones(size(K));
	B = [ones(size(K)) K Zb K.^2 K.*Zb Zb.^2];

	end

Let me explain some of the syntax here.  The main work of the function is done by the line ``B = [ones(size(K)) K Zb K.^2 K.*Zb Zb.^2];`` which constructs the basis matrix and assigns it to the output variable ``B`` that was specified in the function header ``function B  = PolyBasis(K,Z)``.  In constructing the basis matrix, the square brackets indicated concatenation of arrays and the dot-operators ``.^`` and ``.*`` indicate element-by-element operations as opposed to ``^`` and ``*`` which mean matrix power and matrix multiplication.  You might be wondering why we take the array ``Z`` and multiply it by a matrix of ones ``Zb = Z.*ones(size(K));``  The reason is that in some cases we will want to pass an argument ``Z`` that is a scalar and an argument ``K`` that is an array and in this case we need to create an array ``Zb`` that is the same shape as ``K`` that can concatenated with ``K``.  If ``Z`` is already the same shape as ``K`` this step has no effect.


It will also be convenient to have a function that takes values for the capital stock and TFP and the associated function value (e.g. the value function at those points) and produces the basis coefficients for our polynomial approximation.

::

	function b = PolyGetCoef(K,Z,Y)
	% b = PolyGetCoef(Grid,Y)
	%   Fits the polynomial from PolyBasis to the function(s) in column(s) of
	%   Y.
	%
	% inputs
	% K    n x 1   points for K
	% Z    n x 1   points for Z
	% Y    n x 1   valies for function at (K,Z)
	%
	% outputs
	% b    6 x 1   basis coefficients

	b = PolyBasis(K,Z) \ Y;

	end

That was easy.  The operator ``\`` solves a system of linear equations and if it is over-determined it provides the least-squares solution.  You can read more about it in the `Matlab documentation <http://www.mathworks.com/help/matlab/ref/mldivide.html>`_.  An equivalent way of doing this would be

::

	B = PolyBasis(K,Z);
	b = inv(B' * B) * B' * Y;

where ``'`` indicates transpose.


Taking expectations on the computer using a Markov chain
==================================================================

In the Bellman equation and in the Euler equation we have the expectation operator.  This expectation is taken over the value of :math:`\varepsilon'`, which will imply :math:`Z'`.  Put differently, we are taking expectations of :math:`Z'` conditional on the current :math:`Z`.   An expectation is an integral and an integral can be thought of as a the limit of a Riemann sum as the number of rectangles (or other shapes) goes to infinity.  On the computer we will use the same idea as a Riemann sum, just a finite number of them.

So suppose we have a distribution over :math:`x\in\mathcal{X}` with pdf :math:`p` and CDF :math:`P`.  To take the expectation of :math:`f(x)` we will choose a set of *quadrature points* :math:`X \in \mathcal{X}` and we will construct some *quadrature weights* :math:`Q` that represent the probability associated with :math:`X`.  Suppose we have :math:`J` points in :math:`X` and :math:`Q`.  Our expectation of :math:`f(x)` is then given by
 .. math::

 		\mathbb E \left[ f(x) \right] = \int_{x\in\mathcal{X}} f(x) p(x) dx \approx \sum_{j=1}^J f(X_j) Q_j

There are many ways we could choose the quadrature points and weights.  We will use a method that is particularly conveneint for our purposes where we are interested in taking expectations of a stochastic process that follows an AR(1) process.  This is the `Tauchen (1986) <http://dx.doi.org/10.1016/0165-1765(86)90168-0>`_ method for discretizing an AR(1).

By repeated substitution the random variable :math:`Z_t` can be considered as an infinite sum
 .. math::

   Z_t = \sum_{\tau = -\infty}^t\rho^{t-\tau} \varepsilon_\tau.

As the weighted sum of normally distributed variables is it self normal we have that the unconditional distribution of :math:`Z_t` is normal with mean given by
 .. math::

   \sum_{\tau = -\infty}^t\rho^{t-\tau} \mathbb E \left[ \varepsilon_\tau \right] = 0.

and variance given by
 .. math::

   \sum_{\tau = -\infty}^t\rho^{2(t-\tau)} \sigma^2 = \frac{\sigma^2}{1-\rho^2}.

We will put an evenly-spaced grid on :math:`Z_t` on the interval that is centered at zero and spans two standard deviations of its unconditional distribution in each direction.  For each :math:`Z_i` in this grid, we will form the conditional expectation over :math:`Z'` by using a set of quadrature points for :math:`\varepsilon'` such that the resulting values for :math:`Z' = \rho Z_i + \varepsilon'` are exactly the points in our grid on :math:`Z`.  So for each :math:`Z_i` we have a probability of moving to the other discrete values in the grid :math:`Z_j`.  The result is a Markov chain that approximates the continuous AR(1) process.

To create the transition probabilities for the Markov chain, the Tauchen method simply assigns the probability mass in an interval between two grid points to the nearer of the two grid points. The following diagram shows how the Tauchen method creates the transition probability :math:`Q_{i,j}` of moving from :math:`Z_i` to :math:`Z_j`.

.. image:: figs/Tauchen.png
    :width: 676px
    :align: center
    :height: 394px
    :alt: Tauchen diagram

I have chosen to present the Tauchen method because it is easy to understand. There are other `similar methods <http://dx.doi.org/10.1016/j.red.2010.02.002>`_ that have been shown to be more accurate and reliable.

Once we have the Markov chain approximation to the AR(1), taking a conditional expectation of a function of :math:`Z'` conditional on :math:`Z` is very easy.  First, let's store the transition probabilities in a transition matrix for which the columns correspond to the current state :math:`Z` and for which the rows correspond to the future state :math:`Z'`
 .. math::

    Q \equiv \left( \begin{array}{c c c c}
			Q_{1,1} & Q_{2,1} & \cdots Q_{J,1} \\
      Q_{1,2} & Q_{2,2} & \cdots Q_{J,2} \\
			\vdots & \vdots & \ddots & \vdots  \\
			Q_{1,J} & Q_{2,J} & \cdots Q_{J,J}
			\end{array} \right).

Second, we compute the value of the function at each grid point for :math:`Z`.  Let's store that in a row vector :math:`Y`
 .. math::

    Y \equiv \left( \begin{array}{c c c c}
			f(Z_1) & f(Z_2) & \cdots & f(Z_J)
			\end{array} \right).

Then if we take the matrix product :math:`Y Q` we have a new row vector whose elements are :math:`\mathbb E \left[ f(Z') | Z_i \right]`.


Numerical maximization
========================

In the Bellman equation we maximize over :math:`(C,K')`.  By substituting the aggregate resource constraint we can phrase this as a maximization over one choice variable
 .. math::

   V(K,Z) &=  \max_{K'} \left\{ u\left( f(K,Z) - K' \right) + \beta \mathbb E V(K',Z') \right\}.

Optimization is a large field of numerical analysis, but in these notes we will consider just one method that is easy to understand and fairly robust called "golden section search".  Suppose we have a function :math:`f(x)` and we want to maximize it on the interval :math:`[\underline{x},\bar x]`. We start by evaluating the function at four points labeled :math:`A`, :math:`B`, :math:`C`, and :math:`D` in the following diagram.

.. image:: figs/Golden_1.png
    :width: 676px
    :align: center
    :height: 394px
    :alt: Golden section search diagram

To start with we would set :math:`A = \underline{x}` and :math:`D = \bar x`.  We then observe that :math:`f(C)>f(B)` so we use the following reasoning to discard the interval :math:`[A,B]` from containing the maximizer: if :math:`f(x)` is concave then the derivative is weakly decreasing and if the true maximizer :math:`x^*` were in :math:`[A,B]` then the function goes down from :math:`x^*` to B and then up from B to C, which contradicts the supposition that :math:`f` is concave.  You could use this algorithm if :math:`f` is not concave but then there are no guarantees that it will work.

Once we eliminate the interval :math:`[A,B]` we repeat the same steps on the interval :math:`[B,D]` that is smaller than our original interval :math:`[A,D]`.  So point :math:`B` will become the new :math:`A`.  The clever part of the algorithm is that we put the points :math:`B` and :math:`C` in just the right places so that now :math:`C` becomes the new :math:`B`.  This means we don't need to re-evaluate the function at our new points :math:`A` or :math:`B`.  We now introduce a new point :math:`C`, evaluate the function there and repeat as shown in the next diagram.

.. image:: figs/Golden_2.png
    :width: 676px
    :align: center
    :height: 394px
    :alt: Golden section search diagram

In the next iteration we find :math:`f(B)>f(C)` so now we discard :math:`[C,D]`, :math:`C` becomes the new :math:`D`, and :math:`B` becomes the new :math:`C`.

.. image:: figs/Golden_3.png
    :width: 676px
    :align: center
    :height: 394px
    :alt: Golden section search diagram

As you can see in the diagrams, the region under consideration is shrinking towards the true maximizer.  If we repeat this procedure many times over we will have :math:`[A,D]` tightly bracketing :math:`x^*`.

Finally, how do we choose the points :math:`B` and :math:`C`?  We are going to impose some symmetry on the choice so the interval :math:`[A,B]` is the same length as :math:`[C,D]`. This reduces the problem to choosing one number :math:`p \in (0,1)` as shown in the following figure.

.. image:: figs/GoldenSections.png
    :width: 516px
    :align: center
    :height: 198px
    :alt: Golden sections diagram

We want to choose :math:`p` so that when we eliminate :math:`[A,B]` point :math:`C` becomes the new point :math:`B` so we need the ratio of :math:`BC` relative to :math:`BD` to be the same as :math:`AB` relative to :math:`AD`.  To accomplish this we need
 .. math::

    \frac{2p-1}{p} = 1-p.

This equation reduces to a quadratic equation with one root :math:`p = (-1 \+ \sqrt{5})/2`.  As the other root implies :math:`p<0` we disregard it.

Why is the method called "golden section search"? Because :math:`1/p` is the `golden ratio <http://mathworld.wolfram.com/GoldenRatio.html>`_.

Value function iteration
==================================================================

Setting up the problem
-----------------------

We are now ready to solve the model using our first method, which is value function iteration.  To do so, we are going to have to set up the model in the computer, by which I mean we are going to tell the computer what we know about the problem.  We can begin as follows:

::

  % Store the parameters in a structure
  Par.beta = 0.99;
  Par.gamma = 2;
  Par.alpha = 0.36;
  Par.delta = 0.03;
  Par.rho = 0.95;
  Par.sigma = 0.007;


  %% Solve for the steady state
  Zstar = 0;
  Kstar = ((1/Par.beta - 1 + Par.delta)./(Par.alpha * exp(Zstar))).^(1/(Par.alpha-1));


We store the parameters in a Matlab structure, which is a grouping of variables.  This might seem like extra work but it will be useful because then we can pass the whole structure to a function instead of having to specify all of the individual variables.  We then solve for the steady state capital stock from the equation
 .. math::

      f'(K^*) = 1/\beta.

The next step is to set up our grid on :math:`K` and :math:`Z`.  This grid needs to cover a two-dimensional space so we are going to create two one-dimensional grids and then take the Cartesian product of the two.
::

  %% Create a grid for Z
  meanZ = 0;
  Grid.nZ = 7;  % number of points in our grid for Z
  numStdZ = 2;  % number of standard deviations to cover with the grid
  [Grid.Z, Grid.Q]  = tauchen(Grid.nZ, meanZ, Par.rho, Par.sigma, numStdZ);

  Grid.Q = Grid.Q'; % this is a 7 x 7 transition matrix for which the columns sum to 1
  % the (i,j) element is the probability of moving from j to i.

To create the grid for :math:`Z` and the transition matrix we use an implementation of the Tauchen method.  We discussed this above and we won't go further into the details here.  We create a ``Grid`` structure to hold all of our variables associated with the grids.
::

  %% Create a grid for K
  Grid.nK = 20;
  Grid.K = linspace(0.75*Kstar, 1.25*Kstar,Grid.nK)';  % this is a 20 x 1 array of evenly spaced points

The grid for :math:`K` is made up of a set of 20 evenly spaced points that extends from 25 percent below the steady state to 25 percent above.
::

  %% Create a product of the two grids
  [ZZ,KK] =meshgrid(Grid.Z,Grid.K);
  Grid.KK = KK(:);
  Grid.ZZ = ZZ(:);

The Matlab function ``meshgrid`` creates the Cartesian product of two vectors.  The output of ``meshgrid`` is two 20 x 7 arrays with  ``KK(i,j)`` equal to ``K(i)`` and  ``ZZ(i,j)`` equal to ``Z(j)``.  The next line unravels the array ``KK`` into a single column and stores it in the ``Grid`` structure.  The use of the `colon operator <http://www.mathworks.com/help/matlab/ref/colon.html>`_ ``(:)`` means treat the elements of this array as a single column regarless of the shape of the array.  The last line unravels ``ZZ`` in the same way.

How does Matlab store and reference arrays?
-------------------------------------------

A small detour about how Matlab stores arrays.  Suppose you create the array
::

  A = [1 2;
       3 4;
       5 6]

This is an array with three rows and two columns.  If we type ``A(:)`` we find
::

  >> A(:)

  ans =

     1
     3
     5
     2
     4
     6

Notice that Matlab reads down the columns first and then down the rows.    This works the other way around, too.  If we create a vector and reshape it into an array, Matlab will fill the array up the column by column
::

  >> A = 1:4

  A =

       1     2     3     4

  >> reshape(A,2,2)

  ans =

       1     3
       2     4


The Bellman operator
---------------------

We now code the Bellman equation into the computer and we do this in two steps. First, we write a function that evaluates the right hand side of the Bellman equation for a given choice of :math:`K'`. Second, we write a function that will maximize the right hand side of the Bellman equation. Before we begin, it will be convenient to have a function that evaluates our production function :math:`f(K,Z)`.
::

  function y  = f(Par, K,Z )
  % y  = f( K,Z )
  %   Production function gross of undepreciated capital

  y =  exp(Z) .* K.^Par.alpha + (1-Par.delta)*K;
  end

Next we code the Bellman equation given a level of savings.
::

  function V  = Bellman( Par, b, K, Z, Kp )
  % V  = Bellman( Par, b, K, Z, Kp )
  %   Evaluate the RHS of the Bellman equation
  %
  % Inputs
  % Par   Parameter structure
  % b     6 x 1 coefficients in polynomial for E[ V(K',Z') | Z ]
  % K     n x 1 array of current capital
  % Z     n x 1 array of current TFP
  % Kp    n x 1 array of savings
  %
  % Output
  % V     n x 1 array of value function
  %

  C = f(Par,K,Z) - Kp;
  u = C.^(1-Par.gamma) / (1-Par.gamma);
  V = u + Par.beta * PolyBasis(Kp,Z) * b;



  end

The first line calculates consumption from the resource constraint and the next line computes the period utility from that consumption. The third line creates the output variable ``V`` equal to the sum of the period and continuaion utilities.

Notice that we are not actually going to approximate the value function but rather we are approximating :math:`\mathbb E \left[ V \left( K' ,Z' \right) | Z \right]`.  Therefore the function we are approximating is a function of how much we save :math:`K'`, which is known today, and the current level of TFP :math:`Z`.  To evaluate that function we create our polynomial basis matrix and then multiply it against the coefficients of our polynomial.

Notice that the ``Bellman`` function takes the level of savings as an input, but the Bellman equation involves maximizing over this choice variable.


We will perform that maximization

The idea of value function iteration is that we can start with any value function and apply the Bellman operator repeatedly to iterate towards the true value function.  So we need a guess of the value function to start with.  We will do something naive and guess that it is the zero function meaning expected value is also the zero function and the coefficients of our approximating polynomial are all zero.
::

  b = zeros(6,1);


We can now do the maximization in the Bellman equation using golden section search and we will create a function to do this.
::

  function [V, Kp] = MaxBellman(Par,b,Grid)
  % [V, Kp] = MaxBellman(Par,b,Grid)
  %   Maximizes the RHS of the Bellman equation using golden section search
  %
  % Inputs
  % Par       Parameter structure
  % b     6 x 1 coefficients in polynomial for E[ V(K',Z') | Z ]
  % Grid      Grid structure


  p = (sqrt(5)-1)/2;

  A = Grid.K(1) * ones(size(Grid.KK));
  D = min(f(Par,Grid.KK,Grid.ZZ) - 1e-3, Grid.K(end)); % -1e-3 so we always have positve consumption.

  B = p*A+(1-p)*D;
  C = (1-p)*A + p * D;


  fB = Bellman(Par,b,Grid.KK,Grid.ZZ,B);
  fC = Bellman(Par,b,Grid.KK,Grid.ZZ,C);


  MAXIT = 1000;
  for it_inner = 1:MAXIT

      if all(D-A < 1e-6)
          break
      end

      I = fB > fC;

      D(I) = C(I);
      C(I) = B(I);
      fC(I) = fB(I);
      B(I) = p*C(I) + (1-p)*A(I);
      fB(I) = Bellman(Par,b,Grid.KK(I),Grid.ZZ(I),B(I));

      A(~I) = B(~I);
      B(~I) = C(~I);
      fB(~I) = fC(~I);
      C(~I) = p*B(~I) + (1-p)*D(~I);
      fC(~I) = Bellman(Par,b,Grid.KK(~I),Grid.ZZ(~I),C(~I));

  end

  % At this stage, A, B, C, and D are all within a small epsilon of one
  % another.  We will use the average of B and C as the optimal level of
  % savings.
  Kp = (B+C)/2;

  % evaluate the Bellman equation at the optimal policy to find the new
  % value function.
  V = Bellman(Par,b,Grid.KK,Grid.ZZ,Kp);

  end

We start by defining the constant :math:`p` for the golden section search.  We then specify the interval of :math:`K'` that we will search over and we will assume that we never leave the grid on the low side and we never save so much as to have negative consumption on the high side and we will never save so much as to leave the grid.  Notice that these are vectors that corresond to the size of our grid vectors ``KK`` and ``ZZ``.  We then define the points ``B`` and ``C`` using the golden section ratios and evaluate the ``Bellman`` function at those points.  Notice that everything we are doing here is operating on vectors of :math:`K'` that correspond to the choices at the corresponding levels of :math:`K` and :math:`Z` that appear in ``Grid.KK`` and ``Grid.ZZ``.

 Next we enter a loop for 1,000 iterations.  We will continue shrinking the interval until the distance between ``D`` and ``A`` is sufficiently small as shown in the first lines of the loop code.  The comparison ``D - A < 1e-6`` will generate a logical array (an array of true and false values) that checks the distance between each element of ``D`` and each corresponding element of ``A``.  The ``all`` function will then tell us whether all of the values in the logical array are true.   The ``break`` command tells the program to leave the loop it is executing.  One could use a ``while`` loop for this, but a danger with that approach is that if there is a bug the program could run forever in a ``while`` loop whereas with the ``for`` loop it will stop after ``MAXIT`` iterations at most.


The next part of the algorithm is a little tricky.  When we discussed golden section search we were maximizing a scalar function with respect to a single argument.  But in this application we are maximizing different functions (if we think of the variation in the rows of ``Grid.KK`` and ``Grid.ZZ`` giving rise to different ``Bellman`` functions for ``Kp``) with respect to a vector of values ``Kp``.  Each row can be considered a scalar function of a scalar argument but we are doing many independent maximizations at the same time.  In this case we won't necessarily find that ``fB > fC`` for all of the rows, nor vice versa. So we define an indicator function or "logical arrays" for whether ``fB > fC``.  In the lines that follow we will first do the operations for those rows for which the logical array is true and then we will do the operations for which the rows are false using ``~I`` because the ``~`` inverts the logical values in the array.  In either case the operations are straightforward.  For the cases for which ``fB > fC`` we are going to eliminate the interval :math:`CD` so for these rows :math:`C` becomes the new :math:`D` which we accomplish with ``D(I) = C(I)`` then :math:`B` becomes the new :math:`C` which we accomplish with ``C(I) = B(I)`` and the function values at :math:`B` become the function values at :math:`C` which we accomplish with ``fC(I) = fB(I)``.  Next we create a new :math:`B` and evaluate the function at those points. The remaining block of code performs the analogous operations for the cases ``fB < fC`` where we eliminate the interval :math:`AB`.

After we have exited the loop, we just have a little more work to do.  By construction, ``A``, ``B``, ``C``, and ``D`` are all within ``1e-6`` of one another so we will take the average of ``B`` and ``C`` and treat it as the true maximizer ``Kp``.  Next we evaluate the ``Bellman`` function at ``Kp`` to get the maximized value.

Updating our guess of the value function
-----------------------------------------

After applying ``MaxBellman`` we have a vector of values that give us :math:`V(K,Z)` on our grid for :math:`K` and :math:`Z`.  But this value function was calculated for a given guess of :math:`\mathbb E \left[ V \left( K' ,Z' \right) | Z \right]`, which was not necessarily the true expected continuation value.  Now we are going to use our new :math:`V(K,Z)` to update our guess of :math:`\mathbb E \left[ V \left( K' ,Z' \right) | Z \right]`.  To do so is easy given that we have already done the work of discretizing our AR(1) process for :math:`Z`.

The first step is to reshape ``V`` into a 20 by 7 array where the rows correspond to different levels of capital and the columns correspond to different levels of TFP.  Because we were a little bit clever in how we set up the arrays ``Grid.KK`` and ``Grid.ZZ``, this is just a matter of calling Matlab's ``reshape`` command like this ``reshape(V,Grid.nK,Grid.nZ)``.  So a given row of this array is the value function at the same level of capital but at different levels of TFP.  To compute the conditional expectation conditional on the previous level of TFP, we just need to take the dot product of this row of the array with the appropriate column of the Markov chain transition matrix.  We can get all of the conditional expectations with a matrix multiplication like this
::

  EV = reshape(V,Grid.nK,Grid.nZ) * Grid.PZ;

``EV`` is an array where each row corresponds to :math:`K_t`, each column corresponds to :math:`Z_{t-1}` and the entires are :math:`\mathbb E \left[ V( K_t, Z_t) | Z_{t-1} \right]`.

Now that we have computed the expected value on our grid, we can update the coefficients of the polynomial that approximates this function using the ``PolyGetCoef`` function that we wrote earlier.
::

  b1 = PolyGetCoef(Grid.KK,Grid.ZZ,EV(:));

We use ``(:)`` to unravel the ``EV`` array into one column that matches the layout of ``Grid.KK`` and ``Grid.ZZ``.


Putting the *iteration* in value function iteration
--------------------------------------------------------------

We have now completed one step in our algorithm.  We started with a guess of the value function polynomial coefficients and we arrived at a new set of coefficients.  We are now going to repeat the process many times over using a ``for`` loop until the value function has converged.

How do we know when to stop? A good approach is to check for convergence in terms of something that you actually are interested in.  So instead of checking that the polynomial coefficients have converged, let's check that the policy rule has converged.  To do that, let's introduce a copy of the previous policy rule, call it ``Kp0`` and compare the new policy rule to it.
::

  Kp0 = zeros(size(Grid.KK));
  MAXIT = 2000;
  for it = 1:MAXIT

      [V, Kp] = MaxBellman(Par,b,Grid);


      % take the expectation of the value function from the perspective of
      % the previous Z
      EV = reshape(V,Grid.nK,Grid.nZ) * Grid.PZ;

      % update our polynomial coefficients
      b = PolyGetCoef(Grid.KK,Grid.ZZ,EV(:));

      % see how much our coefficients have changed
      test = max(abs(Kp0 - Kp));
      Kp0 = Kp;

      disp(['iteration ' num2str(it) ', test = ' num2str(test)])
      if test < 1e-5
          break
      end
  end


We're done!  Kind of.
----------------------------

We now have all of the pieces we need solve for the value function and policy rule using value function iteration.  If you run the script ``VFI.m`` you should see your computer iterate for about 200 iterations.

Now that we have solved the problem, we can look at the results.  The following commands make a plot of the policy rules
::

  DK = Grid.K/Kstar-1; % Capital grid as percent deviation from steady state

  DKp = reshape(Kp,Grid.nK,Grid.nZ)./reshape(Grid.KK,Grid.nK,Grid.nZ) - 1;
    % savings policy rule as a 20 x 7 array expressed as a percent change from current K

  plot(DK, DKp);  % plot the policy rule

  hold on;        % next plots go on the same figure

  plot(DK, zeros(Grid.nK,1), 'k--'); % add a zero line, k-- means black and dashsed

  xlabel('K in % deviation from steady state')  % label the axes
  ylabel('(K'' - K)/K')


.. image:: figs/VFI_policy_rules.png
      :width: 563px
      :align: center
      :height: 422
      :alt: Policy rules


The horizontal axis is the current capital stock as a percentage difference from the steady state and the vertical axis is the next period's capital stock as a percentage difference from this period's.  There are seven lines on the figure corresponding to the different levels of current TFP.  The lowest line corresponds to the lowest TFP and so on.  As you can see, if current capital is low enough then the optimal policy is to increase the capital stock for all TFP levels and if the current capital is high enough then the optimal policy is to eat down the capital stock for all TFP levels.  In between, the optimal policy is to build capital if TFP is high and reduce capital if TFP is low.

Suppose now we want to know the optimal policy at a point not on our grid.  We can use interpolation to infer this from the value that we computed on the grid.
::

  bKp =  PolyGetCoef(Grid.KK,Grid.ZZ,Kp);
  Kp2903 = PolyBasis(29,0.03) * bKp
  C2903 = f(Par,29,0.03) - Kp2903

First we fit a polynomial to the policy rule and then we evaluate that polynomial at the point :math:`K = 29` :math:`Z = 0.03` to get the interpolated savings at that point.  Finally we can use the aggregate resource constraint to find consumption at that point.
