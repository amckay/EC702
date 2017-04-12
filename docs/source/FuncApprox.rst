.. NumericalAnalysis documentation master file, created by
   sphinx-quickstart on Thu Aug 11 20:18:56 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.



Function approximation
========================================================

Suppose there is a function :math:`F(x)` defined on some domain :math:`\mathcal{X} \subset \mathbb R`.   People sometimes say that a function is "an infinite dimensional object" meaning that one could represent the function as an infinite table of values associated with each :math:`x \in \mathcal{X}`.  Representing a function this way would require an infinite amount of memory on a computer so we use a different approach.  Instead we will represnt :math:`F(x)` as a weighted sum of a finite number of pre-defined functions, which we call "basis functions."  For example we know that any continuous function can be well approximated by a polynomial so we could use polynomials as our basis functions and then look for the coefficients on the polynomials to approximate :math:`F(x)`.  Concretely, let :math:`B_i(x) = x^i` so we have :math:`B_0(x) = 1`, :math:`B_1(x) = x`, :math:`B_2(x) = x^2`, and so on.  Let's now suppose we will use a cubic polynomial to approximate :math:`F(x)`.  We then have
 .. math::

 		F(x) \approx \sum_{i=0}^3 b_i B_i(x)

for some coefficients :math:`\left\{b_i \right\}_{i=0}^3`.  Using this approach the "infinite dimensional" function can now be represented by the four couefficients :math:`\left\{b_i \right\}_{i=0}^3`.

We still need to find the coefficients :math:`b\equiv \left\{b_i \right\}_{i=0}^3` that approximate our function. In our application we are looking for an unknown function that is implied by the Bellman equation (or Euler equation), but for now we will suppose that we know how to evaluate the function :math:`F(x)` at any value of :math:`x`.  We can then choose a few points :math:`X \subset \mathcal{X}` and evaluate :math:`F(x)` to define :math:`Y = \left\{F(x) : x \in X \right\}.`  Suppose there are :math:`J` points in :math:`X` and :math:`Y`.  We can then set up the following system of :math:`J` equations in our four unknown coefficients :math:`b`
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

The matrix :math:`B(X)` is the *basis matrix* in which each row corresponds to one of our points in :math:`X` and each column corresponds to a basis function evaluated at that value of :math:`x`.  If we have four distinct points in :math:`X` then the basis matrix will be invertable and we can solve immediately for :math:`b`.  If  :math:`J>4` the system is overdetermined and there may not be any :math:`b` for which the equations hold exactly.  In this case we can find the :math:`b` that minimizes the sum of squared errors in the equations.  Some matrix calculus leads to the usual least-squares formula
 .. math::

 		b = \left[ B(X)^T B(X)\right]^{-1} B(X)^T Y

where the :math:`T` superscript indicates a matrix transpose.


Interpolation: using the approximation
----------------------------------------

Suppose we only know the value of the function at the points :math:`X` and don't know the value of the function at other points in :math:`\mathcal{X}`.  Suppose there is some set :math:`\tilde X` of points for which we would like to know the value of :math:`F(x)`.  We can use our approximation of the function to *interpolate* the values at :math:`X` to the values  at :math:`\tilde X`.  Doing so is simply a matter of constructing our basis functions at :math:`\tilde X` and weighting them with our coefficients :math:`b`.  We can phrase this in terms of matrix multiplication
 .. math::

 			\tilde Y = B(\tilde X) b

where the column vector :math:`\tilde Y` now contains the function values associated with :math:`\tilde X`.

Assessing accuracy
----------------------------------------

The easiest way to check the accuracy of our approximation is to select a set of points :math:`\tilde X \neq X` and then evaluate the function at those points both directly and through interpolation and check the difference
 .. math::

 			R(\tilde X) = B(\tilde X) b - F(\tilde X)

where :math:`R(\tilde X)` now contains the residuals at :math:`\tilde X`.  You might summarize this vector of residuals with the maximum absolute value, which is called the supremum norm and often written :math:`\| R(\tilde X) \|_\infty`.


Choosing the basis functions and  grid
----------------------------------------

In this class we will keep things simple.  We will work with the simple polynomial basis functions as shown above and we will choose evenly spaced grid points :math:`X \subset \mathcal{X}.` In a real applicaiton you might make different choices.  Using simple polynomials can lead to a basis matrix that is badly behaved in that small deviations in :math:`Y` will lead to very different outcomes for :math:`b`.  This is a problem because even if you knew :math:`Y` exactly, a computer is forced to round it to some finite number of decimal places.  This is not likely to be a problem with a cubic polynomial but if you want a more accurate approximation you will need to include more basis functions and as the degree of the polynomial increases this "ill-conditioning" becomes a problem.  In that case you may want to learn about Chebyshev polynomials.

Similarly, the way we select the grid points will affect the accuracy of our approximation.  First, more grid points provide more information about the function.  However, the more gridpoints we have the longer it will take your computer to evaluate the function at all of the grid points and solve the least-squares fitting problem.  So there is a trade off between speed and accuracy.  Second, grid points near the boundaries of :math:`X` provide more information than those near the center (for a cubic polynomial we can learn a lot about :math:`b_3` by observing the function value at large positive and negative values of :math:`x` while observing :math:`F(0)` provides no information.)  While we will use evenly spaced grid points, you might want to learn about the Chebyshev nodes, which can increase accuracy for a given number of grid points.


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

and our basis coefficient vector :math:`b` will have six elements.

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
