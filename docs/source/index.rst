.. NumericalAnalysis documentation master file, created by
   sphinx-quickstart on Thu Aug 11 20:18:56 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Numerical analysis notes for EC 702
=============================================

.. toctree::
   :maxdepth: 1

   Function approximation <FuncApprox>
   Taking expecations on the computer <Expectations>
   Numerical maximization <Maximize>
   Value function iteration <VFI>
   Iterating on the Euler equation <EulerIt>
   Simulating the solution <Simulation>
   Assessing the accuracy of our solutions <Accuracy>

The programs associated with these notes can be accessed on `GitHub <https://github.com/amckay/EC702>`_.

Overview
=================================================

We started this part of the course asking what does it mean to solve a model?  Our answer was that solving a model means finding a function or stochastic process that satisfies some conditions.  But how do you do that?  As we have seen, in very special cases we can do this analytically but in most cases there is no closed-form solution.  Most macroeconomic models are analyzed numerically.  These notes present two approaches to computing numerical solutions of models. The first approach searches for a value function that satisfies the Bellman equation.  The second approach searches for a policy rules that satisfy the equilibrium conditions of the model (e.g. we solve the functional Euler equation).  The advantages of value function methods are that they are applicable to complex models and they have good convergence properties.  The advantage of working directly with the equilibrium conditions is that they are often faster. In the course of these notes we will cover some fundamental techniques of numerical analysis.  We will discuss how we can represent  function on the computer and we will discuss how we can take expectations on the computer.


In these notes, I will use the programming language Matlab to demonstrate the methods.  For the applications we will consider we could use any number of different languages, but graduate students typically start out with Matlab and so that is why I have chosen it.  I think Matlab's popularity stems in part from the fact that it has good documentation.  There are some disadvantages of Matlab that might lead you to consider alternatives.  First, you need to buy a license in order to install it on your computer while some of the alternatives are free open-source software.  Second, for some applications Matlab is slower than some of the alternatives.  You can find a thorough comparison of programming languages from the perspective of economists `here <http://www.nber.org/papers/w20263>`_.  MathWorks, the company that produces Matlab, has `resources to help you get started with Matlab if it is new to you. <http://www.mathworks.com/help/matlab/examples.html>`_



Application: the stochastic neoclassical growth model
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
        K' &= f(K,Z) - C(K,Z) \\
        f(K,Z) &= e^Z K^\alpha + (1-\delta)K \\
        Z' &= \rho Z + \varepsilon'.

Our second approach to solving the model will be to find a policy rule for consumption that satisfies the Euler equation.

In order to solve the model numerically, we need values for the parameters.  We will discuss how we use data to choose parameters, but for now we will simply assume :math:`\gamma = 2`, :math:`\beta = 0.99`, :math:`\beta = 0.99`, :math:`\delta = 0.03`, and :math:`\rho = 0.95`.  We also need a distribution for the innovations and we will assume :math:`\varepsilon \sim N(0,\sigma^2)` where the standard deviation :math:`\sigma = 0.007`.


Before we get to the discussion of solving the model, we will need to discuss some fundamental techniques of numerical analysis.
