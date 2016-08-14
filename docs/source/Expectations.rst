.. NumericalAnalysis documentation master file, created by
   sphinx-quickstart on Thu Aug 11 20:18:56 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


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
