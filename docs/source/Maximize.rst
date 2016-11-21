.. NumericalAnalysis documentation master file, created by
   sphinx-quickstart on Thu Aug 11 20:18:56 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Numerical maximization
========================

In the Bellman equation we maximize over :math:`(C,K')`.  By substituting the aggregate resource constraint we can phrase this as a maximization over one choice variable
 .. math::

   V(K,Z) &=  \max_{K'} \left\{ u\left( f(K,Z) - K' \right) + \beta \mathbb E V(K',Z') \right\}.

Optimization is a large field of numerical analysis, but in these notes we will consider just one method that is easy to understand and fairly robust called "golden section search".  Suppose we have a function :math:`F(x)` and we want to maximize it on the interval :math:`[\underline{x},\bar x]`. We start by evaluating the function at four points labeled :math:`A`, :math:`B`, :math:`C`, and :math:`D` in the following diagram.

.. image:: figs/Golden_1.png
    :width: 676px
    :align: center
    :height: 394px
    :alt: Golden section search diagram

To start with we would set :math:`A = \underline{x}` and :math:`D = \bar x`.  We then observe that :math:`F(C)>F(B)` so we use the following reasoning to discard the interval :math:`[A,B]` from containing the maximizer: if :math:`F(x)` is concave then the derivative is weakly decreasing and if the true maximizer :math:`x^*` were in :math:`[A,B]` then the function goes down from :math:`x^*` to B and then up from B to C, which contradicts the supposition that :math:`f` is concave.  You could use this algorithm if :math:`f` is not concave but then there are no guarantees that it will work.

Once we eliminate the interval :math:`[A,B]` we repeat the same steps on the interval :math:`[B,D]` that is smaller than our original interval :math:`[A,D]`.  So point :math:`B` will become the new :math:`A`.  The clever part of the algorithm is that we put the points :math:`B` and :math:`C` in just the right places so that now :math:`C` becomes the new :math:`B`.  This means we don't need to re-evaluate the function at our new points :math:`A` or :math:`B`.  We now introduce a new point :math:`C`, evaluate the function there and repeat as shown in the next diagram.

.. image:: figs/Golden_2.png
    :width: 676px
    :align: center
    :height: 394px
    :alt: Golden section search diagram

In the next iteration we find :math:`F(B)>F(C)` so now we discard :math:`[C,D]`, :math:`C` becomes the new :math:`D`, and :math:`B` becomes the new :math:`C`.

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

This equation reduces to a quadratic equation with one root :math:`p = (-1 + \sqrt{5})/2`.  As the other root implies :math:`p<0` we disregard it.

Why is the method called "golden section search"? Because :math:`1/p` is the `golden ratio <http://mathworld.wolfram.com/GoldenRatio.html>`_.
