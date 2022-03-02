======================
Holonomic  Constraints
======================

Four-Bar Linkage
================

Consider the linkage shown below:

.. _configuration-four-bar
.. figure:: figures/configuration-four-bar.svg
   :align: center

This is a planar `four bar linkage`_ with reference frames :math:`N,A,B,C`
attached to each bar. Four bar linkages are used in a wide variety of
mechanisms. One you may be familiar with is this rear supsension on a mountain
bicycle:

.. figure:: https://upload.wikimedia.org/wikipedia/commons/thumb/7/7c/MtbFrameGeometry_FSR.png/320px-MtbFrameGeometry_FSR.png
   :align: center

   Four bar linkage shown in blue, red, orange, and green used in the rear
   suspension mechanism of a mountain bicycle.

   Cartemere, CC BY-SA 3.0 <https://creativecommons.org/licenses/by-sa/3.0>, via Wikimedia Commons

A four bar linkage is an example of a *closed kinematic loop*. The the case of
:numref:`configuration-four-bar` there are two vector paths to point
:math:`P_4`:

.. math::

   \bar{r}^{P_4/P_1} = \bar{r}^{P_2/P_1} + \bar{r}^{P_3/P_2} + \bar{r}^{P_4/P_3}

For the loop to close, the two vector paths must equate. We can resolve this by
disconnecting the loop at some location, :math:`P_4` in our case, and forming
the open loop vector equations to points that should coincide.

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   me.init_vprinting()

Setup the variables, refrence frames, and points:

.. jupyter-execute::

   q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')
   la, lb, lc, ln = sm.symbols('l_a, l_b, l_c, l_n')

   N = me.ReferenceFrame('N')
   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')
   C = me.ReferenceFrame('C')

   A.orient_axis(N, q1, N.z)
   B.orient_axis(A, q2, A.z)
   C.orient_axis(B, q3, B.z)

   P1 = me.Point('P1')
   P2 = me.Point('P2')
   P3 = me.Point('P3')
   P4 = me.Point('P4')

SymPy Mechanics will warn you if you try to establish a closed loop among a set
of points and you should not do that. Instead you will establish positions
among points on one open leg of the chain:

.. jupyter-execute::

   P2.set_pos(P1, la*A.x)
   P3.set_pos(P2, lb*B.x)
   P4.set_pos(P3, lc*C.x)

   P4.pos_from(P1)

Now, declare a vector for the other path to :math:`P_4`:

.. jupyter-execute::

   r_P1_P4 = ln*N.x

Now we can form the left hand side of the following equation:

.. math::

   \bar{r}^{P_4/P_1} - \left( \bar{r}^{P_2/P_1} + \bar{r}^{P_3/P_2} + \bar{r}^{P_4/P_3} \right) = 0

Using :external:py:meth:`~sympy.physics.vector.point.Point.pos_from` for the
open loop leg made of points and the additional vector:

.. jupyter-execute::

   loop = P4.pos_from(P1) - r_P1_P4
   loop

This "loop" vector equation must equate to zero for our linkage to always be a
closed loop. We have a planar mechanism, so we can extract two scalar equations
associated with a pair of unit vectors in the plane of the mechanism:

.. jupyter-execute::

   fhx = sm.trigsimp(loop.dot(N.x))
   fhx

.. jupyter-execute::

   fhy = sm.trigsimp(loop.dot(N.y))
   fhy

For the loop to close, these two expressions must equal zero for all values
:math:`q_1,q_2,q_3`. These are two nonlinear equations in three time varying
varialbes. A solution, sometimes analytically but likely only numerical, can be
found if we solve for two of the time varying variables. For example,
:math:`q_2` and :math:`q_3` can be solved for in terms of :math:`q_1`. We would
then say that :math:`q_2` and :math:`q_3` depend on :math:`q_1`. These two
equations are called holonomic constraints, or configuration constraints
because they constrain the kinematic configuration to be a loop. Holonomic
constraints take the form:

.. math::
   :label: configuration-constraint

   \bar{f}_h(q_1, \ldots, q_n, t) = 0 \textrm{ where } \bar{f} \in \mathbb{R}^N

These constraints are functions of configuration variables: time varying angles
and distances.

.. jupyter-execute::

   q4 = sm.trigsimp(C.x.angle_between(N.x))
   q4

General Holonomic Constraints
=============================

If you consider each of the points :math:`P_1,P_2,P_3,P_4` 

points are located in Euclidean space by three scalars, one scalar
for each Cartesian coordinate.

