Configuration Constraints
=========================

Consider the linkage shown below

.. figure:: https://upload.wikimedia.org/wikipedia/commons/thumb/7/7c/MtbFrameGeometry_FSR.png/320px-MtbFrameGeometry_FSR.png

   Cartemere, CC BY-SA 3.0 <https://creativecommons.org/licenses/by-sa/3.0>, via Wikimedia Commons

.. figure:: figures/configuration-four-bar.svg

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   me.init_vprinting()


.. jupyter-execute::

   q1, q2, q3 = me.dynamicsymbols('q1, q2, q3', real=True)
   la, lb, lc, ln = sm.symbols('l_a, l_b, l_c, l_n', real=True, positive=True)

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

   P2.set_pos(P1, la*A.x)
   P3.set_pos(P2, lb*B.x)
   P4.set_pos(P3, lc*C.x)

   r_P1_P4 = ln*N.x

   loop = P4.pos_from(P1) - r_P1_P4
   loop

.. jupyter-execute::

   fx = sm.trigsimp(loop.dot(N.x))
   fx

.. jupyter-execute::

   fy = sm.trigsimp(loop.dot(N.y))
   fy

These two equations are called configuration constraints. Configuration
constraints take the form:

.. math::
   :label: configuration-constraint

   \bar{f}(q_1, \ldots, q_n) = 0 \textrm{where} \bar{f} \in \mathbb{R}^N

In general, points are located in Euclidean space by three scalars, one scalar
for each Cartesian coordinate.

.. jupyter-execute::

   q4 = sm.trigsimp(C.x.angle_between(N.x))
   q4
