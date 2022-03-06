========================
Nonholonomic Constraints
========================

In the prior chapter, we discussed constraints on the configuration of a
system. The configuration considers where points are and how reference frames
are oriented. In this chapter, we will consider constraints on the motion of a
system.

Take for example parallel parking a car as a motivating example.

.. _motion-parallel:
.. figure:: figures/motion-parallel.svg
   :align: center

We know that car 2 can be in either the left or right location in a), i.e. the
car's configuration permits either location. But the scenario in b) isn't
possible. A car can't move from the left configuration to the right
configuration by simply moving to the right [*]_. Although, this surely would
be nice if we could. A car has wheels and only the front wheels can be steered,
so the scenario in c) is the only way for the car to end up in the right
configuration. The car has to move in a specific way to get from one
configruation to another. This entails that we have some kind of constraint on
the motion but not the configuration. Constraints such as these are called
*nonholonomic constraints* and they take the form:

.. math::
   :label: nonholonomic-constraints

   \bar{f}_n(\bar{u}, \bar{q}, t) = 0 \\
   \textrm{ where } \\
   \bar{f}_n \in \mathbb{R}^m \\
   \bar{u} = \left[ u_1, \ldots, u_n\right]^T \in \mathbb{R}^n\\
   \bar{q} = \left[ q_1, \ldots, q_n\right]^T \in \mathbb{R}^n

Kinematical Differential Equations
==================================

The variables :math:`u_1, \ldots, u_n` are defined as linear functions of the
time derivatives of the generalized coordinates :math:`\dot{q}_1, \ldots,
\dot{q}_n`. These variables are called generalized speeds. They take the form:

.. math::
   :label: generalized-speeds

   \bar{u} := \mathbf{Y}_k \dot{\bar{q}} + \bar{z}_k(\bar{q}, t)

:math:`\bar{u}` must be chosen such that :math:`\mathbf{Y}_k` is invertible.
Eq. :math:numref:`generalized-speeds` are called *kinematical differential
equations*. The most common, and always valid, choice of generalized speeds
is:

.. math::
   :label: generalized-speeds

   \bar{u} = \mathbf{I} \dot{\bar{q}}

where :math:`u_i = \dot{q}_i` for :math:`i=1,\ldots,n`.

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   me.init_vprinting(use_latex='mathjax')

.. jupyter-execute::

   psi, theta, phi = me.dynamicsymbols('psi, theta, phi')

   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')

   B.orient_body_fixed(A, (psi, theta, phi), 'ZXY')

   A_w_B = B.ang_vel_in(A).simplify()
   A_w_B

.. jupyter-execute::

   u1_ = A_w_B.dot(B.x)
   u2_ = A_w_B.dot(B.y)
   u3_ = A_w_B.dot(B.z)

   fk = sm.Matrix([u1_, u2_, u3_])

   qdot = [psi.diff(), theta.diff(), phi.diff()]

   Yk = fk.jacobian([psi.diff(), theta.diff(), phi.diff()])
   Yk

.. jupyter-execute::

   zk = fk.xreplace({psi.diff(): 0, theta.diff(): 0, phi.diff(): 0})
   zk

.. jupyter-execute::

   u1, u2, u3 = me.dynamicsymbols('u1, u2, u3')
   u = sm.Matrix([u1, u2, u3])

   sm.trigsimp(Yk.LUsolve(u - zk))

.. jupyter-execute::

   u1_ = A_w_B.dot(A.x)
   u2_ = A_w_B.dot(A.y)
   u3_ = A_w_B.dot(A.z)

   fk = sm.Matrix([u1_, u2_, u3_])

   qdot = [psi.diff(), theta.diff(), phi.diff()]

   Yk = sm.trigsimp(fk.jacobian(qdot))
   Yk

.. jupyter-execute::

   zk = fk.xreplace({psi.diff(): 0, theta.diff(): 0, phi.diff(): 0})
   zk

.. jupyter-execute::

   u1, u2, u3 = me.dynamicsymbols('u1, u2, u3')
   u = sm.Matrix([u1, u2, u3])

   sm.trigsimp(Yk.LUsolve(u - zk))

Chaplygin Sleigh
================

.. jupyter-execute::

   x, y, theta = me.dynamicsymbols('x, y, theta')

   N = me.ReferenceFrame('N')
   A = me.ReferenceFrame('A')

   A.orient_axis(N, theta, N.z)

   O = me.Point('O')
   P = me.Point('P')

   P.set_pos(O, x*N.x + y*N.y)

   O.set_vel(N, 0)

   P.vel(N).express(A)

.. jupyter-execute::

   fn = P.vel(N).dot(A.y)
   fn

.. todo:: Show how to prove it is nonholonomic.

Snake Board
===========

.. jupyter-execute::

   q1, q2, q3, q4, q5 = me.dynamicsymbols('q1, q2, q3, q4, q5')
   l = sm.symbols('l')

   N = me.ReferenceFrame('N')
   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')
   C = me.ReferenceFrame('C')

   A.orient_axis(N, q3, N.z)
   B.orient_axis(A, q4, A.z)
   C.orient_axis(A, q5, A.z)

   O = me.Point('O')
   Ao = me.Point('A_o')
   Bo = me.Point('B_o')
   Co = me.Point('C_o')

   Ao.set_pos(O, q1*N.x + q2*N.y)
   Bo.set_pos(Ao, l/2*A.x)
   Co.set_pos(Ao, -l/2*A.x)

   O.set_vel(N, 0)

   Bo.vel(N)

.. jupyter-execute::

   u1, u2, u3, u4, u5 = me.dynamicsmbols('u1, u2, u3, u4, u5')

   fn = sm.Matrix([Bo.vel(N).dot(B.y), Co.vel(N).dot(C.y)])
   sm.trigsimp(fn)



.. rubric:: Footnotes

.. [*] Well, we could find a very strong person to push th ecar sideways,
   overcoming the very high resisting friction force.
