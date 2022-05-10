===============================
Constrained Equations of Motion
===============================

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`constrained-eom` or Jupyter Notebook:
   :jupyter-download:notebook:`constrained-eom`.

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   me.init_vprinting(use_latex='mathjax')

When there are holonomic constraints present the equations of motion are
comprised of the kinematical differential equations, dynamical differential
equations, and the holonomic constraint equations. This set of equations are
differential algebraic equations, instead of ordinary differential equations.

N : number of coordinates
n : number of genereralized coordinates
M : number of configuratoin constraints
p : number of independent generalized speeds

Given $N$ coordinates where $n$ of those are independent generalized
coordinates, we cannot, in general, explicitly solve for the independent
coordinates. So we must formulate the kinematical and dynamical equations of
motion with $N$ coordinates.

q_s n indepdentdent generalized coordinates
q_r M dependent coordinates

q = [q_s, q_r]

.. math::

   \bar{q}, \bar{u} \in \mathbb{R}^N

.. math::
   :label: eq-holonomic-constrained-eom

   \bar{f}_d(\dot{\bar{u}}, \bar{u}, \bar{q}, t)  = 0 \\
   \bar{f}_k(\dot{\bar{q}}, \bar{u}, \bar{q}, t)  = 0 \\
   \bar{f}_h(\bar{q}, t) = 0


.. jupyter-execute::

   q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')
   u1, u2, u3 = me.dynamicsymbols('u1, u2, u3')
   la, lb, lc, ln = sm.symbols('l_a, l_b, l_c, l_n')
   m = sm.symbols('m')

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

   P4.pos_from(P1)
   r_P1_P4 = ln*N.x
   loop = P4.pos_from(P1) - r_P1_P4
   fh = sm.Matrix([loop.dot(N.x), loop.dot(N.y)])

   A.set_ang_vel(N, u1*N.z)
   B.set_ang_vel(A, u2*A.z)
   C.set_ang_vel(B, u3*B.z)

   P1.set_vel(N, 0)
   P2.v2pt_theory(P1, N, A)
   P3.v2pt_theory(P2, N, B)
   P4.v2pt_theory(P3, N, C)

   Rs_P2 = -m*P2.acc(N)
   Rs_P3 = -m*P3.acc(N)
   Rs_P4 = -m*P4.acc(N)

   Frs = sm.Matrix([
       P2.vel(N).diff(u1, N).dot(Rs_P2) + P3.vel(N).diff(u1, N).dot(Rs_P3) + P4.vel(N).diff(u1, N).dot(Rs_P4),
       P2.vel(N).diff(u2, N).dot(Rs_P2) + P3.vel(N).diff(u2, N).dot(Rs_P3) + P4.vel(N).diff(u2, N).dot(Rs_P4),
       P2.vel(N).diff(u3, N).dot(Rs_P2) + P3.vel(N).diff(u3, N).dot(Rs_P3) + P4.vel(N).diff(u3, N).dot(Rs_P4),
   ])

   Frs

.. jupyter-execute::

   t = me.dynamicsymbols._t
   sm.trigsimp(fh.diff(t, 2))
