=============================================
Bringing Noncontributing Forces into Evidence
=============================================

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`noncontributing` or Jupyter Notebook:
   :jupyter-download:notebook:`noncontributing`.

This is a simple double pendulum with two masses and each pendulum
section has the same length.

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   me.init_vprinting(use_latex='mathjax')

.. jupyter-execute::

   m1, m2, l1, l2, g = sm.symbols('m1, m2, l1, l2, g')
   q1, q2, u1, u2 = me.dynamicsymbols('q1, q2, u1, u2')
   TP1, TP2 = me.dynamicsymbols('T_{P_1}, T_{P_2}')
   T1, T2 = me.dynamicsymbols('T1, T2')
   t = me.dynamicsymbols._t

   q = sm.Matrix([q1, q2])
   u = sm.Matrix([u1, u2])
   r = sm.Matrix([TP1, TP2])
   ud = u.diff(t)

.. jupyter-execute::

   N = me.ReferenceFrame('N')
   A = N.orientnew('A', 'Axis', (q1, N.z))
   B = N.orientnew('B', 'Axis', (q2, N.z))

.. jupyter-execute::

   A.set_ang_vel(N, u1*N.z)
   B.set_ang_vel(N, u2*N.z)

.. jupyter-execute::

   O = me.Point('O')
   P1 = O.locatenew('P1', -l1*A.y)
   P2 = P1.locatenew('P2', -l2*B.y)

.. jupyter-execute::

   O.set_vel(N, 0)
   P1.v2pt_theory(O, N, A)
   P2.v2pt_theory(P1, N, B)

Verify answer using Newton-Euler equations
==========================================

.. jupyter-execute::

   R_P1 = (TP2*sm.sin(q2) - TP1*sm.sin(q1))*N.x + (TP1*sm.cos(q1)-TP2*sm.cos(q2) - m1*g)*N.y
   R_P2 = (-TP2*sm.sin(q2)*N.x + (TP2*sm.cos(q2)-m2*g)*N.y)

   R_P1

.. jupyter-execute::

   R_P1 = -T1*A.y + T2*B.y - m1*g*N.y
   R_P1.express(N).simplify()

.. jupyter-execute::

   R_P2 = -T2*B.y - m2*g*N.y
   R_P2.express(N).simplify()

.. jupyter-execute::

   veq1 = -m1*P1.acc(N) + R_P1
   veq2 = -m2*P2.acc(N) + R_P2

.. jupyter-execute::

   scalar_eqs = sm.Matrix([
       veq1.dot(N.x),
       veq1.dot(N.y),
       veq2.dot(N.x),
       veq2.dot(N.y),
   ])

.. jupyter-execute::

   newton = [u1.diff(), u2.diff(), TP1, TP2]
   newton_zero = {v: 0 for v in newton}

   M = scalar_eqs.jacobian(newton)
   g = scalar_eqs.xreplace(newton_zero)

.. jupyter-execute::

   newton_sol = M.LUsolve(g)
   newton_sol = sm.trigsimp(newton_sol)
   newton_sol

Introduce fictitious generalized speeds that correspond to components of desired forces and torques
===================================================================================================

Here I introduce the fictitious generalized speed u3 that lets the
particle P1 have a “separation velocity” relative to its fixed location
on the pendulum arm. This is aligned with the desired non-contributing
tension force we want to bring into evidence.

.. jupyter-execute::

   u3, u4 = me.dynamicsymbols('u3, u4')

   P1.set_vel(N, P1.vel(N) + u3*A.y)
   P1.vel(N)

.. jupyter-execute::

   P2.v2pt_theory(P1, N, B)

Add a similar fictitious generalized speed u4 for the second tension
force.

.. jupyter-execute::

   P2.set_vel(N, P2.vel(N) + u4*B.y)
   P2.vel(N)

.. jupyter-execute::

   P1.acc(N)

.. jupyter-execute::

   P2.acc(N)

Introduce unknown force and torques into the resultants
=======================================================

These are the two time varying tension forces we want to bring into
evidence:

For u1 and u2, we use the resultant with only the original contributing
forces.

.. jupyter-execute::

   RP1 = -m1*g*N.y
   RP2 = -m2*g*N.y

For the particle we need to add the non-contributing forces that
correspond to u3 and u4

.. jupyter-execute::

   RP1_aux = RP1 + TP1*A.y
   RP2_aux = RP1 + TP2*B.y

We also need equal and opposite tension forces acting back on the
pendulum arm (but not the force due to gravity):

.. jupyter-execute::

   RP1_aux_neg = -TP1*A.y
   RP2_aux_neg = -TP2*B.y

GAF
===

Calculate the two GAFs for the the real genearlized speeds as normal:

.. jupyter-execute::

   F1 = P1.vel(N).diff(u1, N).dot(RP1) + P2.vel(N).diff(u1, N).dot(RP2)
   F2 = P1.vel(N).diff(u2, N).dot(RP1) + P2.vel(N).diff(u2, N).dot(RP2)

For F3 and F4 you need to use the resultaants that include the tension
forces and they need to be associated with the appropriate velocities
for the equal and opposite forces.

.. jupyter-execute::

   F3 = (P1.vel(N).diff(u3, N).dot(RP1 + TP1*A.y) +  # velocity of the particle which includes u3 
         (P1.vel(N) - u3*A.y).diff(u3, N).dot(-TP1*A.y) +  # velocity of the tip of the pendulum arm (does not include u3)
         P2.vel(N).diff(u3, N).dot(RP2 + TP2*B.y) +  # velocity of the second particle which includes u3 and u4
         (P2.vel(N) - u4*B.y).diff(u3, N).dot(-TP2*B.y))  # velocity of the tip of the second pendulum arm (includes u3 but not u4)
   F3

.. jupyter-execute::

   F4 = (P1.vel(N).diff(u4, N).dot(RP1 + TP1*A.y) +
         (P1.vel(N) - u3*A.y).diff(u4, N).dot(-TP1*A.y) +
         P2.vel(N).diff(u4, N).dot(RP2 + TP2*B.y) +
         (P2.vel(N) - u4*B.y).diff(u4, N).dot(-TP2*B.y))
   F4

GIF
===

Calculate all GIFs with u1, u2, u3, and u4 present in the velocities and
accelerations.

.. jupyter-execute::

   F1s = P1.vel(N).diff(u1, N).dot(-m1*P1.acc(N)) + P2.vel(N).diff(u1, N).dot(-m2*P2.acc(N))
   F2s = P1.vel(N).diff(u2, N).dot(-m1*P1.acc(N)) + P2.vel(N).diff(u2, N).dot(-m2*P2.acc(N))
   F3s = P1.vel(N).diff(u3, N).dot(-m1*P1.acc(N)) + P2.vel(N).diff(u3, N).dot(-m2*P2.acc(N))
   F4s = P1.vel(N).diff(u4, N).dot(-m1*P1.acc(N)) + P2.vel(N).diff(u4, N).dot(-m2*P2.acc(N))

   F1s

.. jupyter-execute::

   F2s

.. jupyter-execute::

   F3s

.. jupyter-execute::

   F4s

Kane’s Equations
================

.. jupyter-execute::

   k1 = F1 + F1s
   k2 = F2 + F2s
   k3 = F3 + F3s
   k4 = F4 + F4s

Substitute zero for all fictitious quantities
=============================================

.. jupyter-execute::

   k1_ = k1.subs({u3.diff(): 0, u4.diff(): 0, u3: 0, u4: 0})
   k2_ = k2.subs({u3.diff(): 0, u4.diff(): 0, u3: 0, u4: 0})
   k3_ = k3.subs({u3.diff(): 0, u4.diff(): 0, u3: 0, u4: 0})
   k4_ = k4.subs({u3.diff(): 0, u4.diff(): 0, u3: 0, u4: 0})

.. jupyter-execute::

   kanes = [k1_, k2_, k3_, k4_]

.. jupyter-execute::

   me.find_dynamicsymbols(k1_)

.. jupyter-execute::

   me.find_dynamicsymbols(k2_)

.. jupyter-execute::

   me.find_dynamicsymbols(k3_)

.. jupyter-execute::

   me.find_dynamicsymbols(k4_)

.. jupyter-execute::

   k1_

.. jupyter-execute::

   k2_

.. jupyter-execute::

   k3_

.. jupyter-execute::

   k4_

Solve for all unknowns
======================

.. jupyter-execute::

   sol = sm.solve(kanes, u1.diff(), u2.diff(), TP1, TP2)

.. jupyter-execute::

   sol[u1.diff()].simplify()

.. jupyter-execute::

   sol[u2.diff()].simplify()

.. jupyter-execute::

   sol[TP1].simplify()

.. jupyter-execute::

   sol[TP2].simplify()

.. jupyter-execute::

   TP1_sol = sol[TP1].simplify()

.. jupyter-execute::

   me.find_dynamicsymbols(TP1_sol)

.. jupyter-execute::

   TP1_sol.free_symbols

Evaluate the force expressions with arrays
==========================================

And compare their results numerically.

.. jupyter-execute::

   eval_TP1 = sm.lambdify((q1, q2, u1, u2, m1, m2, g, l), TP1_sol)

.. jupyter-execute::

   import numpy as np

.. jupyter-execute::

   times = np.linspace(0, 10)
   omega = 0.2

.. jupyter-execute::

   q1_vals = np.sin(omega*times)

.. jupyter-execute::

   u1_vals = omega*np.cos(omega*times)

.. jupyter-execute::

   vals = eval_TP1(q1_vals, q1_vals, u1_vals, u1_vals, 1.0, 1.0, 9.81, 2.0)

.. jupyter-execute::

   me.find_dynamicsymbols(newton_sol[TP1])

.. jupyter-execute::

   eval_TP1_newton = sm.lambdify((q1, q2, u1, u2, m1, m2, g, l), newton_sol[TP1])
   vals_newton = eval_TP1_newton(q1_vals, q1_vals, u1_vals, u1_vals, 1.0, 1.0, 9.81, 2.0)

.. jupyter-execute::

   import matplotlib.pyplot as plt

.. jupyter-execute::

   plt.plot(times, vals, times, vals_newton, '.')

.. jupyter-execute::

   func = lambda x, y: x + y

.. jupyter-execute::

   func(1, 2)

.. jupyter-execute::

   generate_numeric_func = sm.lambdify

.. jupyter-execute::

   eval_TP1 = generate_numeric_func((q1, u1, m1, g, l), TP1_sol)
