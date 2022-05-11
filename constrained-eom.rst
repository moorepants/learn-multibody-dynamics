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

Nonholomic Constraints
======================

Nonholonomic constraint equations are linear in the independent and dependent
generalized speeds, by definition. We have shown that you can explicitly solve
for the independent generalized speeds :math:`\bar{u}_s` as a function of the
dependent generalized speeds :math:`\bar{u}_r`. This means that number of
dynamical differential equations can be reduced to :math:`p` from :math:`n`.
The nonholonomic constraints :math:`\bar{f}_n` and its time derivative
:math:`\dot{\bar{f}}_n` can be used to solve for the dependent generalized
speeds and their time derivatives.

.. math::
   :label: eq-nonholonomic-eom

   \bar{f}_d(\dot{\bar{u}}_s, \bar{u}, \bar{q}, t) = \mathbf{M}_d\dot{\bar{u}}_s + \bar{g}_d = 0 \\
   \bar{f}_k(\dot{\bar{q}}, \bar{u}, \bar{q}, t) = \mathbf{M}_k\dot{\bar{q}} + \bar{g}_k  = 0 \\
   \bar{f}_n(\bar{u}_s, \bar{u}_r, \bar{q}, t) = \mathbf{M}_n\bar{u}_r + \bar{g}_n = 0 \\
   \dot{\bar{f}}_n(\dot{\bar{u}}_s, \dot{\bar{u}}_r, \bar{u}_s, \bar{u}_r, \bar{q}, t) = \mathbf{M}_{nd}\dot{\bar{u}}_r + \bar{g}_{nd}= 0

The following psuedo code shows how the derivatives of the generalized speeds
and generalized coordinates can be calculated from the above equations.

.. code:: python

   def eval_rhs(t, x, p):
      """
      t : float
      x : array_like, shape(2*n,)
      p : array_like
      """

      q = x[:n]  # n generalized coordinates
      us = x[n:n+p]  # p independent generalized speeds
      ur = x[n+p:]  # m dependent generalized speeds

      # solve fn for the dependent generalized speeds using the constraints
      Mn, gn = eval_n(q, us, p)
      ur = solve(Mn, -gn)

      # solve for the n qdots
      Mk, gk = eval_k(q, us, ur)
      qd = solve(Mk, -gk)

      # solve for the p us_dots
      Md, gd = eval_d(q, us, ur)
      usd = solve(Md, -gd)

      # solve for the m ur_dots
      Mnd, gnd = eval_nd(q, us, ur, usd)
      urd = solve(Mnd, -gnd)

      return np.hstack((qd, usd, urd))

Let's revisit the snakeboard example and develop the equations of motion for
that nonholomoic system. For simplicity we will assume that the mass and
moments of inertia of the three bodies are the same.

.. jupyter-execute::

   q1, q2, q3, q4, q5 = me.dynamicsymbols('q1, q2, q3, q4, q5')
   l = sm.symbols('l')

The reference frames are all simple rotations about the axis normal to the
plane:

.. jupyter-execute::

   N = me.ReferenceFrame('N')
   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')
   C = me.ReferenceFrame('C')

   A.orient_axis(N, q3, N.z)
   B.orient_axis(A, q4, A.z)
   C.orient_axis(A, q5, A.z)

.. jupyter-execute::

   A.ang_vel_in(N)
   B.ang_vel_in(N)
   C.ang_vel_in(N)


.. jupyter-execute::

   O = me.Point('O')
   Ao = me.Point('A_o')
   Bo = me.Point('B_o')
   Co = me.Point('C_o')

   Ao.set_pos(O, q1*N.x + q2*N.y)
   Bo.set_pos(Ao, l/2*A.x)
   Co.set_pos(Ao, -l/2*A.x)

.. jupyter-execute::

   O.set_vel(N, 0)
   Ao.vel(N)

.. jupyter-execute::

   Bo.v2pt_theory(Ao, N, A)

.. jupyter-execute::

   Co.v2pt_theory(Ao, N, A)

.. jupyter-execute::

   fn = sm.Matrix([Bo.vel(N).dot(B.y),
                   Co.vel(N).dot(C.y)])
   fn = sm.trigsimp(fn)
   fn

.. jupyter-execute::

   u1, u2, u3, u4, u5 = me.dynamicsymbols('u1, u2, u3, u4, u5')

   u_repl = {
       q1.diff(): u1,
       q2.diff(): u2,
       l*q3.diff()/2: u3,
       q4.diff(): u4,
       q5.diff(): u5
   }

   fn = fn.subs(u_repl)
   fn

.. jupyter-execute::

   q = sm.Matrix([q1, q2, q3, q4, q5])
   qd = q.diff()
   qd_zero = {qdi: 0 for qdi in qd}

   fk = sm.Matrix([rhs - lhs for lhs, rhs in u_repl.items()])
   Mk = fk.jacobian(qd)
   gk = fk.xreplace(qd_zero)
   Mk, gk

.. jupyter-execute::

   solk = Mk.LUsolve(-gk)
   solk

.. jupyter-execute::

   qd_repl = {qdi: solki for qdi, solki in zip(qd, solk)}
   qd_repl

.. jupyter-execute::

   us = sm.Matrix([u3, u4, u5])
   ur = sm.Matrix([u1, u2])

   u = us.col_join(ur)

   ur_zero = {ui: 0 for ui in ur}
   us_zero = {ui: 0 for ui in us}
   u_zero = {ui: 0 for ui in u}

.. jupyter-execute::

   Mn = fn.jacobian(ur)
   gn = fn.xreplace(ur_zero)
   Mn, gn

.. jupyter-execute::

   soln = Mn.LUsolve(-gn)
   soln

.. jupyter-execute::

   t = me.dynamicsymbols._t
   u1_dep = sm.Function('u1')(u3, u4, u5, t)
   u2_dep = sm.Function('u2')(u3, u4, u5, t)
   u1_dep = soln[0]
   u2_dep = soln[1]
   u1_dep.diff(t)

.. jupyter-execute::

   N_w_A = A.ang_vel_in(N).xreplace(qd_repl).xreplace({u1: u1_dep, u2: u2_dep})

.. jupyter-execute::

   N_w_B = B.ang_vel_in(N).xreplace(qd_repl).xreplace({u1: u1_dep, u2: u2_dep})

.. jupyter-execute::

   N_w_C = C.ang_vel_in(N).xreplace(qd_repl).xreplace({u1: u1_dep, u2: u2_dep})

.. jupyter-execute::

   N_v_Ao = Ao.vel(N).xreplace(qd_repl).xreplace({u1: u1_dep, u2: u2_dep})

.. jupyter-execute::

   N_v_Bo = Bo.vel(N).xreplace(qd_repl).xreplace({u1: u1_dep, u2: u2_dep})

.. jupyter-execute::

   N_v_Co = Co.vel(N).xreplace(qd_repl).xreplace({u1: u1_dep, u2: u2_dep})

.. jupyter-execute::

   vels = (N_w_A, N_w_B, N_w_C, N_v_Ao, N_v_Bo, N_v_Co)

   me.partial_velocity(vels, (u3, u4, u5), N)

Holonomic Constraints
=====================

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
   m, g = sm.symbols('m, g')

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

.. jupyter-execute::

   t = me.dynamicsymbols._t
   qd_repl = {q1.diff(t): u1, q2.diff(t): u2, q3.diff(t): u3}
   fhd = fh.diff(t).xreplace(qd_repl)
   me.find_dynamicsymbols(fhd)

.. jupyter-execute::

   res = sm.solve(fhd, u2, u3)
   #{k: sm.trigsimp(v) for k, v in res.items()}

.. jupyter-execute::

   fhdd = fhd.diff(t).xreplace(qd_repl)
   me.find_dynamicsymbols(fhdd)

.. jupyter-execute::

   A.set_ang_vel(N, u1*N.z)
   B.set_ang_vel(A, res[u2]*A.z)
   C.set_ang_vel(B, res[u3]*B.z)

   P1.set_vel(N, 0)
   P2.v2pt_theory(P1, N, A)
   P3.v2pt_theory(P2, N, B)
   P4.v2pt_theory(P3, N, C)

   R_P2 = -m*g*N.y
   R_P3 = -m*g*N.y
   R_P4 = -m*g*N.y

   Fr = sm.Matrix([
       P2.vel(N).diff(u1, N).dot(R_P2) +
       P3.vel(N).diff(u1, N).dot(R_P3) +
       P4.vel(N).diff(u1, N).dot(R_P4),
   ])

   me.find_dynamicsymbols(Fr)

.. jupyter-execute::

   me.find_dynamicsymbols(P2.acc(N).to_matrix(N))

.. jupyter-execute::

   me.find_dynamicsymbols(P3.acc(N).to_matrix(N))

.. jupyter-execute::

   me.find_dynamicsymbols(P4.acc(N).to_matrix(N))

.. jupyter-execute::

   Rs_P2 = -m*P2.acc(N)
   Rs_P3 = -m*P3.acc(N).xreplace(qd_repl).xreplace(res)
   Rs_P4 = -m*P4.acc(N).xreplace(qd_repl).xreplace(res)

   Frs = sm.Matrix([
       P2.vel(N).diff(u1, N).dot(Rs_P2) +
       P3.vel(N).diff(u1, N).dot(Rs_P3) +
       P4.vel(N).diff(u1, N).dot(Rs_P4),
   ])
   me.find_dynamicsymbols(Frs)

.. jupyter-execute::

   q = sm.Matrix([q1, q2, q3])
   u = sm.Matrix([u1])
   p = sm.Matrix([la, lb, lc, ln, m, g])

.. jupyter-execute::

   Md = Frs.jacobian([u1.diff()])
   gd = Frs.xreplace({u1.diff(): 0}) + Fr

.. jupyter-execute::

   eval_Mdgd = sm.lambdify((q, u, p), (Md, gd))
   eval_fh = sm.lambdify((q, p), fh)

.. jupyter-execute::

   import numpy as np

   p_vals = np.array([
       0.8,
       2.0,
       1.0,
       2.0,
       1.0,
       9.81,
   ])

   q1_val = np.deg2rad(10.0)

   from scipy.optimize import fsolve

   eval_fh_fsolve = lambda x, q1, p: np.squeeze(eval_fh(np.hstack((q1, x)), p))

   q2_val, q3_val = fsolve(eval_fh_fsolve, np.deg2rad([-6.0, 150]), args=(q1_val, p_vals))

   q_vals = np.array([q1_val, q2_val, q3_val])
   np.rad2deg(q_vals)

.. jupyter-execute::

   eval_u2u3 = sm.lambdify((q, u, p), (res[u2], res[u3]))
   eval_u2u3(q_vals, [1.0], p_vals)

.. jupyter-execute::


   def eval_rhs(t, x, p):

       q1, q2, q3, u1 = x

       q1d = u1

       u2, u3 = eval_u2u3([q1, q2, q3], [u1], p)

       q2d = u2
       q3d = u3

       Md, gd = eval_Mdgd([q1, q2, q3], [u1], p)

       u1d = -Md[0]/gd[0]

       return np.array([q1d, q2d, q3d, u1d[0]])

.. jupyter-execute::

   x0 = np.hstack((q_vals, 0.1))

   eval_rhs(0.0, x0, p_vals)

.. jupyter-execute::

   from scipy.integrate import solve_ivp

   sol = solve_ivp(eval_rhs, (0.0, 5.0), x0, args=(p_vals,))

   q1_traj, q2_traj, q3_traj, u1_traj = sol.y

   constraint_violation = eval_fh((q1_traj, q2_traj, q3_traj), p_vals)

.. jupyter-execute::

   import matplotlib.pyplot as plt
   plt.plot(sol.t, sol.y.T)
   plt.legend(['q1', 'q2', 'q3', 'u1'])

.. jupyter-execute::

   plt.plot(sol.t, np.squeeze(constraint_violation).T)
