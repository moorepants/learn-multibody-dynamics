===============================
Constrained Equations of Motion
===============================

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`constrained-eom` or Jupyter Notebook:
   :jupyter-download:notebook:`constrained-eom`.

.. jupyter-execute::

   import numpy as np
   import sympy as sm
   import sympy.physics.mechanics as me
   me.init_vprinting(use_latex='mathjax')

We introduced two types of constraints holonomic (configuration) constraints
and nonholonomic (motion) constraints in prior chapters. In general, holonomic
constraints are nonlinear constraints in the coordinates.

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

.. math::
   :label: eq-nonholonomic-steps

   \bar{u} = [\bar{u}_s \ \bar{u}_r] \\
   \bar{u}_r = -\mathbf{M}_n^{-1}(\bar{q}, t)\bar{g}_n(\bar{u}_s, \bar{q}, t) \\
   \dot{\bar{q}} = -\mathbf{M}_k^{-1}(\bar{q}, t)\bar{g}_k(\bar{u}, \bar{q}, t) \\
   \dot{\bar{u}}_s = -\mathbf{M}_d(\bar{q}, t) \bar{g}_d(\bar{u}, \bar{q}, t) \\
   \dot{\bar{u}}_r = -\mathbf{M}_{nd}^{-1}(\dot{\bar{u}}_s, \bar{q}, t) \bar{g}_{nd}(\bar{u}, \bar{q}, t)

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
   u1, u2, u3, u4, u5 = me.dynamicsymbols('u1, u2, u3, u4, u5')
   l, I, m = sm.symbols('l, I, m')
   t = me.dynamicsymbols._t

   p = sm.Matrix([l, I, m])
   q = sm.Matrix([q1, q2, q3, q4, q5])  # coordinates
   us = sm.Matrix([u3, u4, u5])  # independent
   ur = sm.Matrix([u1, u2])  # dependent
   u = us.col_join(ur)

   p, q, us, ur, u

.. jupyter-execute::

   qd = q.diff()
   ud = u.diff(t)
   usd = us.diff(t)
   urd = ur.diff(t)

   qd_zero = {qdi: 0 for qdi in qd}
   ur_zero = {ui: 0 for ui in ur}
   us_zero = {ui: 0 for ui in us}
   usd_zero = {udi: 0 for udi in usd}
   urd_zero = {udi: 0 for udi in urd}

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

   A.ang_vel_in(N)
   B.ang_vel_in(N)
   C.ang_vel_in(N)

   O = me.Point('O')
   Ao = me.Point('A_o')
   Bo = me.Point('B_o')
   Co = me.Point('C_o')

   Ao.set_pos(O, q1*N.x + q2*N.y)
   Bo.set_pos(Ao, l/2*A.x)
   Co.set_pos(Ao, -l/2*A.x)

   O.set_vel(N, 0)
   Ao.vel(N)
   Bo.v2pt_theory(Ao, N, A)
   Co.v2pt_theory(Ao, N, A)

The coordinates (independent and dependent) may be present in all of the
equations.

.. jupyter-execute::

   # n=5 kinematical differential equations
   fk = sm.Matrix([
      u1 - q1.diff(t),
      u2 - q2.diff(t),
      u3 - l*q3.diff(t)/2,
      u4 - q4.diff(t),
      u5 - q5.diff(t),
   ])
   # kinematical differential equation linear coefficients
   Mk = fk.jacobian(qd)
   gk = fk.xreplace(qd_zero)
   eval_k = sm.lambdify((q, u, p), (Mk, gk))

.. jupyter-execute::

   # solve the kinematical differential equations symbollically for substitution
   qd_sol = Mk.LUsolve(-gk)
   qd_repl = dict(zip(qd, qd_sol))
   qd_repl

.. jupyter-execute::

   fn = sm.Matrix([Bo.vel(N).dot(B.y),
                   Co.vel(N).dot(C.y)])
   fn

.. jupyter-execute::

   fn = fn.xreplace(qd_repl)
   fn

.. jupyter-execute::

   # nonholonomic constraints linear coefficients
   Mn = fn.jacobian(ur)
   gn = fn.xreplace(ur_zero)
   Mn, gn

.. jupyter-execute::

   eval_n = sm.lambdify((us, q, p), (Mn, gn))

.. jupyter-execute::

   # solve the nonholonomic constraints for the dependent generalized speeds ur
   ur_sol = Mn.LUsolve(-gn)
   ur_repl = dict(zip(ur, ur_sol))

.. jupyter-execute::

   N_w_A = A.ang_vel_in(N).xreplace(qd_repl).xreplace(ur_repl)
   N_w_B = B.ang_vel_in(N).xreplace(qd_repl).xreplace(ur_repl)
   N_w_C = C.ang_vel_in(N).xreplace(qd_repl).xreplace(ur_repl)
   N_v_Ao = Ao.vel(N).xreplace(qd_repl).xreplace(ur_repl)
   N_v_Bo = Bo.vel(N).xreplace(qd_repl).xreplace(ur_repl)
   N_v_Co = Co.vel(N).xreplace(qd_repl).xreplace(ur_repl)

.. jupyter-execute::

   vels = (N_w_A, N_w_B, N_w_C, N_v_Ao, N_v_Bo, N_v_Co)
   w_A, w_B, w_C, v_Ao, v_Bo, v_Co = me.partial_velocity(vels, us, N)

.. jupyter-execute::

   fnd = fn.diff(t).xreplace(qd_repl)
   Mnd = fnd.jacobian(urd)
   gnd = fnd.xreplace(urd_zero).xreplace(ur_repl)
   usd_dummy = sm.Matrix([sm.Dummy('u3d'), sm.Dummy('u4d'), sm.Dummy('u5d')])
   usd_dummy_repl = dict(zip(usd, usd_dummy))
   eval_nd = sm.lambdify((usd_dummy, u, q, p), (Mnd, gnd.xreplace(usd_dummy_repl)))
   urd_sol = Mnd.LUsolve(-gnd)
   urd_repl = dict(zip(urd, urd_sol))

   qdd_repl = {k.diff(t): v.diff(t) for k, v in qd_repl.items()}
   qdd_repl

.. jupyter-execute::

   Rs_Ao = -m*Ao.acc(N).xreplace(qdd_repl).xreplace(qd_repl).xreplace(urd_repl)
   Rs_Bo = -m*Bo.acc(N).xreplace(qdd_repl).xreplace(qd_repl).xreplace(urd_repl)
   Rs_Co = -m*Co.acc(N).xreplace(qdd_repl).xreplace(qd_repl).xreplace(urd_repl)

   me.find_dynamicsymbols(Rs_Bo, reference_frame=N)

.. jupyter-execute::

   I_A_Ao = me.inertia(A, 0, 0, I)
   I_B_Bo = me.inertia(B, 0, 0, I)
   I_C_Co = me.inertia(C, 0, 0, I)

.. jupyter-execute::

   Ts_A = -A.ang_acc_in(N).dot(I_A_Ao).xreplace(qdd_repl).xreplace(qd_repl).xreplace(urd_repl)
   Ts_B = -B.ang_acc_in(N).dot(I_B_Bo).xreplace(qdd_repl).xreplace(qd_repl).xreplace(urd_repl)
   Ts_C = -C.ang_acc_in(N).dot(I_C_Co).xreplace(qdd_repl).xreplace(qd_repl).xreplace(urd_repl)

.. jupyter-execute::

   F3s = (v_Ao[0].dot(Rs_Ao) + v_Bo[0].dot(Rs_Bo) + v_Co[0].dot(Rs_Co) +
          w_A[0].dot(Ts_A) + w_B[0].dot(Ts_B) + w_C[0].dot(Ts_C))
   F4s = (v_Ao[1].dot(Rs_Ao) + v_Bo[1].dot(Rs_Bo) + v_Co[1].dot(Rs_Co) +
          w_A[1].dot(Ts_A) + w_B[1].dot(Ts_B) + w_C[1].dot(Ts_C))
   F5s = (v_Ao[2].dot(Rs_Ao) + v_Bo[2].dot(Rs_Bo) + v_Co[2].dot(Rs_Co) +
          w_A[2].dot(Ts_A) + w_B[2].dot(Ts_B) + w_C[2].dot(Ts_C))

.. jupyter-execute::

   Frs = sm.Matrix([F3s, F4s, F5s])
   Md = Frs.jacobian(usd)
   gd = Frs.xreplace(usd_zero)

   me.find_dynamicsymbols(Frs)

.. jupyter-execute::

   Md

.. jupyter-execute::

   me.find_dynamicsymbols(Md)

.. jupyter-execute::

   me.find_dynamicsymbols(gd)

.. jupyter-execute::

   eval_d = sm.lambdify((q, us, p), (Md, gd))

.. jupyter-execute::

   def eval_rhs(t, x, p):
      # x = [q1, q2, q3, q4, q5, u3, u4, u5, u1, u2]
      q = x[:5]
      us = x[5:8]

      Mn, gn = eval_n(us, q, p)
      ur = np.linalg.solve(Mn, -gn.squeeze())

      u = np.hstack((us, ur))

      Mk, gk = eval_k(q, u, p)
      qd = np.linalg.solve(Mk, -gk.squeeze())

      Md, gd = eval_d(q, us, p)
      usd = np.linalg.solve(Md, -gd.squeeze())

      Mnd, gnd = eval_nd(usd, u, q, p)
      urd = np.linalg.solve(Mnd, -gnd.squeeze())

      return np.hstack((qd, usd, urd))

   print(eval_rhs(1.0, np.random.random(10), np.random.random(3)))

.. jupyter-execute::

   p_vals = np.array([
       0.3,  # l [m]
       0.1,  # I [kg*m^2]
       1.0,  # m [kg]
   ])

   q0 = np.array([
       0.0,  # q1 [m]
       0.0,  # q2 [m]
       0.0,  # q3 [rad]
       np.deg2rad(5.0),  # q4 [rad]
       -np.deg2rad(5.0),  # q5 [rad]
   ])

   us0 = np.array([
       0.01,  # u3 [m/s]
       0.01,  # u4 [rad/s]
       -0.01,  # u5 [rad/s]
   ])

   Mn_vals, gn_vals = eval_n(us0, q0, p_vals)
   ur0 = np.linalg.solve(Mn_vals, -gn_vals.squeeze())

   x0 = np.hstack((q0, us0, ur0))

   from scipy.integrate import solve_ivp

   t0, tf = 0.0, 50.0

   ts = np.linspace(t0, tf, num=1001)

   sol = solve_ivp(eval_rhs, (ts[0], ts[-1]), x0, args=(p_vals,), t_eval=ts)

.. jupyter-execute::

   import matplotlib.pyplot as plt

   fig, axes = plt.subplots(2, 1, sharex=True)

   axes[0].plot(sol.t, np.rad2deg(sol.y[:3]).T)
   axes[1].plot(sol.t, sol.y[3:5].T)

.. jupyter-execute::

   eval_fn = sm.lambdify((q, u, p), fn)

   con_violations = eval_fn(sol.y[:5], sol.y[5:], p_vals).squeeze()

   fig, ax = plt.subplots()
   ax.plot(sol.t, con_violations.T)

.. jupyter-execute::

   Cl, Cr, Bl, Br = sm.symbols('C_l, C_r, B_l, B_r', cls=me.Point)
   Cl.set_pos(Co, -l/4*C.y)
   Cr.set_pos(Co, l/4*C.y)
   Bl.set_pos(Bo, -l/4*B.y)
   Br.set_pos(Bo, l/4*B.y)

   coordinates = Cl.pos_from(O).to_matrix(N)
   for point in [Co, Cr, Co, Ao, Bo, Bl, Br]:
      coordinates = coordinates.row_join(point.pos_from(O).to_matrix(N))

   eval_point_coords = sm.lambdify((q, p), coordinates)
   eval_point_coords(q0, p_vals)

.. jupyter-execute::

   x, y, z = eval_point_coords(q0, p_vals)

   fig, ax = plt.subplots()

   line_prop = {
      'color': 'black',
      'marker': 'o',
      'markerfacecolor': 'blue',
      'markersize': 10,
   }

   lines, = ax.plot(x, y, **line_prop)
   ax.set_xlim((np.min(sol.y[0]) - 0.5, np.max(sol.y[0]) + 0.5))
   ax.set_ylim((np.min(sol.y[1]) - 0.5, np.max(sol.y[1]) + 0.5))

.. jupyter-execute::

   from matplotlib.animation import FuncAnimation

   def animate(i):
       x, y, z = eval_point_coords(sol.y.T[i, :5], p_vals)
       lines.set_data(x, y)

   ani = FuncAnimation(fig, animate, len(sol.t))

   from IPython.display import HTML

   HTML(ani.to_jshtml(fps=30))

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
