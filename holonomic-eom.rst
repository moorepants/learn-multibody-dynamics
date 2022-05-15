==============================================
Equations of Motion with Holonomic Constraints
==============================================

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`holonomic-eom` or Jupyter Notebook:
   :jupyter-download:notebook:`holonomic-eom`.

.. jupyter-execute::

   from IPython.display import HTML
   from matplotlib.animation import FuncAnimation
   from scipy.integrate import solve_ivp
   import matplotlib.pyplot as plt
   import numpy as np
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
   :label: eq-holonomic-constraints

   \bar{f}_h(\bar{q}, \bar{q}_r, t) = 0 \in \mathbb{R}^M

.. math::
   :label: eq-holonomic-constraints-dot

   \dot{\bar{f}}_h(\bar{u}, \bar{u}_r \bar{q}, \bar{q}_r, t) =
   \mathbf{M}_{hd}\bar{u}_r + \bar{g}_{hd} = 0 \in \mathbb{R}^M

.. math::

   \bar{u}_r = -\mathbf{M}_{hd}^{-1} \bar{g}_{hd} \in \mathbb{R}^M

.. math::
   :label: eq-holonomic-constrained-eom

   \bar{f}_k(\dot{\bar{q}}, \bar{u}, \bar{q}, \bar{q}_r, t)  = 0 \in \mathbb{R}^N \\
   \bar{f}_d(\dot{\bar{u}}, \bar{u}, \bar{q}, \bar{q}_r, t)  = 0 \in \mathbb{R}^n \\

Four-bar Linkage Equations of Motion
====================================

.. figure:: figures/configuration-four-bar.svg
   :align: center
   :width: 600px

   a) Shows four links in a plane :math:`A`, :math:`B`, :math:`C`, and
   :math:`N` with respective lengths :math:`l_a,l_b,l_c,l_n` connected in a
   closed loop at points :math:`P_1,P_2,P_3,P_4`. b) Shows the same linkage
   that has been seperated at point :math:`P_4` to make it an open chain of
   links.

.. jupyter-execute::

   q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')
   u1, u2, u3 = me.dynamicsymbols('u1, u2, u3')
   la, lb, lc, ln = sm.symbols('l_a, l_b, l_c, l_n')
   m, g = sm.symbols('m, g')
   t = me.dynamicsymbols._t

   q = sm.Matrix([q1])
   qr = sm.Matrix([q2, q3])
   qN = q.col_join(qr)
   u = sm.Matrix([u1])
   ur = sm.Matrix([u2, u3])
   uN = u.col_join(ur)

   qdN = qN.diff(t)

   qdN_zero = {qdi: 0 for qdi in qdN}
   uN_zero = {ui: 0 for ui in uN}

   p = sm.Matrix([la, lb, lc, ln, m, g])

.. jupyter-execute::

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

:math:`N=3`, :math:`M=2`, and :math:`n=1`.

.. jupyter-execute::

   loop = P4.pos_from(P1) - ln*N.x

   fh = sm.Matrix([loop.dot(N.x), loop.dot(N.y)])
   fh = sm.trigsimp(fh)
   fh

.. jupyter-execute::

   me.find_dynamicsymbols(fh)

.. jupyter-execute::

   fk = sm.Matrix([
       q1.diff(t) - u1,
       q2.diff(t) - u2,
       q3.diff(t) - u3,
   ])
   Mk = fk.jacobian(qdN)
   gk = fk.xreplace(qdN_zero)
   qdN_sol = -Mk.LUsolve(gk)
   qd_repl = dict(zip(qdN, qdN_sol))
   qd_repl

.. jupyter-execute::

   fhd = fh.diff(t).xreplace(qd_repl)
   me.find_dynamicsymbols(fhd)

.. jupyter-execute::

   ur_zero = {ui: 0 for ui in ur}

   Mhd = fhd.jacobian(ur)
   ghd = fhd.xreplace(ur_zero)

   Mhd, ghd

.. jupyter-execute::

   ur_sol = -Mhd.LUsolve(ghd)
   ur_repl = dict(zip(ur, ur_sol))

.. jupyter-execute::

   fhdd = fhd.diff(t).xreplace(qd_repl).xreplace(ur_repl)
   me.find_dynamicsymbols(fhdd)

.. jupyter-execute::

   A.set_ang_vel(N, u1*N.z)
   B.set_ang_vel(A, ur_repl[u2]*A.z)
   C.set_ang_vel(B, ur_repl[u3]*B.z)

   P1.set_vel(N, 0)
   P2.v2pt_theory(P1, N, A)
   P3.v2pt_theory(P2, N, B)
   P4.v2pt_theory(P3, N, C)

   R_P2 = -m*g*N.y
   R_P3 = -m*g*N.y

.. jupyter-execute::

   Fr = sm.Matrix([
       P2.vel(N).diff(u1, N).dot(R_P2) +
       P3.vel(N).diff(u1, N).dot(R_P3)
   ])

   me.find_dynamicsymbols(Fr)

.. jupyter-execute::

   me.find_dynamicsymbols(P2.acc(N), reference_frame=N)

.. jupyter-execute::

   me.find_dynamicsymbols(P3.acc(N), reference_frame=N)

.. jupyter-execute::

   Rs_P2 = -m*P2.acc(N)
   Rs_P3 = -m*P3.acc(N).xreplace(qd_repl).xreplace(ur_repl)

   Frs = sm.Matrix([
       P2.vel(N).diff(u1, N).dot(Rs_P2) +
       P3.vel(N).diff(u1, N).dot(Rs_P3)
   ])
   me.find_dynamicsymbols(Frs)

.. jupyter-execute::

   ud = u.diff(t)
   ud_zero = {udi: 0 for udi in ud}

   Md = Frs.jacobian(ud)
   gd = Frs.xreplace(ud_zero) + Fr
   me.find_dynamicsymbols(Md), me.find_dynamicsymbols(gd)

.. jupyter-execute::

   gk = gk.xreplace(ur_repl)

   eval_k = sm.lambdify((qN, u, p), (Mk, gk))
   eval_d = sm.lambdify((qN, u, p), (Md, gd))

.. jupyter-execute::

   p_vals = np.array([
       0.8,  # la [m]
       2.0,  # lb [m]
       1.0,  # lc [m]
       2.0,  # ln [m]
       1.0,  # m [kg]
       9.81,  # g [m/s^2]
   ])

.. jupyter-execute::

   q1_val = np.deg2rad(10.0)

   eval_fh = sm.lambdify((qr, q, p), fh)

   from scipy.optimize import fsolve

   q2_val, q3_val = fsolve(lambda qr, q, p: np.squeeze(eval_fh(qr, [q], p)),
                           np.deg2rad([20.0, -150]),
                           args=(q1_val, p_vals))

   qN_vals = np.array([q1_val, q2_val, q3_val])
   np.rad2deg(qN_vals)


.. jupyter-execute::

   def eval_rhs(t, x, p):

       qN = x[:3]
       u = x[3:]

       Mk, gk = eval_k(qN, u, p)
       qNd = -np.linalg.solve(Mk, np.squeeze(gk))

       Md, gd = eval_d(qN, u, p)
       ud = -np.linalg.solve(Md, gd)[0]

       return np.hstack((qNd, ud))

.. jupyter-execute::

   u10 = 0.0
   x0 = np.hstack((qN_vals, u10))
   t0, tf = 0.0, 10.0
   fps = 30
   ts = np.linspace(t0, tf, num=int(fps*(tf - t0)))

   eval_rhs(t0, x0, p_vals)

.. jupyter-execute::

   sol = solve_ivp(eval_rhs, (t0, tf), x0, args=(p_vals,), t_eval=ts)

.. jupyter-execute::

   q1_traj, q2_traj, q3_traj, u1_traj = sol.y

   constraint_violations = []
   for i in range(len(sol.t)):
       constraint_violations.append(
           eval_fh((q2_traj[i], q3_traj[i]), [q1_traj[i]], p_vals)
       )

.. jupyter-execute::

   plt.plot(sol.t, sol.y.T)
   plt.legend(['q1', 'q2', 'q3', 'u1'])

.. jupyter-execute::

   plt.plot(sol.t, np.squeeze(constraint_violations))

.. jupyter-execute::

   coordinates = P2.pos_from(P1).to_matrix(N)
   for point in [P3, P4, P1, P2]:
      coordinates = coordinates.row_join(point.pos_from(P1).to_matrix(N))

   eval_point_coords = sm.lambdify((qN, p), coordinates)
   eval_point_coords(qN_vals, p_vals)

.. jupyter-execute::

   x, y, z = eval_point_coords(qN_vals, p_vals)

   fig, ax = plt.subplots()
   fig.set_size_inches((10.0, 10.0))
   ax.set_aspect('equal')
   ax.grid()

   lines, = ax.plot(x, y, color='black',
                    marker='o', markerfacecolor='blue', markersize=10)

   title_template = 'Time = {:1.2f} s'
   title_text = ax.set_title(title_template.format(t0))
   ax.set_xlim((-1.0, 3.0))
   ax.set_ylim((-1.0, 1.0))
   ax.set_xlabel('$x$ [m]')
   ax.set_ylabel('$y$ [m]');

.. jupyter-execute::

   xs = sol.y.T
   fps = 30

   coords = []
   for xi in xs:
        coords.append(eval_point_coords(xi[:3], p_vals))
   coords = np.array(coords)

   def animate(i):
       title_text.set_text(title_template.format(sol.t[i]))
       lines.set_data(coords[i, 0, :], coords[i, 1, :])

   ani = FuncAnimation(fig, animate, len(sol.t))

   HTML(ani.to_jshtml(fps=fps))

.. jupyter-execute::

   def eval_rhs(t, x, p):

       qN = x[:3]
       u = x[3:]

       # correct the depdendent coordinates
       qN[1:] = fsolve(lambda qr, q, p: np.squeeze(eval_fh(qr, [q], p)),
                       qN[1:],  # guess with current solution
                       args=(qN[0], p_vals))

       Mk, gk = eval_k(qN, u, p)
       qNd = -np.linalg.solve(Mk, np.squeeze(gk))

       Md, gd = eval_d(qN, u, p)
       ud = -np.linalg.solve(Md, gd)[0]

       return np.hstack((qNd, ud))

.. jupyter-execute::

   sol = solve_ivp(eval_rhs, (t0, tf), x0, args=(p_vals,), t_eval=ts)

   q1_traj, q2_traj, q3_traj, u1_traj = sol.y

   constraint_violations = []
   for i in range(len(sol.t)):
       constraint_violations.append(
           eval_fh((q2_traj[i], q3_traj[i]), [q1_traj[i]], p_vals)
       )

   plt.plot(sol.t, np.squeeze(constraint_violations))

.. jupyter-execute::

   x, y, z = eval_point_coords(qN_vals, p_vals)

   fig, ax = plt.subplots()
   fig.set_size_inches((10.0, 10.0))
   ax.set_aspect('equal')
   ax.grid()

   lines, = ax.plot(x, y, color='black',
                    marker='o', markerfacecolor='blue', markersize=10)

   title_template = 'Time = {:1.2f} s'
   title_text = ax.set_title(title_template.format(t0))
   ax.set_xlim((-1.0, 3.0))
   ax.set_ylim((-1.0, 1.0))
   ax.set_xlabel('$x$ [m]')
   ax.set_ylabel('$y$ [m]');

.. jupyter-execute::

   xs = sol.y.T
   fps = 30

   coords = []
   for xi in xs:
        coords.append(eval_point_coords(xi[:3], p_vals))
   coords = np.array(coords)

   def animate(i):
       title_text.set_text(title_template.format(sol.t[i]))
       lines.set_data(coords[i, 0, :], coords[i, 1, :])

   ani = FuncAnimation(fig, animate, len(sol.t))

   HTML(ani.to_jshtml(fps=fps))

https://github.com/bmcage/odes/blob/master/ipython_examples/Planar%20Pendulum%20as%20DAE.ipynb

Simulate using a DAE Solver
===========================

.. math::

   \bar{f}_k(\dot{\bar{q}}, \bar{u}, \bar{q}, \bar{q}_r, t)  = 0 \in \mathbb{R}^n \\
   \bar{f}_d(\dot{\bar{u}}, \bar{u}, \bar{q}, \bar{q}_r, t)  = 0 \in \mathbb{R}^n \\
   \bar{f}_h(\bar{q}, \bar{q}_r, t) = 0 \in \mathbb{R}^M

.. jupyter-execute::

   from scikits.odes import dae

.. jupyter-execute::

   def eval_eom(t, x, xd, residual, p):

       q1, q2, q3, u1 = x
       q1d, q2d, q3d, u1d = xd

       Md, gd = eval_d([q1, q2, q3], [u1], p)

       residual[0] = q1d - u1  # 1 equation
       residual[1] = Md[0]*u1d + gd[0]  # 1 equation
       residual[2:] = eval_fh([q2, q3], [q1], p).squeeze()  # 2 equation

   residual = np.empty(4)
   Md_vals, gd_vals = eval_d(qN_vals, [0.0], p_vals)
   xd0 = np.array([
      0.0,
      0.0,
      0.0,
      -np.linalg.solve(Md_vals, gd_vals)[0],
   ])
   eval_eom(t0, x0, xd0, residual, p_vals)
   residual

Options:

https://github.com/bmcage/odes/blob/1e3b3324748f4665ee5a52ed1a6e0b7e6c05be7d/scikits/odes/sundials/ida.pyx#L848

.. jupyter-execute::

   solver = dae('ida',
                lambda t, x, xd, res: eval_eom(t, x, xd, res, p_vals),
                #compute_initcond='yp0',
                first_step_size=1e-18,
                atol=1e-6,
                rtol=1e-6,
                algebraic_vars_idx=[2, 3],
                #compute_initcond_t0=60,
                old_api=False)
   solution = solver.solve(ts, x0, xd0)

.. jupyter-execute::

   ts_dae = solution.values.t
   xs_dae = solution.values.y

   plt.plot(ts_dae, xs_dae)
   plt.legend(['q1', 'q2', 'q3', 'u1'])

.. jupyter-execute::

   q1_traj, q2_traj, q3_traj, u1_traj = xs_dae.T

   constraint_violations = []
   for i in range(len(sol.t)):
       constraint_violations.append(
           eval_fh((q2_traj[i], q3_traj[i]), [q1_traj[i]], p_vals)
       )

   plt.plot(sol.t, np.squeeze(constraint_violations))

.. jupyter-execute::

   x, y, z = eval_point_coords(qN_vals, p_vals)

   fig, ax = plt.subplots()
   fig.set_size_inches((10.0, 10.0))
   ax.set_aspect('equal')
   ax.grid()

   lines, = ax.plot(x, y, color='black',
                    marker='o', markerfacecolor='blue', markersize=10)

   title_template = 'Time = {:1.2f} s'
   title_text = ax.set_title(title_template.format(t0))
   ax.set_xlim((-1.0, 3.0))
   ax.set_ylim((-1.0, 1.0))
   ax.set_xlabel('$x$ [m]')
   ax.set_ylabel('$y$ [m]');

.. jupyter-execute::

   coords = []
   for xi in xs_dae:
        coords.append(eval_point_coords(xi[:3], p_vals))
   coords = np.array(coords)

   def animate(i):
       title_text.set_text(title_template.format(ts_dae[i]))
       lines.set_data(coords[i, 0, :], coords[i, 1, :])

   ani = FuncAnimation(fig, animate, len(sol.t))

   HTML(ani.to_jshtml(fps=fps))
