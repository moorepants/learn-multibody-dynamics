==========
Simulation
==========

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`simuluation` or Jupyter Notebook:
   :jupyter-download:notebook:`simulation`.

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   me.init_vprinting(use_latex='mathjax')

.. jupyter-execute::

   m, g, kt, kl, l = sm.symbols('m, g, k_t, k_l, l')
   q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')
   u1, u2, u3 = me.dynamicsymbols('u1, u2, u3')

   N = me.ReferenceFrame('N')
   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')

   A.orient_axis(N, q1, N.z)
   B.orient_axis(A, q2, A.x)

   A.set_ang_vel(N, u1*N.z)
   B.set_ang_vel(A, u2*A.x)

   O = me.Point('O')
   Ao = me.Point('A_O')
   Bo = me.Point('B_O')

   Ao.set_pos(O, l/2*A.x)
   Bo.set_pos(O, l*A.x)

   O.set_vel(N, 0)
   Ao.v2pt_theory(O, N, A)
   Bo.v2pt_theory(O, N, A)

   Q = me.Point('Q')
   Q.set_pos(Bo, q3*B.y)
   Q.set_vel(B, u3*B.y)
   Q.v1pt_theory(Bo, N, B)

   t = me.dynamicsymbols._t

   qdot_repl = {q1.diff(t): u1,
                q2.diff(t): u2,
                q3.diff(t): u3}

   Q.set_acc(N, Q.acc(N).xreplace(qdot_repl))

   R_Ao = m*g*N.x
   R_Bo = m*g*N.x + kl*q3*B.y
   R_Q = m*g*N.x - kl*q3*B.y
   T_A = -kt*q1*N.z + kt*q2*A.x
   T_B = -kt*q2*A.x

   I = m*l**2/12
   I_A_Ao = I*me.outer(A.y, A.y) + I*me.outer(A.z, A.z)
   I_B_Bo = I*me.outer(B.x, B.x) + I*me.outer(B.z, B.z)

   points = [Ao, Bo, Q]
   forces = [R_Ao, R_Bo, R_Q]
   masses = [m, m, m/4]

   frames = [A, B]
   torques = [T_A, T_B]
   inertias = [I_A_Ao, I_B_Bo]

   Fr = []
   Frs = []

   for ur in [u1, u2, u3]:

      Fri = 0
      Frsi = 0

      for Pi, Ri, mi in zip(points, forces, masses):
         vr = Pi.vel(N).diff(ur, N)
         Fri += vr.dot(Ri)
         Rs = -mi*Pi.acc(N)
         Frsi += vr.dot(Rs)

      for Bi, Ti, Ii in zip(frames, torques, inertias):
         wr = Bi.ang_vel_in(N).diff(ur, N)
         Fri += wr.dot(Ti)
         Ts = -(Bi.ang_acc_in(N).dot(Ii) +
                me.cross(Bi.ang_vel_in(N), Ii).dot(Bi.ang_vel_in(N)))
         Frsi += wr.dot(Ts)

      Fr.append(Fri)
      Frs.append(Frsi)

   Fr = sm.Matrix(Fr)
   Frs = sm.Matrix(Frs)

   u = sm.Matrix([u1, u2, u3])
   ud = u.diff(t)
   ud_zerod = {udr: 0 for udr in ud}

   Yd = Frs.jacobian(ud)
   zd = Frs.xreplace(ud_zerod) + Fr

Numerical Evaluation
====================

.. jupyter-execute::

   q = sm.Matrix([q1, q2, q3])

   p = [m, g, kt, kl, l]

   eval_Yd_zd = sm.lambdify((q, u, p), [Yd, zd])

.. jupyter-execute::

   import numpy as np

   q_vals = [
       np.deg2rad(5.0),  # rad
       np.deg2rad(5.0),  # rad
       0.1,  # m
   ]

   u_vals = [
       0.1,  # rad/s
       0.2,  # rad/s
       0.3,  # m/s
   ]

   p_vals = [
       1.0,  # kg
       9.81,  # m/s**2
       0.01,  # Nm/rad
       2.0,  # N/m
       0.6,  # m
   ]

   Yd_vals, zd_vals = eval_Yd_zd(q_vals, u_vals, p_vals)
   Yd_vals, zd_vals

.. jupyter-execute::

   ud_vals = np.linalg.solve(Yd_vals, zd_vals)
   ud_vals

Forward Simulation
==================

.. jupyter-execute::

   def eval_rhs(t, x, p):
       """Return the right hand side of the explicit ordinary differential
       equations.

       Parameters
       ==========
       t : float
          Time in seconds.
       x : array_like, shape(6,)
          State at time t: [q1, q2, q3, u1, u2, u3]
       p : array_like, shape(5,)
          Constant parameters: [m, g, kt, kl, l]

       Returns
       =======
       xd : ndarray, shape(6,)
           Derivative of the state with respect to time.

       """

       q = x[:3]
       u = x[3:]

       qd = u
       M, F = eval_Yd_zd(q, u, p)
       ud = np.linalg.solve(M, np.squeeze(F))

       xd = np.empty_like(x)
       xd[:3] = qd
       xd[3:] = ud

       return xd

   x0 = np.empty(6)
   x0[:3] = q_vals
   x0[3:] = u_vals

   eval_rhs(0.1, x0, p_vals)


.. math::

   \bar{x}_i =

.. jupyter-execute::

   def euler_integrate(rhs_func, tspan, initial_cond, p_vals):
       delt = 0.01  # seconds/sample
       num_samples = int((tspan[1] - tspan[0])/delt)
       ts = np.linspace(tspan[0], tspan[1], num=num_samples + 1)

       x = np.empty((len(ts), len(initial_cond)))

       # Set the initial conditions to the first element.
       x[0, :] = initial_cond

       # Use a for loop to sequentially calculate each new x.
       for i, ti in enumerate(ts[:-1]):
           x[i + 1, :] = x[i, :] + delt*rhs_func(ti, x[i, :], p_vals)

       return ts, x

.. jupyter-execute::

   ts, xs = euler_integrate(eval_rhs, (0.0, 2.0), x0, p_vals)

.. jupyter-execute::

   ts

.. jupyter-execute::

   type(ts), ts.shape

.. jupyter-execute::

   xs

.. jupyter-execute::

   type(xs), xs.shape

.. jupyter-execute::

   import matplotlib.pyplot as plt

   plt.plot(ts, xs);

.. jupyter-execute::

   def plot_results(ts, xs):

       fig, axes = plt.subplots(4, 1, sharex=True)

       fig.set_size_inches((10.0, 6.0))

       axes[0].plot(ts, np.rad2deg(xs[:, :2]))
       axes[1].plot(ts, xs[:, 2])
       axes[2].plot(ts, np.rad2deg(xs[:, 3:5]))
       axes[3].plot(ts, xs[:, 5])

       axes[0].legend([me.mlatex(q[0], mode='inline'),
                       me.mlatex(q[1], mode='inline')])
       axes[1].legend([me.mlatex(q[2], mode='inline')])
       axes[2].legend([me.mlatex(u[0], mode='inline'),
                       me.mlatex(u[1], mode='inline')])
       axes[3].legend([me.mlatex(u[2], mode='inline')])

       axes[0].set_ylabel('Angle [deg]')
       axes[1].set_ylabel('Distance [m]')
       axes[2].set_ylabel('Angular Rate [deg/s]')
       axes[3].set_ylabel('Speed [m/s]')

       axes[3].set_xlabel('Time [s]')

       fig.tight_layout()

       return axes

   plot_results(ts, xs)

.. jupyter-execute::

   from scipy.integrate import solve_ivp

   res = solve_ivp(eval_rhs, (0.0, 2.0), x0, args=(p_vals,))

.. jupyter-execute::

   plot_results(res.t, res.y.T)

.. jupyter-execute::

   plt.plot(ts, xs, 'k', res.t, res.y.T, 'b');

How do we know that the equations of motion are correct?
========================================================
