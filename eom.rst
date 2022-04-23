===================
Equations of Motion
===================

In the previous chapter, we introduced the generalized active forces and the
generalized inertia forces. Together, these two pieces give us the dynamical
differential equations. The dynamical differential equations are:

.. math::
   :label: eq-kanes-equations

   \bar{F}_r + \bar{F}^*_r = \bar{f}_d(\dot{\bar{u}}, \bar{u}, \bar{q}, t)  = 0

We also call these equations *Kane's Equations* due to the formulation
presented in [Kane1985]_.

:math:`\bar{F}^*_r` is linear in the time derivatives of the generalized speeds
and also contains velocitiy dependent terms such as centripal and Coriolis
forces and rotational velocity couplings. Dynamics texts will often present
it in this form:

.. math::

   \bar{F}^*_r = \mathbf{M}(q, t) \dot{\bar{u}} + \bar{C}(\bar{u}, \bar{q}, t)

where :math:`\mathbf{M}` is called the *mass matrix* and :math:`\bar{C}` is are
the forces due to velocity effects.

The kinematical and dynamical differential equations constitute the *equations
of motion* for a holonomic multibody system. These equations are ordinary
differential equations in the generalized speeds and generalized coordinates.

.. math::
   :label: eq-equations-of-motion

   \bar{f}_d(\dot{\bar{u}}, \bar{u}, \bar{q}, t)  = 0 \\
   \bar{f}_k(\dot{\bar{q}}, \bar{u}, \bar{q}, t)  = 0

They are also linear in :math:`\dot{\bar{u}}` and :math:`\dot{\bar{q}}`,
respectively.

.. math::

   \begin{bmatrix}
   \bar{M} && 0 \\
   0 && \bar{Y}
   \end{bmatrix}
   \begin{bmatrix}
   \dot{\bar{u}} \\
   \dot{\bar{q}}
   \end{bmatrix}
   +
   \begin{bmatrix}
   \bar{f} \\
   \bar{f}
   \end{bmatrix}
   =
   \begin{bmatrix}
   0 \\
   0
   \end{bmatrix}

.. figure:: figures/eom-double-rod-pendulum.svg
   :align: center

   TODO: Add caption.

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   me.init_vprinting(use_latex='mathjax')

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
   Q = me.Point('Q')

   Ao.set_pos(O, l/2*A.x)
   Bo.set_pos(O, l*A.x)
   Q.set_pos(Bo, q3*B.y)

   O.set_vel(N, 0)
   Ao.v2pt_theory(O, N, A)
   Bo.v2pt_theory(O, N, A)
   Q.set_vel(B, u3*B.y)
   Q.v1pt_theory(Bo, N, B)

   R_Ao = m*g*N.x
   R_Bo = m*g*N.x + kl*q3*B.y
   R_Q = m*g*N.x - kl*q3*B.y
   T_A = -kt*q1*N.z + kt*q2*A.x
   T_B = -kt*q2*A.x

   I = m*l**2/12
   I_A_Ao = I*me.outer(A.y, A.y) + I*me.outer(A.z, A.z)
   I_B_Bo = I*me.outer(B.x, B.x) + I*me.outer(B.z, B.z)

.. jupyter-execute::

   points = [Ao, Bo, Q]
   forces = [R_Ao, R_Bo, R_Q]
   frames = [A, B]
   torques = [T_A, T_B]
   inertias = [I_A_Ao, I_B_Bo]

   t = me.dynamicsymbols._t

   qdot_repl = {q1.diff(t): u1,
                q2.diff(t): u2,
                q3.diff(t): u3}

   Fr = []
   Frs = []
   for ur in [u1, u2, u3]:

      Fri = 0
      Frsi = 0

      for Pi, Ri in zip(points, forces):
         vr = Pi.vel(N).diff(ur, N)
         Fri += vr.dot(Ri)
         Rs = -m*Pi.acc(N).xreplace(qdot_repl)
         Frsi += vr.dot(Rs)

      for Bi, Ti, Ii in zip(frames, torques, inertias):
         wr = Bi.ang_vel_in(N).diff(ur, N)
         Fri += wr.dot(Ti)
         Ts = -(Bi.ang_acc_in(N).dot(Ii) +
                me.cross(Bi.ang_vel_in(N), Ii).dot(
                Bi.ang_vel_in(N)).xreplace(qdot_repl))
         Frsi += wr.dot(Ts)

      Fr.append(Fri)
      Frs.append(Frsi)

.. jupyter-execute::

   Fr = sm.Matrix(Fr)
   Fr

.. jupyter-execute::

   Frs = sm.Matrix(Frs)
   Frs

.. jupyter-execute::

   q = sm.Matrix([q1, q2, q3])
   u = sm.Matrix([u1, u2, u3])

   t = me.dynamicsymbols._t
   ud = u.diff(t)

   ud_zerod = {udr: 0 for udr in ud}

.. jupyter-execute::

   M = Frs.jacobian(ud)
   M

.. jupyter-execute::

   F = Frs.xreplace(ud_zerod) + Fr
   F

.. jupyter-execute::

   Y = sm.eye(len(q))
   Y

.. jupyter-execute::

   Ts

Numerical Evaluation
====================

.. jupyter-execute::

   p = [m, g, kt, kl, l]

   eval_MF = sm.lambdify((q, u, p), [M, F])

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

   M_vals, F_vals = eval_MF(q_vals, u_vals, p_vals)
   M_vals, F_vals

.. jupyter-execute::


   ud_vals = np.linalg.solve(M_vals, F_vals)
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
       M, F = eval_MF(q, u, p)
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

   plt.plot(ts, xs, 'k', res.t, res.y.T, 'b')

How do we know that the equations of motion are correct?
========================================================
