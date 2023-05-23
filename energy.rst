================
Energy and Power
================

.. note::

   You can download this example as a Python script:
   :jupyter-download-script:`loads` or Jupyter Notebook:
   :jupyter-download-notebook:`loads`.

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   import numpy as np
   from scipy.optimize import fsolve
   me.init_vprinting(use_latex='mathjax')

.. container:: invisible

   .. jupyter-execute::

      class ReferenceFrame(me.ReferenceFrame):

          def __init__(self, *args, **kwargs):

              kwargs.pop('latexs', None)

              lab = args[0].lower()
              tex = r'\hat{{{}}}_{}'

              super(ReferenceFrame, self).__init__(*args,
                                                   latexs=(tex.format(lab, 'x'),
                                                           tex.format(lab, 'y'),
                                                           tex.format(lab, 'z')),
                                                   **kwargs)
      me.ReferenceFrame = ReferenceFrame

Learning Objectives
===================

After completing this chapter readers will be able to:

Introduction
============

So far we have invesitaged multibody systems from the perspect of forces and
their relationship to motion. It is also useful to understand these systems
from a power and energy perspective. Power is the time rate of change in work
done. Work is the eenrgy gained, disspated, or exchanged in a system.

.. power:: https://en.wikipedia.org/wiki/Power_(physics)

.. math::

   P = \frac{dW}{dt}

Knowing that work is a force dotted with a change in position, power can be
written as a force dotted with a velocity.

.. math::

   P = \bar{F} \cdot \bar{v}

Power can enter into a system, exit a system, or be exhanged within a system.

The time integral of power is work or energy. Energy of a multibody system can
be classified as kinetic, potential (conservative), or non-conservative.

.. math::

   E_k = m \bar{v} \cdot \bar{v} / 2  + \bar{\omega} \cdot \breve{I} \cdot \bar{\omega} / 2

.. math::

   E_p = 

Jumping
=======

.. jupyter-execute::

   g = sm.symbols('g')
   mu, mt, mc, mf = sm.symbols('m_u, m_t, m_c, m_f')
   It, Ic = sm.symbols('I_t, I_c')
   kc, cc, kk, ck = sm.symbols('k_c, c_c, k_k, c_k')
   lt, lc, dt, dc = sm.symbols('l_t, l_c, d_t, d_c')

   q1, q2, q3 = me.dynamicsymbols('q1, q2, q3', real=True)
   u1, u2, u3 = me.dynamicsymbols('u1, u2, u3', real=True)
   Tk = me.dynamicsymbols('T_k')

   t = me.dynamicsymbols._t

.. jupyter-execute::

   N = me.ReferenceFrame('N')
   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')

   A.orient_axis(N, q2, N.z)
   B.orient_axis(A, q3, N.z)

   A.set_ang_vel(N, u2*N.z)
   B.set_ang_vel(A, u3*N.z)

   O = me.Point('O')
   Pu, Pk, Pf = me.Point('P_u'), me.Point('P_k'), me.Point('P_f')
   Ao, Bo = me.Point('A_o'), me.Point('B_o')

   Pf.set_pos(O, q1*N.x)
   Ao.set_pos(Pf, dc*A.y)
   Pk.set_pos(Pf, lc*A.y)
   Bo.set_pos(Pk, dt*B.y)
   Pu.set_pos(Pk, lt*B.y)

   O.set_vel(N, 0)
   Pf.set_vel(N, u1*N.x)
   Pk.v2pt_theory(Pf, N, A)
   Pu.v2pt_theory(Pk, N, B)

   qd_repl = {q1.diff(t): u1, q2.diff(t): u2, q3.diff(t): u3}
   qdd_repl = {q1.diff(t, 2): u1.diff(t), q2.diff(t, 2): u2.diff(t), q3.diff(t, 2): u3.diff(t)}

   holonomic = Pu.pos_from(O).dot(N.y)
   vel_con = holonomic.diff(t).xreplace(qd_repl)
   acc_con = vel_con.diff(t).xreplace(qdd_repl).xreplace(qd_repl)

   # q2 is dependent

   u2_repl = {u2: sm.solve(vel_con, u2)[0]}
   u2d_repl = {u2.diff(t): sm.solve(acc_con, u2.diff(t))[0].xreplace(u2_repl)}

   R_Pu = -mu*g*N.x
   R_Ao = -mt*g*N.x
   R_Bo = -mc*g*N.x
   zp = (sm.Abs(q1) - q1)/2
   zd = zp.diff(t).xreplace(qd_repl)
   Fc = (kc*zp**(sm.S(3)/2) + cc*zp**(sm.S(3)/2)*zd)*N.x
   R_Pf = -mf*g*N.x + Fc

   T_A = (kk*q3 + ck*u3 + Tk)*N.z
   T_B = -T_A

   I_A_Ao = It*me.outer(N.z, N.z)
   I_B_Bo = Ic*me.outer(N.z, N.z)

   points = [Pu, Ao, Bo, Pf]
   forces = [R_Pu, R_Ao, R_Bo, R_Pf]
   masses = [mu, mt, mc, 0]

   frames = [A, B]
   torques = [T_A, T_B]
   inertias = [I_A_Ao, I_B_Bo]

   Fr_bar = []
   Frs_bar = []

   for ur in [u1, u3]:

      Fr = 0
      Frs = 0

      for Pi, Ri, mi in zip(points, forces, masses):
         N_v_Pi = Pi.vel(N).xreplace(u2_repl)
         vr = N_v_Pi.diff(ur, N)
         Fr += vr.dot(Ri)
         N_a_Pi = Pi.acc(N).xreplace(u2d_repl).xreplace(u2_repl)
         Rs = -mi*N_a_Pi
         Frs += vr.dot(Rs)

      for Bi, Ti, Ii in zip(frames, torques, inertias):
         N_w_Bi = Bi.ang_vel_in(N).xreplace(u2_repl)
         N_alp_Bi = Bi.ang_acc_in(N).xreplace(u2d_repl).xreplace(u2_repl)
         wr = N_w_Bi.diff(ur, N)
         Fr += wr.dot(Ti)
         Ts = -(N_alp_Bi.dot(Ii) + me.cross(N_w_Bi, Ii).dot(N_w_Bi))
         Frs += wr.dot(Ts)

      Fr_bar.append(Fr)
      Frs_bar.append(Frs)

   Fr = sm.Matrix(Fr_bar)
   Frs = sm.Matrix(Frs_bar)
   kane_eq = Fr + Frs

   q = sm.Matrix([q1, q2, q3])
   u = sm.Matrix([u1, u3])
   ud = u.diff(t)
   p = sm.Matrix([Ic, It, cc, ck, dc, dt, g, kc, kk, lc, lt, mc, mf, mt, mu])
   r = sm.Matrix([Tk])

.. todo:: cse fails

.. jupyter-execute::

   eval_kane = sm.lambdify((q, ud, u, r, p), kane_eq) #, cse=True)
   eval_holo = sm.lambdify((q, p), holonomic) #, cse=True)
   eval_vel_con = sm.lambdify((q, u, p), vel_con) #, cse=True)


   def eval_eom(t, x, xd, residual, p):
       """Returns the residual vector of the equations of motion.

       Parameters
       ==========
       t : float
          Time at evaluation.
       x : ndarray, shape(5,)
          State vector at time t: x = [q1, q2, q3, u1, u3].
       xd : ndarray, shape(5,)
          Time derivative of the state vector at time t: xd = [q1d, q2d, q3d, u1d, u3d].
       residual : ndarray, shape(5,)
          Vector to store the residuals in: residuals = [fk, fd, fh1, fh2].
       p : ndarray, shape(6,)
          Constant parameters: p = [la, lb, lc, ln, m, g]

       """

       q1, q2, q3, u1, u3 = x
       q1d, _, q3d, u1d, u3d = xd  # ignore the q2d value

       residual[0] = -q1d + u1
       residual[1] = -q3d + u3
       residual[2:4] = eval_kane([q1, q2, q3], [u1d, u3d], [u1, u3], [0], p).squeeze()
       residual[4] = eval_holo([q1, q2, q3], p)

.. jupyter-execute::

   residual = np.empty(5)
   eval_eom(1.0, np.random.random(5), np.random.random(5), residual, np.random.random(15))
   residual

.. jupyter-execute::

   p_vals = np.array([
     0.101,  # Ic,
     0.282,  # It,
     0.0,  # cc,
     0.0,  # ck,
     0.387,  # dc,
     0.193,  # dt,
     9.81,  # g,
     1000.0,  # kc,
     10.0,  # kk,
     0.611,  # lc,
     0.424,  # lt,
     6.769,  # mc,
     3.0,  # mf,
     17.01,  # mt,
     32.44,  # mu
   ])

   q0 = np.array([
       0.1,
       np.nan,
       np.deg2rad(60.0),
   ])

   q0[1] = fsolve(lambda q2: eval_holo([q0[0], q2, q0[2]], p_vals), np.deg2rad(-5.0))

   u0 = np.array([
       0.0,
       np.nan,
       0.0,
   ])

   u0[1] = fsolve(lambda u2: eval_vel_con(q0, [u0[0], u2, u0[2]], p_vals),  np.deg2rad(-5.0))

   x0 = np.hstack((q0, u0))

   ud0 = np.array([
   0.0,
   0.0,
   0.0])

   xd0 = np.hstack((u0, ud0))


.. jupyter-execute::

   solver = dae('ida',
                eval_eom,
                rtol=1e-3,
                atol=1e-6,
                algebraic_vars_idx=[4],
                user_data=p_vals,
                old_api=False)

.. jupyter-execute::

   ts = np.linspace(0.0, 5.0, num=101)

   solution = solver.solve(ts, x0, xd0)

   ts_dae = solution.values.t
   xs_dae = solution.values.y
