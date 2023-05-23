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
   kc, cc = sm.symbols('k_c, c_c')
   kk, ck = sm.symbols('k_k, c_k')
   lt, lc, dt, dc = sm.symbols('lt, lc, dt, dc')

   q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')
   u1, u2, u3 = me.dynamicsymbols('u1, u2, u3')
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

   O, Pu, Pk, Pf = me.Point('O'), me.Point('P_u'), me.Point('P_k'), me.Point('P_f')
   Ao, Bo = me.Point('A_o'), me.Point('B_o')

   Pu.set_pos(O, q1*N.x)
   Ao.set_pos(Pu, dt*A.y)
   Pk.set_pos(Pu, lt*A.y)
   Bo.set_pos(Pk, dc*B.y)
   Pf.set_pos(Pk, lc*B.y)

   O.set_vel(N, 0)
   Pu.set_vel(N, u1*N.x)
   Pk.v2pt_theory(Pu, N, A)
   Pf.v2pt_theory(Pk, N, B)

   qd_repl = {q1.diff(t): u1, q2.diff(t): u2, q3.diff(t): u3}
   qdd_repl = {q1.diff(t, 2): u1.diff(t), q2.diff(t, 2): u2.diff(t), q3.diff(t, 2): u3.diff(t)}

   holonomic = Pf.pos_from(O).dot(N.y)
   vel_con = holonomic.diff(t).xreplace(qd_repl)
   acc_con = vel_con.diff(t).xreplace(qdd_repl).xreplace(qd_repl)

   # q2 is dependent

   u2_repl = {u2: sm.solve(vel_con, u2)[0]}
   u2d_repl = {u2.diff(t): sm.solve(acc_con, u2.diff(t))[0].xreplace(u2_repl)}

   R_Pu = -mu*g*N.x
   R_Ao = -mt*g*N.x
   R_Bo = -mc*g*N.x
   zh = Pf.pos_from(O).dot(N.x)
   zp = (sm.Abs(zh) - zh)/2
   Fc = (kc*zp**(sm.S(3)/2) + cc*zp**(sm.S(3)/2)*zp.diff(t))*N.x
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
         vr = Pi.vel(N).xreplace(u2_repl).diff(ur, N)
         Fr += vr.dot(Ri)
         Rs = -mi*Pi.acc(N).xreplace(u2d_repl).xreplace(u2_repl)
         Frs += vr.dot(Rs)

      for Bi, Ti, Ii in zip(frames, torques, inertias):
         N_w_Bi = Bi.ang_vel_in(N).xreplace(u2_repl)
         wr = N_w_Bi.diff(ur, N)
         Fr += wr.dot(Ti)
         Ts = -(Bi.ang_acc_in(N).xreplace(u2d_repl).xreplace(u2_repl).dot(Ii) +
                  me.cross(N_w_Bi, Ii).dot(N_w_Bi))
         Frs += wr.dot(Ts)

      Fr_bar.append(Fr)
      Frs_bar.append(Frs)

   Fr = sm.Matrix(Fr_bar)
   Frs = sm.Matrix(Frs_bar)

   q = sm.Matrix([q1, q2, q3])
   u = sm.Matrix([u1, u3])
   p = sm.Matrix([g])
