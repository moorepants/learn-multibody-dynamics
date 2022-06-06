=====================================================
Unconstrained Equations of Motion with the TMT Method
=====================================================

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`tmt` or Jupyter Notebook:
   :jupyter-download:notebook:`tmt`.

.. jupyter-execute::

   import numpy as np
   import sympy as sm
   import sympy.physics.mechanics as me
   me.init_vprinting(use_latex='mathjax')

There are many mathematical methods available to formulate the equations of
motion given a description of a multibody system model. The different methods
offer different advantages and disadvantages over the original methods proposed
by Netwon and Euler. For example, Joseph-Louis Lagrange developed a way to
arrive at the equations of motion from the descriptions of kinetic and
potential energy of the system. Sir William Hamilton then reformulated
Lagrange's approach in terms of generalized momenta. Since then, Gibbs, Appell,
Kane and others have formaulted other methods. In this chapter, we present an
alternative method presented in [Vallery2020]_ call the "TMT Method".

Vallery and Schwab present a matrix that can be used to directly transform the
mass and inertia of a system of rigid bodies expressed in a Cartesian
coordinate system into the mass matrix for a specific set of generalized
coordinates. This matrix is populated by the measure numbers of the partial
velocities expressed in the inertial reference frame. The method has some
advantages:

- Only the position and orientation of the bodies need be determined.
- Contributing forces can be formaluated as generalized forces, or not, or
  both.

Given :math:`\nu` rigid bodies in a multibody system, the velocities of each
mass center and the angular velocities of each body in an inertial reference
frame :math:`N` can be written in column vector :math:`\bar{v}` form by
extracting the measure numbers of each velocity term.

.. math::

   \bar{v}(\dot{\bar{q}}, \bar{q}, t) =
   \begin{bmatrix}
   {}^N\bar{v}^{B_{1o}} \cdot \hat{n}_x \\
   {}^N\bar{v}^{B_{1o}} \cdot \hat{n}_y \\
   {}^N\bar{v}^{B_{1o}} \cdot \hat{n}_z \\
   {}^N\bar{\omega}^{B_1} \cdot \hat{n}_x \\
   {}^N\bar{\omega}^{B_1} \cdot \hat{n}_y \\
   {}^N\bar{\omega}^{B_1} \cdot \hat{n}_z \\
   \vdots \\
   {}^N\bar{v}^{B_{\nu o}} \cdot \hat{n}_x \\
   {}^N\bar{v}^{B_{\nu o}} \cdot \hat{n}_y \\
   {}^N\bar{v}^{B_{\nu o}} \cdot \hat{n}_z \\
   {}^N\bar{\omega}^{B_\nu} \cdot \hat{n}_x \\
   {}^N\bar{\omega}^{B_\nu} \cdot \hat{n}_y \\
   {}^N\bar{\omega}^{B_\nu} \cdot \hat{n}_z \\
   \end{bmatrix}
   \in
   \mathbb{R}^{6\nu}

The measure numbers of the partial velocities with respect to each speed in
:math:`\dot{\bar{q}}` can be efficiently found with the Jacobian and stored in
a matrix :math:`\mathbf{T}`.

.. math::

   \mathbf{T} = J_\bar{v},\dot{\bar{q}}

The mass an inertia of each body can be expressed in a 

.. math::

   \mathbf{M}_{B_1} =
   \begin{bmatrix}
   m_{B_1} & 0 & 0 & 0 & 0 & 0 \\
   0 & m_{B_1} & 0 & 0 & 0 & 0 \\
   0 & 0 & m_{B_1} & 0 & 0 & 0 \\
   0 & 0 & 0 &
   \breve{I}^{B_1/B_{1o}} \cdot \hat{n}_x\hat{n}_x &
   \breve{I}^{B_1/B_{1o}} \cdot \hat{n}_x\hat{n}_y &
   \breve{I}^{B_1/B_{1o}} \cdot \hat{n}_x\hat{n}_z \\
   0 & 0 & 0 &
   \breve{I}^{B_1/B_{1o}} \cdot \hat{n}_y\hat{n}_x &
   \breve{I}^{B_1/B_{1o}} \cdot \hat{n}_y\hat{n}_y &
   \breve{I}^{B_1/B_{1o}} \cdot \hat{n}_y\hat{n}_z \\
   0 & 0 & 0 &
   \breve{I}^{B_1/B_{1o}} \cdot \hat{n}_z\hat{n}_x &
   \breve{I}^{B_1/B_{1o}} \cdot \hat{n}_z\hat{n}_y &
   \breve{I}^{B_1/B_{1o}} \cdot \hat{n}_z\hat{n}_z \\
   \end{bmatrix}

.. math::

   \mathbf{M} =
   \begin{bmatrix}
   \mathbf{M}_{B_1} & \mathbf{0}       & \ldots     & \mathbf{0} \\
   \mathbf{0}       & \mathbf{M}_{B_2} & \ldots     & \vdots \\
   \vdots           & \vdots           & \ddots     & \vdots \\
   \mathbf{0}       & \mathbf{0}       & \ldots     & \mathbf{M}_{B_\nu}
   \end{bmatrix}

Vallery and Schwab show that the the mass matrix of the system is:

.. math::

   \mathbf{M}_d = \mathbf{T}^T \mathbf{M} \mathbf{T}

Now given a vector :math:`\bar{F}` of resultant forces and torques of couples acting on each rigid body:

.. math::

   \bar{F} =
   \begin{bmatrix}
   \bar{R}^{B_{1o}} \cdot \hat{n}_x \\
   \bar{R}^{B_{1o}} \cdot \hat{n}_y \\
   \bar{R}^{B_{1o}} \cdot \hat{n}_z \\
   \bar{T}^{B_1} \cdot \hat{n}_x \\
   \bar{T}^{B_1} \cdot \hat{n}_y \\
   \bar{T}^{B_1} \cdot \hat{n}_z \\
   \vdots \\
   \bar{R}^{B_{2o}} \cdot \hat{n}_x \\
   \bar{R}^{B_{2o}} \cdot \hat{n}_y \\
   \bar{R}^{B_{2o}} \cdot \hat{n}_z \\
   \bar{T}^{B_2} \cdot \hat{n}_x \\
   \bar{T}^{B_2} \cdot \hat{n}_y \\
   \bar{T}^{B_2} \cdot \hat{n}_z \\
   \end{bmatrix}

.. math::

   \bar{g} = \frac{d\bar{v}}{dt}\bigg\rvert_{\ddot{\bar{q}}=\bar{0}}

.. math::

   \mathbf{T}^T \mathbf{M} \mathbf{T} \ddot{\bar{q}} =
   \mathbf{T}^T\left(\bar{F} - \mathbf{M}\bar{g}\right)

Example
=======

.. _fig-eom-double-rod-pendulum:
.. figure:: figures/eom-double-rod-pendulum.svg
   :align: center
   :width: 600px

   Three dimensional pendulum made up of two pinned rods and a sliding mass on
   rod :math:`B`. Each degree of freedom is resisted by a linear spring. When
   the generalized coordinates are all zero, the two rods are perpendicular to
   each other.

The following code is reproduced from the prior chapter and gives the
velocities and angular velocities of :math:`A_o`, :math:`B_o`, :math:`A`, and
:math:`B` in the inertial reference frame :math:`N`.

.. jupyter-execute::

   m, g, kt, kl, l = sm.symbols('m, g, k_t, k_l, l')
   q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')
   t = me.dynamicsymbols._t

   q = sm.Matrix([q1, q2, q3])
   q

.. jupyter-execute::

   N = me.ReferenceFrame('N')
   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')

   A.orient_axis(N, q1, N.z)
   B.orient_axis(A, q2, A.x)

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
   Q.set_vel(B, q3.diff()*B.y)
   Q.v1pt_theory(Bo, N, B)

   Ao.vel(N), A.ang_vel_in(N), Bo.vel(N), B.ang_vel_in(N), Q.vel(N)

.. jupyter-execute::

   R_Ao = m*g*N.x
   R_Bo = m*g*N.x + kl*q3*B.y
   R_Q = m/4*g*N.x - kl*q3*B.y
   T_A = -kt*q1*N.z + kt*q2*A.x
   T_B = -kt*q2*A.x

.. jupyter-execute::

   I = m*l**2/12
   I_A_Ao = I*me.outer(A.y, A.y) + I*me.outer(A.z, A.z)
   I_B_Bo = I*me.outer(B.x, B.x) + I*me.outer(B.z, B.z)

.. jupyter-execute::

   v = sm.Matrix([
       Ao.vel(N).dot(N.x),
       Ao.vel(N).dot(N.y),
       Ao.vel(N).dot(N.z),
       A.ang_vel_in(N).dot(N.x),
       A.ang_vel_in(N).dot(N.y),
       A.ang_vel_in(N).dot(N.z),
       Bo.vel(N).dot(N.x),
       Bo.vel(N).dot(N.y),
       Bo.vel(N).dot(N.z),
       B.ang_vel_in(N).dot(N.x),
       B.ang_vel_in(N).dot(N.y),
       B.ang_vel_in(N).dot(N.z),
       Q.vel(N).dot(N.x),
       Q.vel(N).dot(N.y),
       Q.vel(N).dot(N.z),
   ])
   v


.. jupyter-execute::

   MA = sm.diag(m, m, m).col_join(sm.zeros(3)).row_join(sm.zeros(3).col_join(I_A_Ao.to_matrix(N)))
   MA

.. jupyter-execute::

   MB = sm.diag(m, m, m).col_join(sm.zeros(3)).row_join(sm.zeros(3).col_join(I_B_Bo.to_matrix(N)))
   MB

.. jupyter-execute::

   MQ = sm.diag(m/4, m/4, m/4)
   MQ

.. jupyter-execute::

   M = sm.diag(MA, MB, MQ)
   sm.trigsimp(M)

.. jupyter-execute::

   F = sm.Matrix([
       R_Ao.dot(N.x),
       R_Ao.dot(N.y),
       R_Ao.dot(N.z),
       T_A.dot(N.x),
       T_A.dot(N.y),
       T_A.dot(N.z),
       R_Bo.dot(N.x),
       R_Bo.dot(N.y),
       R_Bo.dot(N.z),
       T_B.dot(N.x),
       T_B.dot(N.y),
       T_B.dot(N.z),
       R_Q.dot(N.x),
       R_Q.dot(N.y),
       R_Q.dot(N.z),
   ])
   F

.. jupyter-execute::

   T = v.jacobian(q.diff(t))
   T

.. jupyter-execute::

   qdd_repl = {qddi: 0 for qddi in q.diff(t, 2)}
   gz = v.diff(t).xreplace(qdd_repl)
   gz

.. jupyter-execute::

   Md = -sm.trigsimp(T.transpose()*M*T)
   Md

.. jupyter-execute::

   gd = sm.trigsimp(T.transpose()*(F - M*gz))
   gd

.. jupyter-execute::

   q_vals = np.array([
       np.deg2rad(25.0),  # q1, rad
       np.deg2rad(5.0),  # q2, rad
       0.1,  # q3, m
   ])

.. jupyter-execute::

   u_vals = np.array([
       0.1,  # u1, rad/s
       2.2,  # u2, rad/s
       0.3,  # u3, m/s
   ])

.. jupyter-execute::

   p_vals = np.array([
       9.81,  # g, m/s**2
       2.0,  # kl, N/m
       0.01,  # kt, Nm/rad
       0.6,  # l, m
       1.0,  # m, kg
   ])

.. todo:: gd is slightly differen than my prior solution

.. jupyter-execute::

   p = sm.Matrix([g, kl, kt, l, m])

   eval_d = sm.lambdify((q.diff(t), q, p), (Md, gd))

   eval_d(u_vals, q_vals, p_vals)

.. jupyter-execute::

   Md_vals, gd_vals = eval_d(u_vals, q_vals, p_vals)
   ud_vals = -np.linalg.solve(Md_vals, np.squeeze(gd_vals))
   ud_vals
