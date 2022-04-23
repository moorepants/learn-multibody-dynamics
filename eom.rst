===========================================
Unconstrained Holonomic Equations of Motion
===========================================

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`eom` or Jupyter Notebook:
   :jupyter-download:notebook:`eom`.

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   me.init_vprinting(use_latex='mathjax')

In the previous chapter, we introduced the generalized active forces and the
generalized inertia forces. Together, these two pieces give us the *dynamical
differential equations*. The dynamical differential equations for a holonomic
system are defined as and they are a function of the generalized coordinates,
the generalized speeds, the time derivatives of the genaralized speeds, and
time:

.. math::
   :label: eq-kanes-equations

   \bar{F}_r + \bar{F}^*_r = \bar{f}_d(\dot{\bar{u}}, \bar{u}, \bar{q}, t)  = 0

We also call these equations *Kane's Equations* as they were introduced in this
form in [Kane1985]_. The dynamical differential equations can only be formed
with respect to an `inertial reference frame`_. An inertial reference frame is
one that is not accelerating, or can be assumed not to be with respect to the
motion of the bodies of interest. An inertial referene frame is one, where
Newton's First Law holds, i.e. objects at rest stay at rest unless an external
force acts on them.

.. _inertial reference frame: https://en.wikipedia.org/wiki/Inertial_frame_of_reference

:math:`\bar{F}^*_r` is always linear in the time derivatives of the generalized
speeds and contains velocity dependent terms such as the centripal and Coriolis
forces and the rotational velocity coupling terms. Dynamics texts will often
present it in this form, showing the linear nature:

.. math::
   :label: eq-canonical-eom-form
   \mathbf{M}(\bar{q}, t) \dot{\bar{u}} + \bar{C}(\bar{u}, \bar{q}, t) &= \bar{F}(\bar{u}, \bar{q}, t) \\
   -\bar{F}^*_r &= \bar{F}_r

where :math:`\mathbf{M}` is called the *mass matrix* and :math:`\bar{C}` is are
the forces due to the various velocity effects.

.. todo:: Same something about how M is always invertible and positive definite
   (I think).

The kinematical and dynamical differential equations constitute the *equations
of motion* for an unconstrained holonomic multibody system. These equations are
ordinary differential equations in the generalized speeds and generalized
coordinates.

.. math::
   :label: eq-equations-of-motion

   \bar{f}_d(\dot{\bar{u}}, \bar{u}, \bar{q}, t)  = 0 \\
   \bar{f}_k(\dot{\bar{q}}, \bar{u}, \bar{q}, t)  = 0

and since they are both linear in :math:`\dot{\bar{u}}` and
:math:`\dot{\bar{q}}`, respectively, they can be written in a combined form:

.. math::
   :label: eq-intermediate-state-form

   \begin{bmatrix}
   \mathbf{Y}_k && 0 \\
   0 && \mathbf{Y}_d \\
   \end{bmatrix}
   \begin{bmatrix}
   \dot{\bar{q}} \\
   \dot{\bar{u}}
   \end{bmatrix}
   +
   \begin{bmatrix}
   \bar{z}_k(\bar{u}, \bar{q}, t) \\
   \bar{z}_d(\bar{u}, \bar{q}, t)
   \end{bmatrix}
   =
   \begin{bmatrix}
   0 \\
   0
   \end{bmatrix}

which we write as:

.. math::
   :label: eq-state-form

   \mathbf{Y}
   \dot{\bar{x}}
   +
   \bar{z}
   = \bar{0}

where :math:`\bar{x}=[\bar{q} \quad \bar{u}]^T` is called the *state* of the
system and is comprised of the generalized coordinates and generalized speeds.

Returning to the example from the previous chapter, I will add a additional
particle of mass :math:`m/4` at point :math:`Q` that can slides along the rod
:math:`B` and is attached to point :math:`B_o` via a linear translational
spring with stiffness :math:`k_l` and located by generalized coordinate
:math:`q_3`. See :numref:`fig-eom-double-rod-pendulum` for a visual
description.

.. _fig-eom-double-rod-pendulum:
.. figure:: figures/eom-double-rod-pendulum.svg
   :align: center
   :width: 600px

   Three dimensional pendulum made up of two pinned rods and a sliding mass on
   rod :math:`B`. Each degree of freedom is resisted by a linear spring.

The following code is reproduced from the prior chapter and gives the
velocities and angular velocities of :math:`A_o`, :math:`B_o`, :math:`A`, and
:math:`B` in the inertial reference frame :math:`N`.

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

   Ao.vel(N), Bo.vel(N), A.ang_vel_in(N), B.ang_vel_in(N)

We now have the partical at :math:`Q` so we need its velocity for its
contribution to  :math:`F_r` and :math:`F_r^*`. :math:`Q` is moving in
:math:`B` so the one point velocity theorem can be used.

.. jupyter-execute::

   Q = me.Point('Q')
   Q.set_pos(Bo, q3*B.y)
   Q.set_vel(B, u3*B.y)
   Q.v1pt_theory(Bo, N, B)

   Q.vel(N)

We will also need the accelerations of the points and frames for the
generalized inertia forces. For points :math:`A_o`, :math:`B_o` and frames
:math:`A` and :math:`B` these are nicely expressed in terms of
:math:`\dot{\bar{u}}, \bar{u}, \bar{q}`:

.. jupyter-execute::

   Ao.acc(N), Bo.acc(N), A.ang_acc_in(N), B.ang_acc_in(N)

but the acceleration of point :math:`Q` contains :math:`\dot{\bar{q}}` terms,
so we need to eliminate those with the kinematical differential equations:

.. jupyter-execute::

   Q.acc(N)

.. jupyter-execute::

   t = me.dynamicsymbols._t

   qdot_repl = {q1.diff(t): u1,
                q2.diff(t): u2,
                q3.diff(t): u3}

   Q.set_acc(N, Q.acc(N).xreplace(qdot_repl))
   Q.acc(N)

Now we formulate the resultant forces and torques on each relevant point and
frame:

.. jupyter-execute::

   R_Ao = m*g*N.x
   R_Bo = m*g*N.x + kl*q3*B.y
   R_Q = m*g*N.x - kl*q3*B.y
   T_A = -kt*q1*N.z + kt*q2*A.x
   T_B = -kt*q2*A.x

The inertia dyadics of the two rods are:

.. jupyter-execute::

   I = m*l**2/12
   I_A_Ao = I*me.outer(A.y, A.y) + I*me.outer(A.z, A.z)
   I_B_Bo = I*me.outer(B.x, B.x) + I*me.outer(B.z, B.z)

With all of the necessary elements present for forming :math:`\bar{F}_r` and
:math:`\bar{F}_r^*` we can take advantage of Python for loops to systematically
formulate the generalized forces and inertia forces:

.. jupyter-execute::

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

The generalized forces are:

.. jupyter-execute::

   Fr = sm.Matrix(Fr)
   Fr

The generalized inertia forces are:

.. jupyter-execute::

   Frs = sm.Matrix(Frs)
   Frs

Forward Simulation
==================

Eq. :math:numref:`eq-state-form` is written in an *implicit form*, meaning that
the derivatives are not explicitly solved for. The *explicit form* is found by
inverting :math:`\mathbf{Y}`:

.. math::
   :label: eq-state-form

   \dot{\bar{x}}
   =
   -\mathbf{Y}^{-1}
   \bar{z}

To determine how the state changes over time, these explicit differential
equations can be solved:

.. math::
   :label: eq-eom-integral

   \bar{x}(t) = \int^{t_f}_{t_0} -\mathbf{Y}^{-1} \bar{z} dt = \int^{t_f}_{t_0} \bar{f}_m(\bar{u}, \bar{q}, t) dt

:math:`\bar{f}_m` is, in general, nonlinear in time, thus analytical solutions
are impossible to find. To solve this integral we must numerically integrate
:math:`\bar{f}_m`.

Our example problem has a simple definition of the kinematical differential
equations:

.. math::

   \begin{bmatrix}
   \dot{q}_1 \\
   \dot{q}_2 \\
   \dot{q}_3
   \end{bmatrix}
   =
   \begin{bmatrix}
   u_1 \\
   u_2 \\
   u_3
   \end{bmatrix}

so :math:`\mathbf{Y}_k` is the identity matrix and need not be formed. But we
will need :math:`\mathbf{Y}_d` to solve explicitly for :math:`\dot{\bar{u}`.

.. jupyter-execute::

   u = sm.Matrix([u1, u2, u3])
   ud = u.diff(t)

   Yd = Frs.jacobian(ud)
   Yd

.. jupyter-execute::

   ud_zerod = {udr: 0 for udr in ud}

   zd = Frs.xreplace(ud_zerod) + Fr
   zd
