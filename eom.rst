=================================
Unconstrained Equations of Motion
=================================

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`eom` or Jupyter Notebook:
   :jupyter-download:notebook:`eom`.

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   me.init_vprinting(use_latex='mathjax')

Dynamical Differential Equations
================================

In the previous chapter, we introduced the generalized active forces and the
generalized inertia forces. Together, these two pieces give us the *dynamical
differential equations*. The dynamical differential equations for a holonomic
system are defined as and they are a function of the generalized coordinates,
the generalized speeds, the time derivatives of the generalized speeds, and
time:

.. math::
   :label: eq-kanes-equations

   \bar{F}_r + \bar{F}^*_r = \bar{f}_d(\dot{\bar{u}}, \bar{u}, \bar{q}, t)  = 0

These are the `Newton-Euler equations`_ for a multibody system in the form
presented in [Kane1985]_, thus we also call these equations *Kane's Equations*.
The dynamical differential equations can only be formed with respect to an
`inertial reference frame`_. An inertial reference frame is one that is not
accelerating, or can be assumed not to be with respect to the motion of the
bodies of interest. An inertial reference frame is one, where Newton's First
Law holds, i.e. objects at rest stay at rest unless an external force acts on
them.

.. _Newton-Euler equations: https://en.wikipedia.org/wiki/Newton%E2%80%93Euler_equations
.. _inertial reference frame: https://en.wikipedia.org/wiki/Inertial_frame_of_reference

:math:`\bar{F}^*_r` is always linear in the time derivatives of the generalized
speeds and contains velocity dependent terms such as the centripetal and Coriolis
forces and the rotational velocity coupling terms. Dynamics texts will often
present it in this form, showing the linear nature:

.. math::
   :label: eq-canonical-eom-form

   -\bar{F}^*_r &= \bar{F}_r \rightarrow
   \mathbf{M}(\bar{q}, t) \dot{\bar{u}} + \bar{C}(\bar{u}, \bar{q}, t) &= \bar{F}(\bar{u}, \bar{q}, t) \\

where :math:`\mathbf{M}` is called the *mass matrix*,  :math:`\bar{C}` is are
the forces due to the various velocity effects (also call the fictitious
forces), and :math:`\bar{F}` are the externally applied forces.

.. todo:: Same something about how M is always invertible and positive definite
   (I think).

Body Fixed Newton-Euler Equations
==================================

To show that Kane's Equations are equivalent to the Newton-Euler equations you
may have seen before, we can find the dynamical differential equations for a
single rigid body using Kane's method and then show the results in the
canonical form. For a rigid body :math:`B` moving in an inertial reference
frame :math:`A` with its velocity and angular velocity expressed in body fixed
coordinates and acted upon by a resultant force at the mass center and a moment
about the mass center.

.. jupyter-execute::

   m, Ixx, Iyy, Izz = sm.symbols('m, I_{xx}, I_{yy}, I_{zz}')
   Ixy, Iyz, Ixz = sm.symbols('I_{xy}, I_{yz}, I_{xz}')
   Fx, Fy, Fz, Mx, My, Mz = me.dynamicsymbols('F_x, F_y, F_z, M_x, M_y, M_z')
   u1, u2, u3, u4, u5, u6 = me.dynamicsymbols('u1, u2, u3, u4, u5, u6')

   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')

   Bo = me.Point('Bo')

Define the angular velocity of the body and the velocity of the mass center,
but express them simply in the body fixed coordinates.

.. jupyter-execute::

   A_w_B = u4*B.x + u5*B.y + u6*B.z
   B.set_ang_vel(A, A_w_B)

   A_v_Bo = u1*B.x + u2*B.y + u3*B.z
   Bo.set_vel(A, A_v_Bo)

Now find the six partial velocities and partial angular velocities. Note that
we use the ``var_in_dcm=False`` keyword argument. We do this because the
generalized speeds are not present in the unspecified direction cosine matrix
relating :math:`A` and :math:`B`. This allows the derivative in :math:`A` to be
formed without use of a direction cosine matrix. Generalized speeds will never
be present in a direction cosine matrix.

.. jupyter-execute::

   v_Bo_1 = A_v_Bo.diff(u1, A, var_in_dcm=False)
   v_Bo_2 = A_v_Bo.diff(u2, A, var_in_dcm=False)
   v_Bo_3 = A_v_Bo.diff(u3, A, var_in_dcm=False)
   v_Bo_4 = A_v_Bo.diff(u4, A, var_in_dcm=False)
   v_Bo_5 = A_v_Bo.diff(u5, A, var_in_dcm=False)
   v_Bo_6 = A_v_Bo.diff(u6, A, var_in_dcm=False)

   v_Bo_1, v_Bo_2, v_Bo_3, v_Bo_4, v_Bo_5, v_Bo_6

.. jupyter-execute::

   w_B_1 = A_w_B.diff(u1, A, var_in_dcm=False)
   w_B_2 = A_w_B.diff(u2, A, var_in_dcm=False)
   w_B_3 = A_w_B.diff(u3, A, var_in_dcm=False)
   w_B_4 = A_w_B.diff(u4, A, var_in_dcm=False)
   w_B_5 = A_w_B.diff(u5, A, var_in_dcm=False)
   w_B_6 = A_w_B.diff(u6, A, var_in_dcm=False)

   w_B_1, w_B_2, w_B_3, w_B_4, w_B_5, w_B_6

The ``partial_velocity()`` function does this same thing. Notice that due to
our velocity definitions, we get a very simple set of partial velocities.

.. jupyter-execute::

   par_vels = me.partial_velocity([A_v_Bo, A_w_B], [u1, u2, u3, u4, u5, u6], A)

   par_vels

Now form the generalized active forces:

.. jupyter-execute::

   T_B = Mx*B.x + My*B.y + Mz*B.z
   R_Bo = Fx*B.x + Fy*B.y + Fz*B.z

   F1 = v_Bo_1.dot(R_Bo) + w_B_1.dot(T_B)
   F2 = v_Bo_2.dot(R_Bo) + w_B_2.dot(T_B)
   F3 = v_Bo_3.dot(R_Bo) + w_B_3.dot(T_B)
   F4 = v_Bo_4.dot(R_Bo) + w_B_4.dot(T_B)
   F5 = v_Bo_5.dot(R_Bo) + w_B_5.dot(T_B)
   F6 = v_Bo_6.dot(R_Bo) + w_B_6.dot(T_B)

   Fr = sm.Matrix([F1, F2, F3, F4, F4, F6])
   Fr

and the generalized inertia forces:

.. jupyter-execute::

   I = me.inertia(B, Ixx, Iyy, Izz, Ixy, Iyz, Ixz)

   Rs = -m*Bo.acc(A)
   Ts = -(B.ang_acc_in(A).dot(I) + me.cross(A_w_B, I).dot(A_w_B))

   F1s = v_Bo_1.dot(Rs) + w_B_1.dot(Ts)
   F2s = v_Bo_2.dot(Rs) + w_B_2.dot(Ts)
   F3s = v_Bo_3.dot(Rs) + w_B_3.dot(Ts)
   F4s = v_Bo_4.dot(Rs) + w_B_4.dot(Ts)
   F5s = v_Bo_5.dot(Rs) + w_B_5.dot(Ts)
   F6s = v_Bo_6.dot(Rs) + w_B_6.dot(Ts)

   Frs = sm.Matrix([F1s, F2s, F3s, F4s, F5s, F6s])
   Frs

and finally Kane's Equations:

.. jupyter-execute::

   Fr + Frs

We can put this in canonical form (Eq. :math:numref:`eq-canonical-eom-form`)by
extracting the mass matrix, which is the linear coefficient matrix of
:math:`\bar{u}`:

.. jupyter-execute::

   u = sm.Matrix([u1, u2, u3, u4, u5, u6])
   t = me.dynamicsymbols._t
   ud = u.diff(t)

The mass matrix is:

.. jupyter-execute::

   M = -Frs.jacobian(ud)
   M

The fictitious forces vector is:

.. jupyter-execute::

   C = -Frs.xreplace({udi: 0 for udi in ud})
   C

And the forcing vector is:

.. jupyter-execute::

   F = Fr
   F

This example may seem overly complicated when using Kane's method, but it is a
systematic method that works for any number of rigid bodies and particles in a
system.

Equations of Motion
===================

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

Example of Kane's Equations
===========================

Returning to the example from the previous chapter, I will add an additional
particle of mass :math:`m/4` at point :math:`Q` that can slide along the rod
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

We now have the particle at :math:`Q` so we need its velocity for its
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
   R_Q = m/4*g*N.x - kl*q3*B.y
   T_A = -kt*q1*N.z + kt*q2*A.x
   T_B = -kt*q2*A.x

Note the equal and opposite spring forces that act on the pairs of points and
pairs of reference frames. We ignored the reaction torque on :math:`N` from
:math:`A` because :math:`N` is our inertial reference frame.

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

Notice that the dynamical differential equations are only functions of the time
varying variables :math:`\dot{\bar{u}},\bar{u},\bar{q}`:

.. jupyter-execute::

   me.find_dynamicsymbols(Fr)

.. jupyter-execute::

   me.find_dynamicsymbols(Frs)

Implicit and Explicit Form
==========================

Eq. :math:numref:`eq-state-form` is written in an *implicit form*, meaning that
the derivatives are not explicitly solved for. The *explicit form* is found by
inverting :math:`\mathbf{Y}`:

.. math::
   :label: eq-state-form-explicit

   \dot{\bar{x}}
   =
   -\mathbf{Y}^{-1}
   \bar{z}
   =\bar{f}_m(\bar{x}, t)

To determine how the state changes over time, these explicit differential
equations can be solved by integrating them with respect to time:

.. math::
   :label: eq-eom-integral

   \bar{x}(t) = \int^{t_f}_{t_0} \bar{f}_m(\bar{x}, t) dt

:math:`\bar{f}_m` is, in general, nonlinear in time, thus analytical solutions
are impossible to find. To solve this integral we must numerically integrate
:math:`\bar{f}_m`. To do so, it will be useful to extract the symbolic forms of
:math:`\mathbf{Y}_k`, :math:`\bar{z}_k`, :math:`\mathbf{Y}_d`, and
:math:`\bar{z}_d`.

Our example problem has a simple definition of the kinematical differential
equations:

.. math::
   :label: eq-qdot-equals-u

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

so :math:`\mathbf{Y}_k` is the identity matrix and need not be formed:

.. math::
   :label: eq-yk-identity

   \mathbf{Y}_k \dot{\bar{q}} = \bar{u}
   \rightarrow
   \begin{bmatrix}
   1 & 0 & 0 \\
   0 & 1 & 0 \\
   0 & 0 & 1 \\
   \end{bmatrix}
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

But we will need :math:`\mathbf{Y}_d` to solve explicitly for
:math:`\dot{\bar{u}}`. Recall that we can use the Jacobian to extract the
linear coefficients of :math:`\dot{\bar{u}}` and then find the terms that
aren't functions of :math:`\dot{\bar{u}}` by substitution (See Sec.
:ref:`sec-solving-linear-systems`).

Form the column vector :math:`\dot{\bar{u}}`:

.. jupyter-execute::

   u = sm.Matrix([u1, u2, u3])
   ud = u.diff(t)
   ud

Extract the coefficients of :math:`\dot{\bar{u}}`:

.. jupyter-execute::

   Yd = Frs.jacobian(ud)
   Yd

Make a substitution dictionary to set :math:`\dot{\bar{u}}=\bar{0}`:

.. jupyter-execute::

   ud_zerod = {udr: 0 for udr in ud}
   ud_zerod

Find :math:`\bar{z}_d` with :math:`\bar{z}_d =
\bar{F}_r^* |_{\dot{\bar{u}}=\bar{0}} + \bar{F}_r`:

.. jupyter-execute::

   zd = Frs.xreplace(ud_zerod) + Fr
   zd

Check that neither are functions of :math:`\dot{\bar{u}}`:

.. jupyter-execute::

   me.find_dynamicsymbols(Yd)

.. jupyter-execute::

   me.find_dynamicsymbols(zd)
