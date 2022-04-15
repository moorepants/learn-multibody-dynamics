==================
Generalized Forces
==================

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`generalized-forces` or Jupyter Notebook:
   :jupyter-download:notebook:`generalized-forces`.

At this point we have developed the three primary ingredients to formulate the
equations of motion of a multibody system:

1. Angular and Translational Kinematics
2. Mass and Mass Distribution
3. Forces, Moments, and Torques

For a single rigid body :math:`B` the `Newton-Euler Equations of Motion`_ can
be written as follows:

.. math::
   :label: eq-newton-euler

   \bar{F} = & \frac{{}^N d \bar{p}}{dt} \textrm{ where } \bar{p} = m_B{}^N\bar{v}^{B_o} \\
   \bar{T} = & \frac{{}^N d\bar{H}}{dt} \textrm{ where }
   \bar{H} = \breve{I}^{B/B_o} \cdot {}^N\bar{\omega}^{B}

.. _Newton-Euler Equations of Motion: https://en.wikipedia.org/wiki/Newton%E2%80%93Euler_equations

For a set of rigid bodies and particles that make up a multibody system,
defined with generalized coordinates and speeds, the generalized speeds
characterize completely the motion of the system. The velocities and angular
velocities of every particle and rigid body in the system are a function of
these generalized speeds. The time rate of change in the generalized speeds
will then play a criticale role of the right hand side of the equations of
motion, as shown in Eq. :math:numref:`eq-newton-euler`.

Take for example the multibody system shown in
:numref:`fig-generalized-forces-multi-pendulum`. A force :math:`\bar{F}`
applied at point :math:`Q` will cause all three of the lower particles to move.
The motion of the particles are described by the velocities, which are
functions of the generalized speeds. Thus :math:`\bar{F}` will, in general,
cause all of the generalized speeds to change. But how much does each
generalized speed change? The so called *partial velocites* of :math:`Q` in
:math:`N` will provide the answer to this question.

.. _fig-generalized-forces-multi-pendulum:
.. figure:: figures/generalized-forces-multi-pendulum.svg
   :align: center

   Four particles attached by massless links making up a 3 link planar simple
   pendulum.

Partial Velocities
==================

Recall that all translational and angular velocities of a multibody system can
be written in terms of the generalized speeds. These velocities can be
expressed uniquely as linear functions of the generalized speeds. For a
holonomic system with :math:`n` degrees of freedom any velocity or angular
velocity can be written as:

.. math::
   :label: eq-holonomic-partial-velocities

   \bar{v} = \sum_{r=1}^n \bar{v}_r u_r + \bar{v}_t \\
   \bar{\omega} = \sum_{r=1}^n \bar{\omega}_r u_r + \bar{\omega}_t

:math:`\bar{v}_r` and :math:`\bar{\omega}_r` are called the r\ :sup:`th`
holonomic partial velocity and angular velocity, respectively.
:math:`\bar{v}_t` and :math:`\bar{\omega}_t` are the remainders that are not
linear in a generalized speed. Since the velocities are linear in the
generalized speeds the partial velocities are equivalent to the partial
derivatives with respect to the generalized speeds:

.. math::
   :label: eq-partial-vel-partial-deriv

   \bar{v}_r = \frac{\partial v}{\partial u_r} \quad
   \bar{v}_t = \frac{\partial v}{\partial t} \\
   \bar{\omega}_r = \frac{\partial \omega}{\partial u_r} \quad
   \bar{\omega}_t = \frac{\partial \omega}{\partial t}


This means that the partial velocities can be considered being the sensitivity
of :math:`\bar{v}` and :math:`\bar{\omega}` to changes in :math:`u_r`.

Figure :numref:`fig-generalized-forces-partial-velocities` gives a graphical
interpretation of how a velocity of :math:`P` in :math:`N` is made up of
partial velocities and a remainder.

.. _fig-generalized-forces-partial-velocities:
.. figure:: figures/generalized-forces-partial-velocities.svg
   :align: center

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   me.init_vprinting(use_latex='mathjax')

.. todo:: Add figure for this example.

.. jupyter-execute::

   L = sm.symbols('L')

   q1, q2, u1, u2 = me.dynamicsymbols('q1, q2, u1, u2')

   N = me.ReferenceFrame('N')
   R = me.ReferenceFrame('R')

   R.orient_axis(N, q2, N.z)

   N_v_A = u1*N.x
   N_v_A

.. jupyter-execute::

   N_w_R = u2*N.z
   r_A_B = -L*R.x
   N_v_B = N_v_A + me.cross(N_w_R, r_A_B)

   N_v_B.express(N)

The partial velocities of :math:`P` in :math:`N` can be found by inspection or
using partial derivatives:

.. jupyter-execute::

   N_v_A.diff(u1, N), N_v_A.diff(u2, N)

.. jupyter-execute::

   N_v_B.diff(u1, N), N_v_B.diff(u2, N).express(N)

.. jupyter-execute::

   N_w_R.diff(u1, N), N_w_R.diff(u2, N)

When a system is nonholonomic, it is also true that every translational and
angular velocity can be expressed uniquely in terms of the :math:`p`
independent generalized speeds (see Eq
:math:numref:`eq-contraint-linear-form-solve`). Thus we can also defined the
*nonholonomic partial velocities* :math:`\tilde{v}_r` and *nonholonomic partial
angular velocities* :math:`\tilde{v}_r`.

.. math::
   :label: eq-holonomic-partial-velocities

   \bar{v} = \sum_{r=1}^p \tilde{v}_r u_r + \tilde{v}_t \\
   \bar{\omega} = \sum_{r=1}^p \tilde{\omega}_r u_r + \tilde{\omega}_t

.. todo:: Show how :math:`\bar{v}` and :math:`\tilde{v}` are related.

Generalized Active Forces
=========================

Suppose we have a holonomic multibody system make up of :math:`\nu` particles
with :math:`n` degrees of freedom in a reference frame :math:`A` that are
described by generalized speeds :math:`u_1,\ldots,u_n`. Each particle may have
a force :math:`\bar{R}` applied to it. We can define the *generalized active
forces*. The rth generalized active force for this system in A is defined as:

.. math::
   :label: eq-rth-gen-active-force

   F_r := \sum_{i=1}^\nu {}^A\bar{v}^{P_i} \cdot \bar{R}_i

This is the sum of the projections of the force acting on each particle in the
direction of the partial velocity. It projects the forces acting on the system
into contribution directions associatd with the rth generalized speed.

Notice that the r\ :sup:`th` geenralized active force is:

1. a scalar value
2. has contributions from all particles except if :math:`{}^N\bar{v}^{P_i}
   \perp \bar{R}_i`
3. associated with the r\ :sup:`th` generalized speed

We will typically collect all of the generalized active forces in a column
vector:

.. math::
   :label: eq-rth-gen-active-force

   \bar{F}_r = \begin{bmatrix}
   \sum_{i=1}^\nu {}^A\bar{v}_1^{P_i} \cdot \bar{R}_i \\
   \vdots \\
   \sum_{i=1}^\nu {}^A\bar{v}_r^{P_i} \cdot \bar{R}_i \\
   \vdots \\
   \sum_{i=1}^\nu {}^A\bar{v}_n^{P_i} \cdot \bar{R}_i
   \end{bmatrix}

.. todo:: Add point O to the figure.

.. _fig-generalized-forces-double-pendulum:
.. figure:: figures/generalized-forces-double-pendulum.svg
   :align: center

   Kinematic schematic of a simple double pendulum.

.. jupyter-execute::

   q1, q2, u1, u2 = me.dynamicsymbols('q1, q2, u1, u2')

   L = sm.symbols('L')

   N = me.ReferenceFrame('N')
   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')

   A.orient_axis(N, q1, N.z)
   B.orient_axis(N, q2, N.z)

   O = me.Point('O')
   P1 = me.Point('P1')
   P2 = me.Point('P2')

   O.set_vel(N, 0)

   P1.set_pos(O, -L*A.y)
   P2.set_pos(P1, -L*B.y)

   P1.v2pt_theory(O, N, A)
   P2.v2pt_theory(P1, N, B)

   P1.vel(N), P2.vel(N)

.. jupyter-execute::

   repl = {q1.diff(): u1, q2.diff(): u2}

   N_v_P1 = P1.vel(N).xreplace(repl)
   N_v_P2 = P2.vel(N).xreplace(repl)

   N_v_P1, N_v_P2

Now make free body diagram of each of the particles that shows all of the
forces acting on them, including both contributing and non-contributing forces.

.. jupyter-execute::

   T1, T2 = me.dynamicsymbols('T1, T2')
   m1, m2, g = sm.symbols('m1, m2, g')

   R1 = -m1*g*N.y + T1*A.y - T2*B.y
   R1

.. jupyter-execute::

   R2 = -m2*g*N.y + T2*B.y
   R2

.. jupyter-execute::

   v_P1_1 = N_v_P1.diff(u1, N)
   v_P1_2 = N_v_P1.diff(u2, N)
   v_P2_1 = N_v_P2.diff(u1, N)
   v_P2_2 = N_v_P2.diff(u2, N)
   v_P1_1, v_P1_2, v_P2_1, v_P2_2

.. jupyter-execute::

   F1 = me.dot(v_P1_1, R1) + me.dot(v_P2_1, R2)
   F1

.. jupyter-execute::

   F2 = me.dot(v_P1_2, R1) + me.dot(v_P2_2, R2)
   F2

Notice that the distance forces :math:`\bar{T}_1,\bar{T}_2` are not present in
the generalized active forces :math:`F_1` or :math:`F_2`. This is not by
coicendence, but will always be true for noncontributing forces. They are
infact named "nonconstributing" because they do not contribute to the
geenralized active forces.

If a holonmic multibody system with :math:`n` degrees of freedom includes a
rigid body :math:`B` then the loads acting on :math:`B` can be described by a
resultant force at an arbitrary point :math:`Q` of :math:`B` and a couple with
torque :math:`\bar{T}`.

.. math::
   :label: eq-generalized-active-force-rigid-body

   (F_r)_B = {}^A\bar{v}^Q_r \cdot \bar{R} + {}^A\bar{\omega}^B_r \cdot \bar{T}

:numref:`fig-generalized-forces-3d-rods` shows two thin rods of length
:math:`L` that are connected at points :math:`O` and :math:`B_o`.

.. _fig-generalized-forces-3d-rods:
.. figure:: figures/generalized-forces-3d-rods.svg
   :align: center

The first step is to define the necessary velocities we'll need: translational
velocities of the two mass centers and the angular velocities of each body. We
use the simple definition of the generalized speeds :math:`u_i=\dot{q}_i`.

.. jupyter-execute::

   m, g, k, L = sm.symbols('m, g, k, L')
   q1, q2, u1, u2 = me.dynamicsymbols('q1, q2, u1, u2')

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

   Ao.set_pos(O, L/2*A.x)
   Bo.set_pos(O, L*A.x)

   O.set_vel(N, 0)
   Ao.v2pt_theory(O, N, A)
   Bo.v2pt_theory(O, N, A)

   Ao.vel(N), Bo.vel(N), A.ang_vel_in(N), B.ang_vel_in(N)

Now determine the holonomic partial velocities:

.. jupyter-execute::

   v_Ao_1 = Ao.vel(N).diff(u1, N)
   v_Ao_2 = Ao.vel(N).diff(u2, N)
   v_Bo_1 = Bo.vel(N).diff(u1, N)
   v_Bo_2 = Bo.vel(N).diff(u2, N)

   v_Ao_1, v_Ao_2, v_Bo_1, v_Bo_2

and the holonomic partial angular velocities:

.. jupyter-execute::

   w_A_1 = A.ang_vel_in(N).diff(u1, N)
   w_A_2 = A.ang_vel_in(N).diff(u2, N)
   w_B_1 = B.ang_vel_in(N).diff(u1, N)
   w_B_2 = B.ang_vel_in(N).diff(u2, N)

   w_A_1, w_A_2, w_B_1, w_B_2

The resultant forces on the two bodies are simply the gravitional forces that
act at each mass center:

.. jupyter-execute::

   R_Ao = m*g*N.x
   R_Bo = m*g*N.x

   R_Ao, R_Bo

With linear torsion springs between A and N and A and B the torques acting on
each body are:

.. jupyter-execute::

   T_A = k*q1*N.z - k*q2*A.x
   T_B = k*q2*A.x

   T_A, T_B

A generalized active force component can be found for each body and each
generalized speed using :math:numref:`eq-generalized-active-force-rigid-body`:

.. jupyter-execute::

   F1_A = v_Ao_1.dot(R_Ao) + w_A_1.dot(T_A)
   F1_B = v_Bo_1.dot(R_Bo) + w_B_1.dot(T_B)
   F2_A = v_Ao_2.dot(R_Ao) + w_A_2.dot(T_A)
   F2_B = v_Bo_2.dot(R_Bo) + w_B_2.dot(T_B)

   F1_A, F1_B, F2_A, F2_B

Summing for each generalized speed and then stacking the two scalars in a
column vector gives the generalized active forces for the system:

.. jupyter-execute::

   F1 = F1_A + F1_B
   F2 = F2_A + F2_B

   Fr = sm.Matrix([F1, F2])
   Fr


For nonholonomic systems with :math:`p` degrees of freedom, the :math:`p`
generalized active forces can be formed instead:

.. math::
   :label: eq-nonholonomic-gaf

   (\tilde{F}_r)_P = {}^A\bar{v}^{P} \cdot \bar{R} \\
   (\tilde{F}_r)_B = {}^A\tilde{v}^Q \cdot \bar{R} + {}^A\tilde{\omega}^B \cdot \bar{T}

Generalized active forces are analagous to the left hand side of the
Newton-Euler equations but for a multibody system.

Generalized Inertia Forces
==========================

*Generalized inertia forces* map the right hand side of the Newton-Euler
equations to the generalized speeds for a multibody system.

For a multibody system made up of a set of particles the generalized inertia
forces are defined as

.. math::

   F_r^* := \sum_{i=1}^\nu {}^A\bar{v}^{P_i}_r \cdot \bar{R}^*_i

where

.. math::

   \bar{R}^*_i := -m_i {}^A\bar{a}^{P_i}_i

The generalized inertia force for a single rigid body :math:`B` with mass
center :math:`B_o` and central inerta dyadic :math:`\breve{I}^{B/Bo}` is
defined as:

.. math::

   (F_r^*)_B = {}^A\bar{v}^{B_o}_r \cdot \bar{R}^* + {}^A\bar{\omega}^B_r \cdot \bar{T}^*

.. math::

   \bar{R}^* := -m_i {}^A\bar{a}^{B_o}

.. math::

   \bar{T}^* := -\left(
   {}^A\bar{\alpha}^B \cdot \breve{I}^{B/Bo} +
   {}^A\bar{\omega}^B \times \breve{I}^{B/Bo} \cdot {}^A\bar{\omega}^B
   \right)


.. jupyter-execute::

   m, g, k, L = sm.symbols('m, g, k, L')
   q1, q2, u1, u2 = me.dynamicsymbols('q1, q2, u1, u2')

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

   Ao.set_pos(O, L/2*A.x)
   Bo.set_pos(O, L*A.x)

   O.set_vel(N, 0)
   Ao.v2pt_theory(O, N, A)
   Bo.v2pt_theory(O, N, A)

   v_Ao_1 = Ao.vel(N).diff(u1, N)
   v_Ao_2 = Ao.vel(N).diff(u2, N)
   v_Bo_1 = Bo.vel(N).diff(u1, N)
   v_Bo_2 = Bo.vel(N).diff(u2, N)

   w_A_1 = A.ang_vel_in(N).diff(u1, N)
   w_A_2 = A.ang_vel_in(N).diff(u2, N)
   w_B_1 = B.ang_vel_in(N).diff(u1, N)
   w_B_2 = B.ang_vel_in(N).diff(u2, N)

.. jupyter-execute::

   A.ang_acc_in(N), B.ang_acc_in(N)

.. jupyter-execute::

   Ao.acc(N), Bo.acc(N)

.. jupyter-execute::

   I = m*L**2/12

   I_A_Ao = I*me.outer(A.y, A.y) + I*me.outer(A.z, A.z)
   I_B_Bo = I*me.outer(B.x, B.x) + I*me.outer(B.z, B.z)

.. jupyter-execute::

   Rs_Ao = -m*Ao.acc(N)
   Rs_Bo = -m*Bo.acc(N)

   Rs_Ao, Rs_Bo

.. jupyter-execute::

   Ts_A = -(A.ang_acc_in(N).dot(I_A_Ao) + me.cross(A.ang_vel_in(N), I_A_Ao).dot(A.ang_vel_in(N)))
   Ts_B = -(B.ang_acc_in(N).dot(I_B_Bo) + me.cross(B.ang_vel_in(N), I_B_Bo).dot(B.ang_vel_in(N)))

   Ts_A, Ts_B

.. jupyter-execute::

   F1s_A = v_Ao_1.dot(Rs_Ao) + w_A_1.dot(Ts_A)
   F1s_B = v_Bo_1.dot(Rs_Bo) + w_B_1.dot(Ts_B)
   F2s_A = v_Ao_2.dot(Rs_Ao) + w_A_2.dot(Ts_A)
   F2s_B = v_Bo_2.dot(Rs_Bo) + w_B_2.dot(Ts_B)

   F1s = F1s_A + F1s_B
   F2s = F2s_A + F2s_B

   Frs = sm.Matrix([F1s, F2s])
   Frs

