==================
Generalized Forces
==================

At this point we have the three primary ingredients to formulate the equations
of motion of a multibody system:

1. Angular and Translational Kinematics
2. Mass and Mass Distribution
3. Forces, Moments, and Torques

For a single rigid body with constant mass and inertia, its equations of motion
can be written as:

.. math::
   :label: eq-newton-euler

   \bar{F} = \frac{{}^N d\bar{p}}{dt} \\
   \bar{M} = \frac{{}^N d\bar{H}}{dt} \\
   \bar{p} = m\bar{v} \\
   \bar{H} = \breve{I} \cdot \bar{\omega}

For a set of rigid bodies and particles that make up a multibody system,
defined with generalized coordinates and speeds, the change in the generalized
speeds. The generalized speeds characterize completely the motion of the
system.

Take for example the multibody system shown in
:numref:`fig-generalized-forces-multi-pendulum`. A force :math:`\bar{F}`
applied at point :math:`Q` will cause all of the particles to move. The motion
of the particles are described by the velocities, which are functions of the
generalized speeds. Thus :math:`\bar{F}` will cause all of the generalized
speeds to change. But how much does each :math:`u_i` change? The so called
*partial velocites* of :math:`Q` in :math:`N` provide the answer to this
question.

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
holonomic system with :math:`n` degrees of freedom the velocities can be
written as so:

.. math::
   :label: eq-holonomic-partial-velocities

   \bar{v} = \sum_{r=1}^n \bar{v}_r u_r + \bar{v}_t \\
   \bar{\omega} = \sum_{r=1}^n \bar{\omega}_r u_r + \bar{\omega}_t

:math:`\bar{v}_r` and :math:`\bar{\omega}_r` are called the rth partial
velocity and angular velocity of :math:`B` in :math:`N`. :math:`\bar{v}_t` and
:math:`\bar{\omega}_t` are the remainders. Since the velocities are linear in
the generalized speeds:

.. math::

   \bar{\omega}_r = \frac{\partial \omega}{\partial u_r} \\
   \bar{v}_r = \frac{\partial v}{\partial u_r} \\

.. math::

   \bar{\omega}_t = \frac{\partial \omega}{\partial t} \\
   \bar{v}_t = \frac{\partial v}{\partial t} \\

Thus the partial velocities can be considered as the sensitivity of
:math:`\bar{v}` and :math:`\bar{\omega}` to changes in :math:`u_r`.

Figure :numref:`fig-generalized-forces-partial-velocities` gives a graphical
interpretation of how a velocity of :math:`P` in :math:`N` is made up of
partial velocities and a remainder.

.. _fig-generalized-forces-partial-velocities:
.. figure:: figures/generalized-forces-partial-velocities.svg
   :align: center

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   #me.init_vprinting(use_latex='mathjax')

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

Generalized Active Forces
=========================

Generalized Inertia Forces
==========================
