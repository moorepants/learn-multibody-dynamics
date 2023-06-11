=====================================================
Unconstrained Equations of Motion with the TMT Method
=====================================================

.. note::

   You can download this example as a Python script:
   :jupyter-download-script:`tmt` or Jupyter Notebook:
   :jupyter-download-notebook:`tmt`.

.. jupyter-execute::

   import numpy as np
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

There are `several mathematical methods`_ available to formulate the equations
of motion of a multibody system. These different methods offer various
advantages and disadvantages over Newton and Euler's original formulations and
among each other. For example, `Joseph-Louis Lagrange`_ developed a way to
arrive at the equations of motion from the descriptions of kinetic and
potential energy of the system. `Sir William Hamilton`_ then reformulated
Lagrange's approach in terms of *generalized momenta* instead of energy. Since
then, `Gibbs & Appell`_, Kane [Kane1985]_, and others have proposed more
methods. In this chapter, we present one of these alternative methods called
the "TMT Method".  The details and derivation of the TMT Method can be found in
[Vallery2020]_ .

.. _several mathematical methods: https://en.wikipedia.org/wiki/Classical_mechanics
.. _Joseph-Louis Lagrange: https://en.wikipedia.org/wiki/Lagrangian_mechanics
.. _Sir William Hamilton: https://en.wikipedia.org/wiki/Hamiltonian_mechanics
.. _Gibbs & Appell: https://en.wikipedia.org/wiki/Appell%27s_equation_of_motion

Vallery and Schwab show how the collection of Newton-Euler equations for each
individual rigid body can be transformed into the reduced dynamical
differential equations associated with the generalized coordinates, speeds, and
accelerations using the :math:`\mathbf{T}` matrix. This :math:`\mathbf{T}`
matrix is populated by the measure numbers of the partial velocities expressed
in the inertial reference frame.

Given :math:`\nu` rigid bodies in a multibody system described by :math:`n`
generalized coordinates and generalized speeds, the velocities of each mass
center and the angular velocities of each body in an inertial reference frame
:math:`N` can be written in column vector :math:`\bar{v}` form by extracting
the measure numbers in the inertial reference frame :math:`N` of each velocity
term.

.. math::

   \bar{v}(\bar{u}, \bar{q}, t) =
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

The measure numbers of the partial velocities with respect to each of the
:math:`n` generalized speeds in :math:`\bar{u}` can be efficiently found by
taking the Jacobian of :math:`\bar{v}` with respect to the generalized speeds,
which we will name matrix :math:`\mathbf{T}`.

.. math::

   \mathbf{T} = \mathbf{J}_{\bar{v},\bar{u}} \in \mathbb{R}^{6\nu \times n}
   \quad
   \textrm{where}
   \quad
   \bar{v} = \mathbf{T} \bar{u}

For each set of six rows in :math:`\bar{v}` tied to a single rigid body, the
associated mass and inertia for that body can be written as:

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

Multiplying the velocities with this matrix gives the momenta of each rigid
body.

.. math::

   \mathbf{M}_{B_1} \bar{v}_{B_1} =
   \begin{bmatrix}
   \bar{p}^{B_{1o}} \cdot \hat{n}_x \\
   \bar{p}^{B_{1o}} \cdot \hat{n}_y \\
   \bar{p}^{B_{1o}} \cdot \hat{n}_z \\
   \bar{H}^{B_1/B_{1o}} \cdot \hat{n}_x \\
   \bar{H}^{B_1/B_{1o}} \cdot \hat{n}_y \\
   \bar{H}^{B_1/B_{1o}} \cdot \hat{n}_z \\
   \end{bmatrix}

The matrices for each rigid body can then be assembled into a matrix for the
entire set of rigid bodies.

.. math::

   \mathbf{M} =
   \begin{bmatrix}
   \mathbf{M}_{B_1} & \mathbf{0}       & \ldots     & \mathbf{0} \\
   \mathbf{0}       & \mathbf{M}_{B_2} & \ldots     & \vdots \\
   \vdots           & \vdots           & \ddots     & \vdots \\
   \mathbf{0}       & \mathbf{0}       & \ldots     & \mathbf{M}_{B_\nu}
   \end{bmatrix}

Allowing the momenta of all the rigid bodies to be found by matrix
multiplication of :math:`\mathbf{M} \bar{v}`.

A vector :math:`\bar{F}` of resultant forces and torques of couples acting on
each rigid body can be formed in a similar manner as :math:`\bar{v}`, by
extracting the measure numbers in the inertial reference frame.

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

The dynamical differential equations for the entire Newton-Euler system are
then:

.. math::

   \frac{d \mathbf{M} \bar{v}}{dt} = \bar{F} \in \mathbb{R}^{6\nu}

We know that selecting :math:`n` generalized coordinates for such a system
allows us to write the dynamical differential equations as a set of :math:`n`
equations which is, in general, much smaller than :math:`6\nu` equations due to
the large number of holonomic constraints that represent the connections of all
the bodies in the system. Vallery and Schwab show that the mass matrix
:math:`\mathbf{M}_d` for this reduced set of equations can be efficiently
calculated using the :math:`\mathbf{T}` matrix ([Vallery2020]_, pg. 349):

.. math::

   \mathbf{M}_d = -\mathbf{T}^T \mathbf{M} \mathbf{T}

and that the forces not proportional to the generalized accelerations is found
with:

.. math::

   \bar{g}_d = \mathbf{T}^T\left(\bar{F} - \bar{g}\right)

where [#]_:

.. math::

   \bar{g} = \frac{d\mathbf{M}\bar{v}}{dt}\bigg\rvert_{\dot{\bar{u}}=\bar{0}}

.. [#] Note that my :math:`\bar{g}` is slightly different than the one
   presented in [Vallery2020]_ to make sure the time derivative of the angular
   momenta are properly calculated.

The equations of motion then take this form:

.. math::

   \bar{0} =
   \mathbf{M}_d\dot{\bar{u}} + \bar{g}_d =
   -\mathbf{T}^T \mathbf{M} \mathbf{T} \dot{\bar{u}} +
   \mathbf{T}^T\left(\bar{F} - \bar{g}\right)

These equations are equivalent to Kane's Equations.

Example Formulation
===================

Let us return once again to the holonomic system introduced in :ref:`Example of
Kane's Equations`.

.. _fig-eom-double-rod-pendulum-again:
.. figure:: figures/eom-double-rod-pendulum.svg
   :align: center
   :width: 600px

   Three dimensional pendulum made up of two pinned rods and a sliding mass on
   rod :math:`B`. Each degree of freedom is resisted by a linear spring. When
   the generalized coordinates are all zero, the two rods are perpendicular to
   each other.

Start by introducing the variables.

.. jupyter-execute::

   m, g, kt, kl, l = sm.symbols('m, g, k_t, k_l, l')
   q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')
   u1, u2, u3 = me.dynamicsymbols('u1, u2, u3')
   t = me.dynamicsymbols._t

   q = sm.Matrix([q1, q2, q3])
   u = sm.Matrix([u1, u2, u3])
   p = sm.Matrix([g, kl, kt, l, m])
   q, u, p

The derivation of the kinematics is done in the same way as before.

.. jupyter-execute::

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

   Ao.vel(N), A.ang_vel_in(N), Bo.vel(N), B.ang_vel_in(N), Q.vel(N)

Only the contributing forces need be declared (noncontributing would cancel out
in the TMT transformation if included). Do not forget Newton's Third Law and be
sure to include the equal and opposite reactions.

.. jupyter-execute::

   R_Ao = m*g*N.x
   R_Bo = m*g*N.x + kl*q3*B.y
   R_Q = m/4*g*N.x - kl*q3*B.y
   T_A = -kt*q1*N.z + kt*q2*A.x
   T_B = -kt*q2*A.x

The inertia dyadics of each body will be needed.

.. jupyter-execute::

   I = m*l**2/12
   I_A_Ao = I*me.outer(A.y, A.y) + I*me.outer(A.z, A.z)
   I_B_Bo = I*me.outer(B.x, B.x) + I*me.outer(B.z, B.z)

Create the TMT Components
=========================

The vector :math:`\bar{v}` is formed from the velocities and angular velocities
of each rigid body or particle.

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

The inertial matrices for each body and the particle :math:`Q` are:

.. jupyter-execute::

   MA = sm.diag(m, m, m).col_join(sm.zeros(3)).row_join(sm.zeros(3).col_join(I_A_Ao.to_matrix(N)))
   MA

.. jupyter-execute::

   MB = sm.diag(m, m, m).col_join(sm.zeros(3)).row_join(sm.zeros(3).col_join(I_B_Bo.to_matrix(N)))
   sm.trigsimp(MB)

.. jupyter-execute::

   MQ = sm.diag(m/4, m/4, m/4)
   MQ

Note that these matrices change with time because we've expressed the inertia
scalars in the inertial reference frame :math:`N`. The matrices for all of the
bodies can be assembled into :math:`\mathbf{M}`:

.. jupyter-execute::

   M = sm.diag(MA, MB, MQ)

:math:`\bar{F}` is constructed to match the order of :math:`\bar{v}`:

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

These are the components we need to form the reduced dynamical differential
equations.

Formulate the reduced equations of motion
=========================================

First find :math:`\mathbf{T}` using the Jacobian:

.. jupyter-execute::

   T = v.jacobian(u)
   T

and then compute :math:`\bar{g}`:

.. jupyter-execute::

   qd_repl = dict(zip(q.diff(t), u))
   ud_repl = {udi: 0 for udi in u.diff(t)}
   gbar = (M*v).diff(t).xreplace(qd_repl).xreplace(ud_repl)
   sm.trigsimp(gbar)

The reduced mass matrix is then formed with
:math:`-\mathbf{T}^T\mathbf{M}\mathbf{T}`:

.. jupyter-execute::

   Md = sm.trigsimp(-T.transpose()*M*T)
   Md

and the reduced remainder term is formed with :math:`\mathbf{T}^T(\bar{F} -
\bar{g})`:

.. jupyter-execute::

   gd = sm.trigsimp(T.transpose()*(F - gbar))
   gd

Evaluate the equations of motion
================================

Now we can check to see if these dynamical differential equations are the same
as the ones we found with Kane's Method by evaluating them with the same set of
numbers we used in :ref:`Numerical Evaluation`. The input values were:

.. jupyter-execute::

   u_vals = np.array([
       0.1,  # u1, rad/s
       2.2,  # u2, rad/s
       0.3,  # u3, m/s
   ])

   q_vals = np.array([
       np.deg2rad(25.0),  # q1, rad
       np.deg2rad(5.0),  # q2, rad
       0.1,  # q3, m
   ])

   p_vals = np.array([
       9.81,  # g, m/s**2
       2.0,  # kl, N/m
       0.01,  # kt, Nm/rad
       0.6,  # l, m
       1.0,  # m, kg
   ])

We can lambdify ``Md`` and ``gq`` to see if these give the same values as those
found with Kane's Equations:

.. jupyter-execute::

   eval_d = sm.lambdify((u, q, p), (Md, gd))

   Md_vals, gd_vals = eval_d(u_vals, q_vals, p_vals)
   Md_vals, gd_vals

These numerical arrays are identical to our prior results. The state
derivatives then should also be identical:

.. jupyter-execute::

   eval_d(u_vals, q_vals, p_vals)
   ud_vals = -np.linalg.solve(Md_vals, np.squeeze(gd_vals))
   ud_vals

which they are. We can be fairly confident that Kane's method and the TMT
method result in the same equations of motion for this system.
