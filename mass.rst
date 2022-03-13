====
Mass
====

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   me.init_vprinting(use_latex='mathjax')

In the prior chapters we have developed the tools to formulate the kinematics
of points and reference frames. The kinematics are the first of three essential
parts needed to form the equations of motion of a multibody system.

When a point is associated with a particle of mass :math:`m` or a reference
frame is associated with a body that has a mass distribution, Newton's and
Euler's Second Laws of motion show that the time rate of change of the linear
and angular momenta must be equal to the forces and torques acting on the
particle or body. The momentum of the particle is determined by its mass and
velocity and the angular momentum is determined by the distribution of mass and
the angular velocity.

Mass
====

Given a set of :math:`\nu` particles with masses :math:`m_1,\ldots,m_\nu` the
total mass, or *zeroth moment of mass*, is defined as:

.. math::
   :label: eq-zeroth-moment

   m := \sum_{i=1}^\nu m_i

For a body with a density :math:`\rho` defined at each point within its
volumetric boundary, the total mass becomes an integral of the general form:

.. math::
   :label: eq-zeroth-moment-rigid-body

   m := \int_{\textrm{volume}} \rho dV

Mass Center
===========

If each particle is located at positions :math:`\bar{r}_i,\ldots,\bar{r}_\nu`
relative to a point :math:`O` and the first mass moment, :math:`\sum_{i=1}^\nu
m_i \bar{r}_i` is equal to 0 then the point :math:`S_o` is refered to as the
mass center of the set of particles. The mass center is defined as:

.. math::
   :label: mass-center-particles

   \bar{r}^{S_o/O} = \frac{ \sum_{i=1}^\nu m_i \bar{r}_i }{\sum_{i=1}^\nu m_i}

.. math::
   :label: mass-center-rigid-body

   \bar{r}^{S_o/O} = \frac{ \int_{\textrm{mass}} \bar{r} dm }{ \int_{\textrm{mass}} dm }

.. jupyter-execute::

   m1, m2, m3 = sm.symbols("m1, m2, m3")
   x1, x2, x3 = me.dynamicsymbols('x1, x2, x3')
   y1, y2, y3 = me.dynamicsymbols('y1, y2, y3')
   z1, z2, z3 = me.dynamicsymbols('z1, z2, z3')

   A = me.ReferenceFrame('A')

   r_O_So = (m1*(x1*A.x + y1*A.y + z1*A.z) +
             m2*(x2*A.x + y2*A.y + z2*A.z) +
             m3*(x3*A.x + y3*A.y + z3*A.z)) / (m1 + m2 + m3)
   r_O_So

If :math:`m_2=2m_1` and :math:`m_3=3m_1` then:

.. jupyter-execute::

   r_O_So.xreplace({m2: 2*m1, m3: 3*m1}).simplify()

Mass Distribution
=================

The inertia, or second moment of mass, describes the distribution of mass
relative to a point.

For a set of particles :math:`P_1,\ldots,P_\nu` with locations
:math:`\bar{r}^{P_i/O}` for :math:`i=1,\ldots,\nu` relative to :math:`O` the
inertia vector for a unit vector :math:`\hat{n}_a` is defined as:

.. math::
   :label: inertia-vector-particles

   \bar{I}_a := \sum_{i=1}^\nu m_i \bar{r}^{P_i/O} \times \left( \hat{n}_a \times
   \bar{r}^{P_i/O}  \right)

.. todo:: Show the inertia vector for a rigid body.

For a single particle, the magnitude of :math:`\bar{I}_a` is:

.. math::
   :label: intertia-vector-magnitude

   \left| \bar{I}_a \right| = m \left| \bar{r}^{P/O} \right| ^2 \sin\theta

where :math:`\theta` is angle between :math:`\bar{r}^{P/O}` and
:math:`\hat{n}_a`. If :math:`\hat{n}_a` is perpendicular to
:math:`\bar{r}^{P/O}` then the magnitude is:

.. math::
   :label: intertia-vector-magnitude-perp

   \left| \bar{I}_a \right| = m \left| \bar{r}^{P/O} \right| ^2

The inertia vector for the set of :math:`\nu` particles relative to :math:`O`

A component of :math:`\bar{I}_a` in the :math:`\hat{n}_b` direction is called
an inertia scalar and is defined as:

.. math::
   :label: inertia-scalar

   I_{ab} := \hat{I}_{a} \cdot \hat{n}_b

The inertia scalar can be rewritten using Eq.
:math:numref:`inertia-vector-particles`:

.. math::

   I_{ab} =
   \sum_{i=1}^\nu m_i
   \left( \bar{r}^{P_i/O} \times \hat{n}_a \right)
   \cdot
   \left( \bar{r}^{P_i/O} \times \hat{n}_b \right)

This form implies that:

.. math::

   I_{ab} = I_{ba}

If :math:`\hat{n}_a = \hat{n}_b` then this inertia scalar is called a "moment
of inertia" and if :math:`\hat{n}_a \neq \hat{n}_b` it is called a "product of
inertia". When :math:`\hat{n}_a = \hat{n}_b`

.. math::

   I_{aa} =
   \sum_{i=1}^\nu m_i
   \left( \bar{r}^{P_i/O} \times \hat{n}_a \right)^2

The *radius of gyration* about a line through :math:`O` parallel to
:math:`\hat{n}_a` is defined as:

.. math::

   k_{aa} := \sqrt{\frac{I_{aa}}{m}}

Inertia Dyadic
==============

For mutually perpendicular unit vectors in :math:`A` the moments of inertia
with repsect to :math:`O` about each unit vector can be computed and the
products of inertia among the pairs of perpendiculr vectors can also be
computed. This, in general, generates six unique inertia scalars. These scalars
are typically presented as a symmetric inertia matrix (or inertia tensor):

.. math::

   \begin{bmatrix}
    I_{xx} & I_{xy} & I_{xz} \\
    I_{yx} & I_{yy} & I_{yz} \\
    I_{zx} & I_{zy} & I_{zz}
   \end{bmatrix}_A

where :math:`I_{xy}=I_{yx}`, :math:`I_{xz}=I_{zx}`, and :math:`I_{yz}=I_{zy}`
and the subscript :math:`A` indicates that these are relative to the unit
vectors of :math:`A`.

This second order tensor is similar to the the first order tensors (vectors)
we've already worked with:

.. math::

   \begin{bmatrix}
   a_1 \\
   0 \\
   0
   \end{bmatrix}_A

but we have notation for writing such a vector, i.e. :math:`a_1\hat{a}_x`.

If we introduce the outer product operator between two vectors:

.. jupyter-execute::

   a1, a2, a3 = sm.symbols('a1, a2, a3')
   b1, b2, b3 = sm.symbols('b1, b2, b3')

   a1*me.outer(A.x, A.x)

.. jupyter-execute::

   a1*me.outer(A.x, A.x).to_matrix(A)

.. jupyter-execute::

   v = a1*A.x + a2*A.y + a3*A.z

   B = me.ReferenceFrame('B')

   w = b1*B.x + b2*B.y + b3*B.z

.. jupyter-execute::

   theta = sm.symbols("theta")
   B.orient_axis(A, theta, A.x)

   Q = v.outer(w)
   Q

.. jupyter-execute::

   Q.to_matrix(A)

.. jupyter-execute::

   Q.to_matrix(B)

This is convenident because we can create dyadics, just like vectors, which are
make up of components in different reference frames:

.. jupyter-execute::

   P = 2*me.outer(B.x, B.x) + 3*me.outer(A.x, B.y) + 4*me.outer(B.z, A.z)
   P

.. jupyter-execute::

   P.express(A)

.. jupyter-execute::

   P.to_matrix(A)

.. todo:: Show dyadic properties including pre and post multiplication.

The unit dyadic is defined as:

.. jupyter-execute::

   U = me.outer(A.x, A.x) + me.outer(A.y, A.y) + me.outer(A.z, A.z)
   U

and it represents the identity matrix:

.. jupyter-execute::

   U.to_matrix(A)

Using the triple vector product identify the inertia vector can be written as:

.. math::

   \bar{I}_a = \sum_{i=1}^\nu m_i \left(\hat{n}_a \bar{r}^2 -
   \hat{n}_a\cdot\bar{r}\bar{r}\right)

.. math::

   \bar{I}_a =
   \sum_{i=1}^\nu m_i \left(
   \hat{n}_a \cdot \breve{U} \left(\bar{r}^{P_i/O}\right)^2 -
   \hat{n}_a \cdot \bar{r}^{P_i/O} \bar{r}^{P_i/O}
   \right)

We can now remove the dependence on :math:`\hat{n}_a` and arrive at the inertia
dyadic of the set of particles relative to :math:`O`:

.. math::

   \breve{I} :=
   \sum_{i=1}^\nu m_i \left(
   \breve{U} \left(\bar{r}^{P_i/O}\right)^2 -
   \bar{r}^{P_i/O} \bar{r}^{P_i/O}
   \right)


This form is basis independent and the position vectors can be express in any
coordinate system in any reference frame.

.. jupyter-execute::

   Ixx, Iyy, Izz, Ixy, Iyz, Ixz = sm.symbols('I_{xx}, I_{yy}, I_{zz}, I_{xy}, I_{yz}, I_{xz}')

   I = me.inertia(A, Ixx, Iyy, Izz, Ixy, Iyz, Ixz)
   I

.. jupyter-execute::

   I.to_matrix(A)

.. jupyter-execute::

   sm.trigsimp(I.to_matrix(B))

Angular Momentum
================



Parallel Axis
=============

Composite Inertias
==================

Principal Moments of Inertia
============================
