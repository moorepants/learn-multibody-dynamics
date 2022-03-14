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

   m := \int_{\textrm{solid}} \rho dV

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

   \bar{r}^{S_o/O} = \frac{ \int_{\textrm{solid}} \bar{r} dm }{ \int_{\textrm{solid}} dm }

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

Inertia Matrix
==============

For mutually perpendicular unit vectors fixed in reference frame :math:`A` the
moments of inertia with repsect to :math:`O` for each unit vector can be
computed and the products of inertia among the pairs of perpendiculr unit
vectors can also be computed. This, in general, generates nine unique inertia
scalars. These scalars are typically presented as a symmetric inertia matrix
(or inertia tensor):

.. math::

   \begin{bmatrix}
    I_{xx} & I_{xy} & I_{xz} \\
    I_{yx} & I_{yy} & I_{yz} \\
    I_{zx} & I_{zy} & I_{zz}
   \end{bmatrix}_A

where :math:`I_{xy}=I_{yx}`, :math:`I_{xz}=I_{zx}`, and :math:`I_{yz}=I_{zy}`
and the subscript :math:`A` indicates that these scalars are relative to unit
vectors :math:`\hat{a}_x,\hat{a}_y,\hat{a}_z`.

This second order tensor is similar to the the first order tensors (vectors)
we've already worked with:

.. math::

   \begin{bmatrix}
   v_1 \\
   v_2 \\
   v_3
   \end{bmatrix}_A

but we have notation for writing such a vector that allows us to combine
components expressed in different reference frames:

.. math::

   v_1\hat{a}_x + v_2\hat{a}_y + v_3\hat{a}_z

There also exists an analgous form for second order tensors.

Dyadics
=======

If we introduce the `outer product`_ operator between two vectors we see that
it generates a matrix akin to the inertia matrix above.

.. math::

   \begin{bmatrix}
   v_1 \\ v_2 \\ v_3
   \end{bmatrix}_A
   \otimes
   \begin{bmatrix}
     w_1 \\ w_2 \\ w_3
   \end{bmatrix}_A
   =
   \begin{bmatrix}
   v_1w_1 & v_1w_2 & v_1w_3 \\
   v_2w_1 & v_2w_2 & v_2w_3 \\
   v_3w_1 & v_3w_2 & v_3w_3 \\
   \end{bmatrix}_A

.. _outer product: https://en.wikipedia.org/wiki/Outer_product

In SymPy Mechanics outer products can be taken between two vectors:

.. jupyter-execute::

   v1, v2, v3 = sm.symbols('v1, v2, v3')
   w1, w2, w3 = sm.symbols('w1, w2, w3')

   A = me.ReferenceFrame('A')

   v = v1*A.x + v2*A.y + v3*A.z
   w = w1*A.x + w2*A.y + w3*A.z

   Q = me.outer(v, w)
   Q

but the result is not the matrix form show above, but instead the result is a
dyadic_. The dyadic is the analogous form for second order tensors as that
we've been using for first order tensors. The matrix form can be found with
:external:py:meth:`~sympy.physics.vector.dyadic.Dyadic.to_matrix`:

.. _dyadic: https://en.wikipedia.org/wiki/Dyadics

.. jupyter-execute::

   Q.to_matrix(A)

The dyadic is made up of scalars multiplied by unit dyads. Examples of unit
dyads are:

.. jupyter-execute::

   me.outer(A.x, A.x)

.. jupyter-execute::

   me.outer(A.x, A.x).to_matrix(A)

.. jupyter-execute::

   me.outer(A.y, A.z)

.. jupyter-execute::

   me.outer(A.y, A.z).to_matrix(A)

These unit dyads can be formed from unit vectors that are fixed in different
reference frames. This is convenident because we can create dyadics, just like
vectors, which are make up of components in different reference frames:

.. jupyter-execute::

   theta = sm.symbols("theta")
   B = me.ReferenceFrame('B')
   B.orient_axis(A, theta, A.x)

.. jupyter-execute::

   P = 2*me.outer(B.x, B.x) + 3*me.outer(A.x, B.y) + 4*me.outer(B.z, A.z)
   P

.. jupyter-execute::

   P.express(A)

.. jupyter-execute::

   P.to_matrix(A)

The unit dyadic is defined as:

.. jupyter-execute::

   U = me.outer(A.x, A.x) + me.outer(A.y, A.y) + me.outer(A.z, A.z)
   U

and it represents the identity matrix:

.. jupyter-execute::

   U.to_matrix(A)

.. todo:: ReferenceFrame should have an attribute that returns the unit dyadic
   (or dyads).

Properties of dyadics
=====================

- Scalar multiplication: :math:`\alpha(\bar{u}\otimes\bar{v}) = \alpha\bar{u}\otimes\bar{v} = \bar{u}\otimes\alpha\bar{v}`
- Distributive: :math:`\bar{u}\otimes(\bar{v} + \bar{w}) = \bar{u}\otimes\bar{v} + \bar{u}\otimes\bar{w}`
- Left and right dot product with a vector (results in a vector):

  - :math:`\bar{u}\cdot(\bar{v}\otimes\bar{w}) = (\bar{u}\cdot\bar{v})\bar{w}`
  - :math:`(\bar{u}\otimes\bar{v})\cdot\bar{w} = \bar{u}(\bar{v}\cdot\bar{w})`

- Left and right cross product with a vector (results in a dyadic):

  - :math:`\bar{u}\times(\bar{v}\otimes\bar{w}) = (\bar{u}\times\bar{v})\otimes\bar{w}`
  - :math:`(\bar{u}\otimes\bar{v})\times\bar{w} = \bar{u}\otimes(\bar{v}\times\bar{w})`

- Not commutative: :math:`\breve{U}\bar{u} \neq \bar{u}\breve{U}`

Inertia Dyadic
==============

Using the triple vector product identify
:math:`\bar{a}\times(\bar{b}\times\bar{c}) = \bar{b}(\bar{a}\cdot\bar{c}) -
\bar{c}(\bar{a}\cdot\bar{b}` the inertia vector can be written as:

.. math::

   \bar{I}_a = \sum_{i=1}^\nu m_i \left(\hat{n}_a \bar{r}^{P_i/O}\cdot\bar{r}^{P_i/O} -
   \hat{n}_a\cdot\bar{r}^{P_i/O}\otimes\bar{r}^{P_i/O}\right)

.. math::

   \bar{I}_a =
   \sum_{i=1}^\nu m_i \left(
   \hat{n}_a \cdot \breve{U} \left(\bar{r}^{P_i/O}\right)^2 -
   \hat{n}_a \cdot \bar{r}^{P_i/O} \bar{r}^{P_i/O}
   \right)

.. math::

   \bar{I}_a =
   \hat{n}_a \cdot\sum_{i=1}^\nu m_i \left(
    \breve{U} \left(\bar{r}^{P_i/O}\right)^2 -
    \bar{r}^{P_i/O} \bar{r}^{P_i/O}
   \right) = \hat{n}_a \cdot \breve{I}

We can now remove the dependence on :math:`\hat{n}_a` and arrive at the inertia
dyadic of the set of particles relative to :math:`O`:

.. math::

   \breve{I} :=
   \sum_{i=1}^\nu m_i \left(
   \breve{U} \left(\bar{r}^{P_i/O}\right)^2 -
   \bar{r}^{P_i/O} \otimes \bar{r}^{P_i/O}
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

.. math::

   {}^A \mathbf{H}^{S/S_o} = \breve{I} \cdot {}^A\bar{\omega}^B

.. jupyter-execute::

   w1, w2, w3 = me.dynamicsymbols('omega1, omega2, omega3')
   I = me.inertia(B, Ixx, Iyy, Izz, Ixy, Iyz, Ixz)
   A_w_B = w1*B.x + w2*B.y + w3*B.z

   I.dot(A_w_B)
