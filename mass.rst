=================
Mass Distribution
=================

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   me.init_vprinting(use_latex='mathjax')

In the prior chapters, we have developed the tools to formulate the kinematics
of points and reference frames. The kinematics are the first of three essential
parts needed to form the equations of motion of a multibody system; the other
two being the mass distribution of and the forces acting on the system.

When a point is associated with a particle of mass :math:`m` or a reference
frame is associated with a rigid body that has some mass distribution,
`Newton's`_ and `Euler's`_ second laws of motion show that the time rate of
change of the linear and angular momenta must be equal to the forces and
torques acting on the particle or rigid body, respectively. The momentum of a
particle is determined by its mass and velocity and the angular momentum of a
rigid body is determined by the distribution of mass and its angular velocity.
In this chapter, we will introduce mass and its distribution.

.. _Newton's: https://en.wikipedia.org/wiki/Newton%27s_laws_of_motion
.. _Euler's: https://en.wikipedia.org/wiki/Euler%27s_laws_of_motion

Particles and Rigid Bodies
==========================

We will introduce and use the concepts of particles and rigid bodies in this
chapter. Both are abstractions of real translating and rotating objects.
Particles are points, that have a location in Euclidean space, which have a
volumetrically infinitesimal mass associated with them. Rigid bodies are made
up of a reference frame, that have orientation, which have an associated
continuous distribution of mass. Each point in the distribution of mass can
translate.

For example, an airplane can be modeled as a rigid body when one is concerned
with both its translation and orientation. This could be useful when
investigating its minimum turn radii and banking angle. But it could also be
modeled as a particle when one is only concerned with its translation; for
example when you observing the motion of the airplane from a location outside
the Earth's atmosphere.

Mass
====

Given a set of :math:`\nu` particles with masses :math:`m_1,\ldots,m_\nu` the
total mass, or *zeroth moment of mass*, of the set is defined as:

.. math::
   :label: eq-zeroth-moment

   m := \sum_{i=1}^\nu m_i

For a solid body with a density :math:`\rho` defined at each point within its
volumetric :math:`V` boundary, the total mass becomes an integral of the
general form:

.. math::
   :label: eq-zeroth-moment-rigid-body

   m := \int_{\textrm{solid}} \rho dV

Mass Center
===========

If each particle in a set of :math:`S` particles is located at positions
:math:`\bar{r}^{P_i/O},\ldots,\bar{r}^{P_\nu/O}` the *first moment of mass* can
be defined as:

.. math::
   :label: eq-first-moment

   \sum_{i=1}^\nu m_i \bar{r}^{P_i/O}

There is then a point :math:`S_o` in which the first mass moment equal to zero:

.. math::
   :label: eq-first-moment-zero

   \sum_{i=1}^\nu m_i \bar{r}^{P_i/S_o} = 0

This point is referred to as the *mass center* (or *center of mass*) of the set
of particles. The mass center can be found with:

.. math::
   :label: mass-center-particles

   \bar{r}^{S_o/O} = \frac{ \sum_{i=1}^\nu m_i \bar{r}^{P_i/O} }{\sum_{i=1}^\nu m_i}

which is the first moment divided by the zeroth moment. For a solid body, this
takes on the integral form:

.. math::
   :label: mass-center-rigid-body

   \bar{r}^{S_o/O} = \frac{ \int_{\textrm{solid}} \rho \bar{r} dV }{ \int_{\textrm{solid}} \rho dV }

The particle form (Eq. :math:numref:`mass-center-particles`) can be calculated
using vectors and scalars in SymPy Mechanics. Here is an example of three
particles each at an arbitrary location relative to :math:`O`:

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

Then, for example, if :math:`m_2=2m_1` and :math:`m_3=3m_1`:

.. jupyter-execute::

   r_O_So.xreplace({m2: 2*m1, m3: 3*m1}).simplify()

Mass Distribution
=================

The inertia, or second moment of mass, describes the distribution of mass
relative to a point about an axis. Inertia characterizes the resistance to
angular acceleration in the same way that mass characterizes the resistance to
linear acceleration. For a set of particles :math:`P_1,\ldots,P_\nu` with
positions :math:`\bar{r}^{P_1/O},\ldots,\bar{r}^{P_\nu/O}` all relative to a
point :math:`O`, the *inertia vector* about the unit vector :math:`\hat{n}_a`
is defined as ([Kane1985]_, pg. 61):

.. math::
   :label: inertia-vector-particles

   \bar{I}_a := \sum_{i=1}^\nu m_i \bar{r}^{P_i/O} \times \left( \hat{n}_a \times
   \bar{r}^{P_i/O}  \right)

.. todo:: Add the rigid body form of the inertia vector.

This vector describes the sum of each mass's contribution to the mass
distribution of the set about a line that is parallel to :math:`\hat{n}_a` and
passes through :math:`O`. Figure :numref:`fig-mass-inertia-vector` shows a
visual representation of this vector for a single particle :math:`P` with mass
:math:`m`.

.. _fig-mass-inertia-vector:
.. figure:: figures/mass-inertia-vector.svg
   :align: center

   Inertia vector for a single particle :math:`P` and its relationship to
   :math:`\hat{n}_a`.

For this single particle the magnitude of :math:`\bar{I}_a` is:

.. math::
   :label: inertia-vector-magnitude

   \left| \bar{I}_a \right| = m \left| \bar{r}^{P/O} \right| ^2 \sin\theta

where :math:`\theta` is angle between :math:`\bar{r}^{P/O}` and
:math:`\hat{n}_a`. We see that :math:`\bar{I}_a` is always perpendicular to
:math:`\bar{r}^{P/O}` and scales with :math:`m`, :math:`| \bar{r}^{P/O} |^2`,
and :math:`\sin\theta`.

If :math:`\hat{n}_a` happens to be parallel to :math:`\bar{r}^{P/O}` then the
magnitude of :math:`\bar{I}_a` is zero. If :math:`\hat{n}_a` is perpendicular
to :math:`\bar{r}^{P/O}` then the magnitude is:

.. math::
   :label: intertia-vector-magnitude-perp

   \left| \bar{I}_a \right| = m \left| \bar{r}^{P/O} \right| ^2

The inertia vector fully describes the distribution of the particles with
respect to :math:`O` about :math:`\hat{n}_a`.

A component of :math:`\bar{I}_a` in the :math:`\hat{n}_b` direction is called
an *inertia scalar* and is defined as ([Kane1985]_, pg 62):

.. math::
   :label: inertia-scalar

   I_{ab} := \bar{I}_{a} \cdot \hat{n}_b

The inertia scalar can be rewritten using Eq.
:math:numref:`inertia-vector-particles`:

.. math::
   :label: eq-product-of-inertia

   I_{ab} =
   \sum_{i=1}^\nu m_i
   \left( \bar{r}^{P_i/O} \times \hat{n}_a \right)
   \cdot
   \left( \bar{r}^{P_i/O} \times \hat{n}_b \right)

This form implies that:

.. math::
   :label: eq-product-of-inertia-equivalence

   I_{ab} = I_{ba}

If :math:`\hat{n}_a = \hat{n}_b` then this inertia scalar is called a *moment
of inertia* and if :math:`\hat{n}_a \neq \hat{n}_b` it is called a *product of
inertia*. Moments of inertia describe the mass distribution about a single axis
whereas products of inertia describe the mass distribution relative to two
axes.

When :math:`\hat{n}_a = \hat{n}_b` Eq. :math:numref:`eq-product-of-inertia`
reduces to the moment of inertia:

.. math::
   :label: eq-moment-of-inertia

   I_{aa} =
   \sum_{i=1}^\nu m_i
   \left( \bar{r}^{P_i/O} \times \hat{n}_a \right)^2

It is common to define the *radius of gyration* :math:`k_{aa}`, which is the
radius of a ring that has the same moment of inertia as the set of particles or
rigid body. The radius of gyration about a line through :math:`O` parallel to
:math:`\hat{n}_a` is defined as:

.. math::

   k_{aa} := \sqrt{\frac{I_{aa}}{m}}

Inertia Matrix
==============

For mutually perpendicular unit vectors fixed in reference frame :math:`A`, the
moments of inertia with respect to :math:`O` about each unit vector and the
products of inertia among the pairs of perpendicular unit vectors can be
computed. This, in general, results in nine inertia scalars that describe the
mass distribution of a set of particles or a rigid body in 3D space. These
scalars are typically presented as a symmetric *inertia matrix* (also called an
*inertia tensor*) that takes this form:

.. math::
   :label: eq-inertia-matrix

   \begin{bmatrix}
    I_{xx} & I_{xy} & I_{xz} \\
    I_{yx} & I_{yy} & I_{yz} \\
    I_{zx} & I_{zy} & I_{zz}
   \end{bmatrix}_A

where the moments of inertia are on the diagonal and the products of inertia
are the off diagonal entries. Eq.
:math:numref:`eq-product-of-inertia-equivalence` holds for the products of
inertia, i.e. :math:`I_{xy}=I_{yx}`, :math:`I_{xz}=I_{zx}`, and
:math:`I_{yz}=I_{zy}`, and the subscript :math:`A` indicates that these scalars
are relative to unit vectors :math:`\hat{a}_x,\hat{a}_y,\hat{a}_z`.

This matrix (or second order tensor) is similar to the vectors (or first order
tensors) we've already worked with:

.. math::
   :label: eq-column-vector

   \begin{bmatrix}
   v_1 \\
   v_2 \\
   v_3
   \end{bmatrix}_A

Recall that we have a notation for writing such a vector that allows us to
combine components expressed in different reference frames:

.. math::

   v_1\hat{a}_x + v_2\hat{a}_y + v_3\hat{a}_z

There also exists an analogous form for second order tensors that are
associated with different reference frames called a dyadic_.

.. _dyadic: https://en.wikipedia.org/wiki/Dyadics

Dyadics
=======

If we introduce the `outer product`_ operator between two vectors we see that
it generates a matrix akin to the inertia matrix above.

.. math::
   :label: eq-vector-outer-product

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

In SymPy Mechanics outer products can be taken between two vectors to create
the dyadic :math:`\breve{Q}` using
:external:py:func:`~sympy.physics.vector.functions.outer`:

.. jupyter-execute::

   v1, v2, v3 = sm.symbols('v1, v2, v3')
   w1, w2, w3 = sm.symbols('w1, w2, w3')

   A = me.ReferenceFrame('A')

   v = v1*A.x + v2*A.y + v3*A.z
   w = w1*A.x + w2*A.y + w3*A.z

   Q = me.outer(v, w)
   Q

The result is not the matrix form shown in Eq.
:math:numref:`eq-vector-outer-product`, but instead the result is a dyadic. The
dyadic is the analogous form for second order tensors as what we've been using
for first order tensors. If the matrix form is needed, it can be found with
:external:py:meth:`~sympy.physics.vector.dyadic.Dyadic.to_matrix`:

.. jupyter-execute::

   Q.to_matrix(A)

The dyadic is made up of scalars multiplied by unit dyads. Examples of unit
dyads are:

.. jupyter-execute::

   me.outer(A.x, A.x)

Unit dyads correspond to unit entries in the 3x3 matrix:

.. jupyter-execute::

   me.outer(A.x, A.x).to_matrix(A)

Unit dyads are analogous to unit vectors. Here is another example:

.. jupyter-execute::

   me.outer(A.y, A.z)

.. jupyter-execute::

   me.outer(A.y, A.z).to_matrix(A)

These unit dyads can be formed from any unit vectors. This is convenient
because we can create dyadics, just like vectors, which are made up of
components in different reference frames. For example:

.. jupyter-execute::

   theta = sm.symbols("theta")

   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')

   B.orient_axis(A, theta, A.x)

   P = 2*me.outer(B.x, B.x) + 3*me.outer(A.x, B.y) + 4*me.outer(B.z, A.z)
   P

The dyadic :math:`\breve{P}` can be expressed in unit dyads of :math:`A`

.. jupyter-execute::

   P.express(A)

.. jupyter-execute::

   P.to_matrix(A)

or :math:`B`: :

.. jupyter-execute::

   P.express(B)

.. jupyter-execute::

   P.to_matrix(B)

The *unit dyadic* is defined as:

.. math::
   :label: eq-unit-dyadic

   \breve{U} :=
   \hat{a}_x \otimes \hat{a}_x +
   \hat{a}_y \otimes \hat{a}_y +
   \hat{a}_z \otimes \hat{a}_z

.. todo:: I need a notation to distinguish a unit dyadic like we do with unit
   vectors and vectors.

The unit dyadic can be created with SymPy:

.. jupyter-execute::

   U = me.outer(A.x, A.x) + me.outer(A.y, A.y) + me.outer(A.z, A.z)
   U

and it represents the identity matrix in :math:`A`:

.. jupyter-execute::

   U.to_matrix(A)

Note that the unit dyadic is the same when expressed in any reference frame:

.. jupyter-execute::

   U.express(B).simplify()

.. todo:: ReferenceFrame should have an attribute that returns the unit dyadic
   (or dyads).

Properties of Dyadics
=====================

Dyadics have similar properties as vectors but are not necessarily commutative.

- Scalar multiplication: :math:`\alpha(\bar{u}\otimes\bar{v}) = \alpha\bar{u}\otimes\bar{v} = \bar{u}\otimes\alpha\bar{v}`
- Distributive: :math:`\bar{u}\otimes(\bar{v} + \bar{w}) = \bar{u}\otimes\bar{v} + \bar{u}\otimes\bar{w}`
- Left and right dot product with a vector (results in a vector):

  - :math:`\bar{u}\cdot(\bar{v}\otimes\bar{w}) = (\bar{u}\cdot\bar{v})\bar{w}`
  - :math:`(\bar{u}\otimes\bar{v})\cdot\bar{w} = \bar{u}(\bar{v}\cdot\bar{w})`

- Left and right cross product with a vector (results in a dyadic):

  - :math:`\bar{u}\times(\bar{v}\otimes\bar{w}) = (\bar{u}\times\bar{v})\otimes\bar{w}`
  - :math:`(\bar{u}\otimes\bar{v})\times\bar{w} = \bar{u}\otimes(\bar{v}\times\bar{w})`

- Dot products between arbitrary vectors and arbitrary dyadics are not
  commutative: :math:`\breve{V}\cdot\bar{u} \neq \bar{u}\cdot\breve{V}`
- Dot products between arbitrary vectors and the unit dyadic are commutative
  and result in the vector itself: :math:`\breve{U}\cdot\bar{v} =
  \bar{v}\cdot\breve{U} = \bar{v}`

Inertia Dyadic
==============

Previously we defined the inertia vector. Using the `vector triple product`_
identity: :math:`\bar{a}\times(\bar{b}\times\bar{c}) =
\bar{b}(\bar{a}\cdot\bar{c}) - \bar{c}(\bar{a}\cdot\bar{b})`, the inertia
vector can be written as ([Kane1985]_, pg 68):

.. _vector triple product: https://en.wikipedia.org/wiki/Triple_product#Vector_triple_product

.. math::
   :label: eq-apply-triple-vec-product

   \bar{I}_a & = \sum_{i=1}^\nu m_i \bar{r}^{P_i/O} \times \left( \hat{n}_a \times \bar{r}^{P_i/O}  \right) \\
   \bar{I}_a & = \sum_{i=1}^\nu m_i
   \left[\hat{n}_a \left( \bar{r}^{P_i/O} \cdot \bar{r}^{P_i/O} \right) -
   \bar{r}^{P_i/O} \left( \bar{r}^{P_i/O} \cdot \hat{n}_a \right) \right]

Now by introducing a unit dyadic, we can write:

.. math::

   \bar{I}_a =
   \sum_{i=1}^\nu m_i \left[
   \left|\bar{r}^{P_i/O}\right|^2 \hat{n}_a \cdot \breve{U}  -
   \hat{n}_a \cdot \left(\bar{r}^{P_i/O} \otimes \bar{r}^{P_i/O}\right)
   \right]

:math:`\hat{n}_a` can be pulled out of the summation:

.. math::

   \bar{I}_a =
   \hat{n}_a \cdot
   \sum_{i=1}^\nu m_i \left(
   \left|\bar{r}^{P_i/O}\right|^2 \breve{U}  -
   \bar{r}^{P_i/O} \otimes \bar{r}^{P_i/O}
   \right)

The *inertia dyadic* :math:`\breve{I}` of a set of :math:`S` particles relative
to :math:`O` is now defined as:

.. math::
   :label: eq-inertia-dyadic

   \breve{I}^{S/O} :=
   \sum_{i=1}^\nu m_i \left(
   \left|\bar{r}^{P_i/O}\right|^2 \breve{U}  -
   \bar{r}^{P_i/O} \otimes \bar{r}^{P_i/O}
   \right)

where:

.. math::

   \bar{I}_a = \hat{n}_a \cdot \breve{I}^{S/O}

.. todo:: Add the rigid body form of the inertia dyadic equation.

Note that we have now described the inertia of the set of particles without
needing to specify a vector :math:`\hat{n}_a`. This inertia dyadic contains the
complete description of inertia with respect to point :math:`O` about any axis.
The vectors and dyadics in Eq. :math:numref:`eq-inertia-dyadic` can be written
in terms of any reference frame unit vectors or unit dyads, respectively.

In SymPy Mechanics, simple inertia dyadics in terms of the unit vectors of a
single reference frame can quickly be generated with
:external:py:func:`~sympy.physics.mechanics.functions.inertia`. For example:

.. jupyter-execute::

   Ixx, Iyy, Izz = sm.symbols('I_{xx}, I_{yy}, I_{zz}')
   Ixy, Iyz, Ixz = sm.symbols('I_{xy}, I_{yz}, I_{xz}')

   I = me.inertia(A, Ixx, Iyy, Izz, ixy=Ixy, iyz=Iyz, izx=Ixz)
   I

.. jupyter-execute::

   I.to_matrix(A)

This inertia dyadic can easily be expressed relative to another reference frame
if the orientation is defined:

.. jupyter-execute::

   sm.trigsimp(I.to_matrix(B))

This is equivalent to the matrix transform to express an inertia matrix in
other reference frame:

.. math::
   :label: eq-inertia-transform

   {}^B\mathbf{C}^A \ \mathbf{I} \ {}^A\mathbf{C}^B

.. jupyter-execute::

   sm.trigsimp(B.dcm(A)*I.to_matrix(A)*A.dcm(B))

.. note:: Angular Momentum

   The angular momentum of a rigid body :math:`B` in reference frame :math:`A`
   about point :math:`O` is defined as:

   .. math::
      :label: eq-angular-momentum

      {}^A \mathbf{H}^{B/O} := \breve{I}^{B/O} \cdot {}^A\bar{\omega}^B

   If the point is instead the mass center of :math:`B`, point :math:`B_o`,
   then the inertia dyadic is the *central inertia dyadic* and the result is
   the *central angular momentum* in :math:`A` is:

   .. math::
      :label: eq-central-angular-momentum

      {}^A \mathbf{H}^{B/B_o} = \breve{I}^{B/B_o} \cdot {}^A\bar{\omega}^B

   Here is an example of calculating the angular momentum in SymPy Mechanics:

   .. jupyter-execute::

      Ixx, Iyy, Izz = sm.symbols('I_{xx}, I_{yy}, I_{zz}')
      Ixy, Iyz, Ixz = sm.symbols('I_{xy}, I_{yz}, I_{xz}')
      w1, w2, w3 = me.dynamicsymbols('omega1, omega2, omega3')

      B = me.ReferenceFrame('B')

      I = me.inertia(B, Ixx, Iyy, Izz, Ixy, Iyz, Ixz)

      A_w_B = w1*B.x + w2*B.y + w3*B.z

      I.dot(A_w_B)

Parallel Axis Theorem
=====================

If you know the central inertia dyadic of a rigid body :math:`B` (or
equivalently a set of particles) about its mass center :math:`B_o` then it is
possible to calculate the inertia dyadic about any other point :math:`O`. To do
so, you must account for the inertial contribution due to the distance between
the points :math:`O` and :math:`B_o`. This is done with the `parallel axis
theorem`_ ([Kane1985]_, pg. 70):

.. math::
   :label: eq-parallel-axis-theorem

   \breve{I}^{B/O} = \breve{I}^{B/B_o} + \breve{I}^{B_o/O}

.. _parallel axis theorem: https://en.wikipedia.org/wiki/Parallel_axis_theorem

The last term is the inertia of a particle with mass :math:`m` (total mass of
the body or set of particles) located at the mass center about point :math:`O`.

.. math::
   :label: eq-parallel-axis-theorem

   \breve{I}^{B_o/O} = m \left(
   \left|\bar{r}^{B_o/O}\right|^2 \breve{U}  -
   \bar{r}^{B_o/O} \otimes \bar{r}^{B_o/O}
   \right)

When :math:`B_o` is displaced from point :math:`O` by three Cartesian distances
:math:`d_x,d_y,d_z` the general form of the last term in Eq.
:math:numref:`eq-parallel-axis-theorem` can be found:

.. jupyter-execute::

   dx, dy, dz, m = sm.symbols('d_x, d_y, d_z, m')

   N = me.ReferenceFrame('N')

   r_O_Bo = dx*N.x + dy*N.y + dz*N.z

   U = me.outer(N.x, N.x) + me.outer(N.y, N.y) + me.outer(N.z, N.z)

   I_Bo_O = m*(me.dot(r_O_Bo, r_O_Bo)*U - me.outer(r_O_Bo, r_O_Bo))
   I_Bo_O

The matrix form of this dyadic shows the typical presentation of the parallel
axis addition term:

.. jupyter-execute::

   I_Bo_O.to_matrix(N)

Principal Axes and Moments of Inertia
=====================================

If the inertia vector :math:`\bar{I}_a` with respect to point :math:`O` is
parallel to its unit vector :math:`\hat{n}_a` then the line through :math:`O`
and parallel to :math:`\hat{n}_a` is called a *principal axis* of the set of
particles or rigid body. The plane that is normal to :math:`\hat{n}_a` is
called a *principal plane*. The moment of inertia about this principal axis is
called a *principal moment of inertia*. The consequence of :math:`\bar{I}_a`
being parallel to :math:`\hat{n}_a` is that the products of inertia are all
zero. The *principal inertia dyadic* can then be written as so:

.. math::
   :label: eq-principal-inertia-dyadic

   \breve{I}^{B/O} =
   I_{11} \hat{b}_1 \otimes \hat{b}_1 +
   I_{22} \hat{b}_2 \otimes \hat{b}_2 +
   I_{33} \hat{b}_3 \otimes \hat{b}_3

where :math:`\hat{b}_1,\hat{b}_2,\hat{b}_3` are mutually perpendicular unit
vectors in :math:`B` that are each parallel to a principal axis and
:math:`I_{11},I_{22},I_{33}` are all principal moments of inertia.

Geometrically symmetric objects with uniform mass density have principal planes
that are perpendicular with the planes of symmetry of the geometry. But there
also exist unique principal axes for non-symmetric and non-uniform density
objects.

The principal axes and their associated principal moments of inertia can be
found by solving the eigenvalue problem. The eigenvalues of an arbitrary
inertia matrix are the principal moments of inertia and the eigenvectors are
the unit vectors parallel to the mutually perpendicular principal axes.
Recalling that the inertia matrix is a symmetric matrix of real numbers, we
know then that it is Hermitian and therefore all its eigenvalues are real.
Symmetric matrices are also diagonalizable and the eigenvectors will then be
orthonormal.

.. warning::

   Finding the eigenvalues of a 3x3 matrix require finding the roots of the
   `cubic equation`_. It is possible to find the symbolic solution, but it is
   not a simple result. Unless you really need the symbolic result, it is best
   to solve for principal axes and moments of inertia numerically.

.. _cubic equation: https://en.wikipedia.org/wiki/Cubic_equation

Here is an example of finding the principal axes and associated moments of
inertia with SymPy:

.. jupyter-execute::

   I = sm.Matrix([[1.0451, 0.0, -0.1123],
                  [0.0, 2.403, 0.0],
                  [-0.1123, 0.0, 1.8501]])
   I

The :external:py:meth:`~sympy.matrices.matrices.MatrixEigen.eigenvects` method on a
SymPy matrix returns a list of tuples that each contain ``(eigenvalue,
multiplicity, eigenvector)``:

.. jupyter-execute::

   I.eigenvects()
