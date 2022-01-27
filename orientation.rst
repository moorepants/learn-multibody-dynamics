===============================
Orientation of Reference Frames
===============================

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`orientation` or Jupyter notebook:
   :jupyter-download:notebook:`orientation`.

.. jupyter-execute::

   import sympy as sm
   sm.init_printing(use_latex='mathjax')

Reference Frames
================

In the study of mutlibody dynamics, we are interested in observing motion of
objects in three dimensinal space. This nescitates the concept of a frame of
reference, or :term:`reference frame`, which is an abstraction that defines the
set of all points in `Euclidean space`_ that is carried by the observer of any
given state of motion. Practically, it is useful to image your eye as an
observer of motion. Your eye can orient itself in 3D space to view from motion
of objects from any direction and the motion of objects will appear differently
in the set of points associated with the reference frame depending on its
orientation.

.. _Euclidean space: https://en.wikipedia.org/wiki/Euclidean_space

It is important to note that a reference frame is not equivalent to a
coordinate system. Any number of coordinates systems can be affixed to a
reference frame. We will characterize a reference frame by a dextral set of
orthonormal vectors that relate to each other by the right hand rule.

Unit Vectors
============

Vectors have a magnitude, direction, and sense (:math:`\pm`). Unit vectors have
a magnitude of unity. Unit vectors can be fixed, orientation-wise, to a
reference frame. For a reference frame named :math:`N` we will define the three
mutually perpendicular dextral unit vectors as :math:`\hat{n}_x, \hat{n}_y,
\hat{n}_z` where these `cross products`_ hold:

.. _cross products: https://en.wikipedia.org/wiki/Cross_product

.. math::

   \hat{n}_x \times \hat{n}_y & = \hat{n}_z \\
   \hat{n}_y \times \hat{n}_z & = \hat{n}_x \\
   \hat{n}_z \times \hat{n}_x & = \hat{n}_y

.. note::

   Unit vectors will be designate using the "hat", e.g. :math:`\hat{v}`.

These unit vectors are fixed in the reference frame :math:`N`. If a second
reference frame :math:`A` is defined, also with its set of dextral mutually
perpendicular unit vectors :math:`\hat{a}_x, \hat{a}_y, \hat{a}_z` then we can
establish the relative orietnation of these two reference frames based on the
angles among the two frames' unit vectors.

.. note::

   .. _orientation-vector-position:

   .. figure:: orientation-vector-position.svg

      The image on the left and right represent the same set of mutually
      perpendicular dextral unit vectors. Vectors, in general, do not have a
      position and can be drawn anywhere in the reference frame. Drawing them
      with their tails coincident is simply for convenient representation.

Simple Rotations
================

Starting with two reference frames :math:`N` and :math:`A` in which their sets
of unit vectors are aligned, the :math:`A` frame can then be simply rotated
about the common parallel :math:`z` unit vectors of the two frames. We say
"reference frame :math:`A` is rotated with respect to reference frame :math:`N`
about the shared :math:`z` unit vectors through an angle :math:`\theta`.

.. _orientation-simple:

.. figure:: orientation-simple.svg

   View of the parallel :math:`xy` planes of the simply rotated reference
   frames.

From the above figure these relationships between the :math:`a` and :math:`n`
unit vectors can be deduced:

.. math::

   \hat{a}_x & = \cos{\theta} \hat{n}_x + \sin{\theta} \hat{n}_y + 0 \hat{n}_z \\
   \hat{a}_y & = -\sin{\theta} \hat{n}_x + \cos{\theta} \hat{n}_y + 0 \hat{n}_z \\
   \hat{a}_z & = 0 \hat{n}_x + 0 \hat{n}_y + 1 \hat{n}_z

These equations can also be written in matrix form:

.. math::

   \begin{bmatrix}
     \hat{a}_x \\
     \hat{a}_y \\
     \hat{a}_z
   \end{bmatrix}
   =
   \begin{bmatrix}
     \cos{\theta} & \sin{\theta} & 0 \\
     -\sin{\theta} & \cos{\theta} & 0 \\
     0 &  0  & 1
   \end{bmatrix}
   \begin{bmatrix}
     \hat{n}_x \\
     \hat{n}_y \\
     \hat{n}_z
   \end{bmatrix}

This matrix uniquely describes the orientation between the two reference frames
and so we can give it its own variable:

.. math::

   \begin{bmatrix}
     \hat{a}_x \\
     \hat{a}_y \\
     \hat{a}_z
   \end{bmatrix}
   =
   {}^A\mathbf{C}^N
   \begin{bmatrix}
     \hat{n}_x \\
     \hat{n}_y \\
     \hat{n}_z
   \end{bmatrix}

This matrix has an important property, which we will demonstrate with SymPy.
Start by creating the matrix:

.. jupyter-execute::

   theta = sm.symbols('theta')

   A_C_N = sm.Matrix([[sm.cos(theta), sm.sin(theta), 0],
                      [-sm.sin(theta), sm.cos(theta), 0],
                      [0, 0, 1]])
   A_C_N

If we'd like the inverse relationship between the two sets of unit vectors,
then:

.. math::

   \begin{bmatrix}
     \hat{n}_x \\
     \hat{n}_y \\
     \hat{n}_z
   \end{bmatrix}
   =
   \left({}^A\mathbf{C}^N\right)^{-1}
   \begin{bmatrix}
     \hat{a}_x \\
     \hat{a}_y \\
     \hat{a}_z
   \end{bmatrix}

SymPy can find this matrix inverse:

.. jupyter-execute::

   sm.trigsimp(A_C_N.inv())

SymPy can also find the transpose of this matrix;

.. jupyter-execute::

   A_C_N.transpose()

Interestingly, the inverse and the transpose are the same here. It turns out
that this will be generally true for these matrices that describe the
orientation between reference frames. Following the notation convention, this
holds:

.. math::

   {}^N\mathbf{C}^A = \left({}^A\mathbf{C}^N\right)^{-1} = \left({}^A\mathbf{C}^N\right)^T

Direction Cosine Matrix
=======================

If now :math:`A` is oriented relative to :math:`N` and the pairwise angles
between each :math:`A` and :math:`N` mutually perpendicular unit vectors are
measured, an orientation matrix for an arbitrary orientation can be defined.
For example the figure below shows the three angles
:math:`\alpha_{xx},\alpha_{xy},\alpha_{xz}` relating :math:`\hat{a}_x` to each
:math:`N` unit vector.

.. _orientation-three-angles:

.. figure:: orientation-three-angles.svg

   Three angles relating :math:`\hat{a}_x` to the :math:`N` unit vectors.

Similarly to the simple example above, we can write these equations:

.. math::

  \hat{a}_x & = \cos\alpha_{xx} \hat{n}_x +\cos\alpha_{xy} \hat{n}_y + \cos\alpha_{xz} \hat{n}_z \\
  \hat{a}_y & = \cos\alpha_{yx} \hat{n}_x +\cos\alpha_{yy} \hat{n}_y + \cos\alpha_{yz} \hat{n}_z \\
  \hat{a}_z & = \cos\alpha_{yx} \hat{n}_x +\cos\alpha_{yy} \hat{n}_y + \cos\alpha_{yz} \hat{n}_z

Since we are working with mutually perpendicular unit vectors the cosine of the
angle between each pair of unit vectors is equivalent to the dot product
between the two vectors, so this also holds:

.. math::

  \hat{a}_x = (\hat{a}_x \cdot \hat{n}_x) \hat{n}_x + (\hat{a}_x \cdot \hat{n}_y) \hat{n}_y + (\hat{a}_x \cdot \hat{n}_z) \hat{n}_z \\
  \hat{a}_y = (\hat{a}_y \cdot \hat{n}_x) \hat{n}_x + (\hat{a}_y \cdot \hat{n}_y) \hat{n}_y + (\hat{a}_y \cdot \hat{n}_z) \hat{n}_z \\
  \hat{a}_x = (\hat{a}_z \cdot \hat{n}_x) \hat{n}_x + (\hat{a}_z \cdot \hat{n}_y) \hat{n}_y + (\hat{a}_z \cdot \hat{n}_z) \hat{n}_z \\

Now the general :term:`direction cosine matrix` of :math:`A` with respect to
:math:`N` is defined as:

.. math::

   \begin{bmatrix}
     \hat{a}_x \\
     \hat{a}_y \\
     \hat{a}_z
   \end{bmatrix}
   =
   \begin{bmatrix}
     \hat{a}_x \cdot \hat{n}_x &\hat{a}_x \cdot \hat{n}_y & \hat{a}_x \cdot \hat{n}_z \\
     \hat{a}_y \cdot \hat{n}_x &\hat{a}_y \cdot \hat{n}_y & \hat{a}_y \cdot \hat{n}_z \\
     \hat{a}_z \cdot \hat{n}_x &\hat{a}_z \cdot \hat{n}_y & \hat{a}_z \cdot \hat{n}_z
   \end{bmatrix}
   \begin{bmatrix}
     \hat{n}_x \\
     \hat{n}_y \\
     \hat{n}_z
   \end{bmatrix}

where the general direction cosine matrix is then:

.. math::

   {}^A\mathbf{C}^N
   =
   \begin{bmatrix}
     \hat{a}_x \cdot \hat{n}_x &\hat{a}_x \cdot \hat{n}_y & \hat{a}_x \cdot \hat{n}_z \\
     \hat{a}_y \cdot \hat{n}_x &\hat{a}_y \cdot \hat{n}_y & \hat{a}_y \cdot \hat{n}_z \\
     \hat{a}_z \cdot \hat{n}_x &\hat{a}_z \cdot \hat{n}_y & \hat{a}_z \cdot \hat{n}_z
   \end{bmatrix}

This matrix uniquely defines the relative orientation between reference frames
:math:`N` and :math:`A` and inverse is equal to the transpose, as shown above
in the simple example.

Successive orientations
=======================

Successive orientations of a series of reference frames provides a convenient
way to manage orientation among more than a pair. Below an additional reference
frame :math:`B` is shown that is simply rotated with respect to :math:`A` in
the same way that :math:`A` is from :math:`N`.

.. _orientation-simple-successive:

.. figure:: orientation-simple-successive.svg

   Two successive simple rotations.

We know that we can define these two relationships:

.. math::

   \begin{bmatrix}
     \hat{a}_x \\
     \hat{a}_y \\
     \hat{a}_z
   \end{bmatrix}
   =
   {}^A\mathbf{C}^N
   \begin{bmatrix}
     \hat{n}_x \\
     \hat{n}_y \\
     \hat{n}_z
   \end{bmatrix}

.. math::

   \begin{bmatrix}
     \hat{b}_x \\
     \hat{b}_y \\
     \hat{b}_z
   \end{bmatrix}
   =
   {}^B\mathbf{C}^A
   \begin{bmatrix}
     \hat{a}_x \\
     \hat{a}_y \\
     \hat{a}_z
   \end{bmatrix}

Now if you'd like to know the relationship between :math:`B` and :math:`N`, the
above two equatiosn can be manipulated to form:

.. math::

   {}^A\mathbf{C}^B
   \begin{bmatrix}
     \hat{b}_x \\
     \hat{b}_y \\
     \hat{b}_z
   \end{bmatrix}
   =
   {}^A\mathbf{C}^N
   \begin{bmatrix}
     \hat{n}_x \\
     \hat{n}_y \\
     \hat{n}_z
   \end{bmatrix}

Finally we can write:

.. math::

   \begin{bmatrix}
     \hat{b}_x \\
     \hat{b}_y \\
     \hat{b}_z
   \end{bmatrix}
   =
   {}^B\mathbf{C}^A
   {}^A\mathbf{C}^N
   \begin{bmatrix}
     \hat{n}_x \\
     \hat{n}_y \\
     \hat{n}_z
   \end{bmatrix}

showing that the direction cosine matrix between :math:`B` and :math:`N`
results from matrix multiplying the intermediate direction cosine matrices.

.. math::

   {}^B\mathbf{C}^N
   =
   {}^B\mathbf{C}^A
   {}^A\mathbf{C}^N

This holds for any series of successive rotations.

.. math::

   {}^Z\mathbf{C}^A
   =
   {}^Z\mathbf{C}^Y
   {}^Y\mathbf{C}^X
   \ldots
   {}^C\mathbf{C}^B
   {}^B\mathbf{C}^A

Using :numref:`orientation-simple-successive` as an explicit example of this property, we start with
the already defined :math:`{}^A\mathbf{C}^N`:

.. jupyter-execute::

   A_C_N

:math:`{}^B\mathbf{C}^A` can then be defined similarly:

.. jupyter-execute::

   alpha = sm.symbols('alpha')

   B_C_A = sm.Matrix([[sm.cos(alpha), sm.sin(alpha), 0],
                      [-sm.sin(alpha), sm.cos(alpha), 0],
                      [0, 0, 1]])

   B_C_A

Finally, :math:`{}^B\mathbf{C}^N` can be found by matrix multiplication:

.. jupyter-execute::

   B_C_N = B_C_A*A_C_N
   B_C_N

Simplifying these trigonometric expressions shows the likely expected result:

.. jupyter-execute::

   sm.trigsimp(B_C_N)

SymPy Mechanics
===============

As shown above, SymPy nicely handles the formulation of direction cosine
matrices, but SymPy offers a more useful abstraction for these things. The
``sympy.physics.mechanics`` modules includes numerous objects and functions
that ease the bookeeping and mental models needed to manage various aspects of
multibody dynamics. We will import the module consistently as:

.. jupyter-execute::

   import sympy.physics.mechanics as me

``mechanics`` includes a way to define and orient reference frames. To create a
reference frame, use ``ReferenceFrame()`` and provide a name for your frame.

.. jupyter-execute::

   N = me.ReferenceFrame('N')

The dextral mutually perpendicular unit vectors associate with a reference
frame are accessed with ``.x``, ``.y``, and ``.z``, like so:

.. jupyter-execute::

   N.x, N.y, N.z

Using :numref:`orientation-simple-successive` again as an example, we can
define the three reference frames:

.. jupyter-execute::

   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')

   N, A, B

We have already defined the direction cosine matrices for these two successive
rotations. For example:

.. jupyter-execute::

   A_C_N

relates :math:`A` and :math:`N`. ``ReferenceFrame`` objects can be oriented wrt
respect to one another. The ``.orient_explicit()`` method

.. jupyter-execute::

   N.orient_explicit(A, A_C_N)

Now you can ask for the direction cosine matrix of :math:`A` with respect to
:math:`N`, i.e. :math:`{}^A\mathbf{C}^N`, using the ``.dcm()`` method:

.. jupyter-execute::

   A.dcm(N)

.. warning::

   Note very carefully what version of the direction cosine matrix you pass to
   ``.orient_explicit()``. Check its docstring with ``N.orient_explicit?``.

But even better is the ``.orient_axis()`` method. This method allows you do
define simple rotations between reference frames more naturally. You provide
the frame to rotate from, the angle to rotate, and the vector to rotate about.
For example, rotate :math:`B` with respect to :math:`A` through :math:`\alpha`
about :math:`\hat{a}_z` by:

.. jupyter-execute::

   B.orient_axis(A, alpha, A.z)

Now the direction cosine matrix is automatically calculated and is returned
with the ``.dcm()`` method:

.. jupyter-execute::

   B.dcm(A)

The inverse is also defined on ``A``:

.. jupyter-execute::

   A.dcm(B)

So each pair of reference frames are aware of its orientation partner (or
partners). Now that we've established orientaitons between :math:`N` and
:math:`A` and :math:`A` and :math:`B`, we might want to know the relationships
between :math:`B` and :math:`N`. Remember that matrix multiplication of the two
successive direction cosine matrices provides the answer:

.. jupyter-execute::

   sm.trigsimp(B.dcm(A)*A.dcm(N))

But, the answer can also be found by calling ``.dcm()`` with the two reference
frames in question. As long as there is a succussive path between the two
reference frames, this is sufficient for obtaining the desired direction cosine
matrix:

.. jupyter-execute::

   sm.trigsimp(B.dcm(N))

Lastly, recall the general definition of the direction cosine matrix. We showed
that the dot product of pairs of unit vectors give the entries to the direction
cosine matrix. ``mechanics`` has a ``dot()`` function that can calcualte dot
products of two vectors. Using it one two of the unit vector pairs returns the
expected direction cosine matrix entry:

.. jupyter-execute::

   sm.trigsimp(me.dot(B.x, N.x))
