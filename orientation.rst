===========
Orientation
===========

Reference Frames
================

In the study of mutlibody dynamics, we are interested in observing motion of
objects in three dimensinal space. This nescitates the concept of a frame of
refenrece, or *Reference Frame*, which is an abstraction that defines a set of
all points in Euclidean space that is carried by the observer of any given
state of motion. Practically, it is useful to image your eye as an observer of
motion. Your eye can orient itself in 3D space to view from motion of objects
from any direction and the motion of objects will appear differently in the
Euclidean space afixed to your eyeball depending on its orientation.

.. _Euclidean space: https://en.wikipedia.org/wiki/Euclidean_space

It is important to note that a reference frame is not equivalent to a
coordinate system. Any number of coordinates systems can be affixed to a
reference frame. We will characterize a reference frame by a dextral set of
orthonormal vectors that relate to each other by the right hand rule.

For a reference frame named :math:`N` we will define the three vectors as
:math:`\hat{n}_x, \hat{n}_y, \hat{n}_z` where these cross products hold:

- :math:`\hat{n}_x \times \hat{n}_y = \hat{n}_z`
- :math:`\hat{n}_y \times \hat{n}_z = \hat{n}_x`
- :math:`\hat{n}_z \times \hat{n}_x = \hat{n}_y`

.. note::

   Unit vectors will be designate using the "hat", e.g. :math:`\hat{v}`. Unit
   vectors have a magnitude, sense, and direction like normal vectors, but
   their magintude is always equal to unity.

These unit vectors are fixed in the reference frame :math:`N`. If a second
reference frame :math:`A` is defined, also with its set of dextral orthonormal
vectors :math:`\hat{a}_x, \hat{a}_y, \hat{a}_z` then we can establish the
relative orietnation of these two reference frames, which may be arbitrarly
oriented with respect to each other.

Take this simple rotation about the parallel :math:`z` unit vectors of
:math:`N` and :math:`A`:

.. todo:: Add a figure

If :math:`A` is rotated with respect to :math:`N` by a simple rotation about
the :math:`\hat{a}_x` through angle :math:`\theta` then these relationships
hold true:

.. math::

   \hat{a}_x & = \cos{\theta} \hat{n}_x + \sin{\theta} \hat{n}_y + 0 \hat{n}_z \\
   \hat{a}_y & = -\sin{\theta} \hat{n}_x + \cos{\theta} \hat{n}_y + 0 \hat{n}_z \\
   \hat{a}_z & = 0 \hat{n}_x + 0 \hat{n}_y + 1 \hat{n}_z

These equations can be written in matrix form:

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

or:

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

The matrix :math:`{}^A\mathbf{C}^N` is the direction cosine matrix, or rotation
matrix, that defines the orientation between two reference frames.

There are some helpful properties of this matrix.

.. jupyter-execute::

   import sympy as sm
   sm.init_printing(use_latex='mathjax')

.. jupyter-execute::

   theta = sm.symbols('theta')

   A_C_N = sm.Matrix([[sm.cos(theta), sm.sin(theta), 0],
                      [-sm.sin(theta), sm.cos(theta), 0],
                      [0, 0, 1]])
   A_C_N

The transpose equals the inverse:

.. jupyter-execute::

   A_C_N.transpose()

.. jupyter-execute::

   sm.trigsimp(A_C_N.inv())

If now :math:`A` is oriented relative to :math:`N` and the three angles
:math:`\alpha_1,\alpha_2,\alpha_3` between :math:`\hat{a}_x` and
:math:`\hat{n}_x,\hat{n}_y,\hat{n}_z`, respectively

.. figure:: fig.svg

Three angles can be found for each unit vector in :math:`A` and then we can
write:

.. math::

  \hat{a}_x & = \cos\alpha_{xx} \hat{n}_x +\cos\alpha_{xy} \hat{n}_y + \cos\alpha_{xz} \hat{n}_z \\
  \hat{a}_y & = \cos\alpha_{yx} \hat{n}_x +\cos\alpha_{yy} \hat{n}_y + \cos\alpha_{yz} \hat{n}_z \\
  \hat{a}_z & = \cos\alpha_{yx} \hat{n}_x +\cos\alpha_{yy} \hat{n}_y + \cos\alpha_{yz} \hat{n}_z

Since we are workign with unit vectors the cosine of the angles between the two
vectors are also equivalent to the dot product between the two vectors.

.. math::

  \hat{a}_x = (\hat{a}_x \cdot \hat{n}_x) \hat{n}_x + (\hat{a}_x \cdot \hat{n}_y) \hat{n}_y + (\hat{a}_x \cdot \hat{n}_z) \hat{n}_z \\
  \hat{a}_y = (\hat{a}_y \cdot \hat{n}_x) \hat{n}_x + (\hat{a}_y \cdot \hat{n}_y) \hat{n}_y + (\hat{a}_y \cdot \hat{n}_z) \hat{n}_z \\
  \hat{a}_x = (\hat{a}_z \cdot \hat{n}_x) \hat{n}_x + (\hat{a}_z \cdot \hat{n}_y) \hat{n}_y + (\hat{a}_z \cdot \hat{n}_z) \hat{n}_z \\

Similar to the simple rotation example above, these can be written in matrix
form:

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
:math:`N` and :math:`A`.

Successive orientations of a series of reference frames provides a convenient
way to manage orientation among more than a pair.


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

then


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

that is

.. math::

   {}^B\mathbf{C}^N
   =
   {}^B\mathbf{C}^A
   {}^A\mathbf{C}^N

Take for example two simple rotatiosn about the colinear :math:`z` unit
vectors.

.. jupyter-execute::

   A_C_N

.. jupyter-execute::

   alpha = sm.symbols('alpha')

   B_C_A = sm.Matrix([[sm.cos(alpha), sm.sin(alpha), 0],
                      [-sm.sin(alpha), sm.cos(alpha), 0],
                      [0, 0, 1]])

   B_C_A

.. jupyter-execute::

   B_C_N = B_C_A*A_C_N
   B_C_N

.. jupyter-execute::

   sm.trigsimp(B_C_N)


.. jupyter-execute::

   import sympy.physics.mechanics as me

.. jupyter-execute::

   N = me.ReferenceFrame('N')

The unit vectors associate with a reference frame are accessed like so:

.. jupyter-execute::

   N.x, N.y, N.z

.. jupyter-execute::

   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')

   N, A, B

.. jupyter-execute::

   A_C_N

.. jupyter-execute::

   A.orient_explicit(N, A_C_N.T)

   A.dcm(N)

.. jupyter-execute::

   B.orient_axis(A, alpha, A.z)

   B.dcm(A)

.. jupyter-execute::

   sm.trigsimp(B.dcm(A)*A.dcm(N))

.. jupyter-execute::

   sm.trigsimp(B.dcm(N))

.. jupyter-execute::

   sm.trigsimp(me.dot(B.x, N.x))

Euler Angles
============
