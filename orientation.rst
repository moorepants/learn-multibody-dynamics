===============================
Orientation of Reference Frames
===============================

.. note::

   You can download this example as a Python script:
   :jupyter-download-script:`orientation` or Jupyter Notebook:
   :jupyter-download-notebook:`orientation`.

.. jupyter-execute::

   import sympy as sm
   sm.init_printing(use_latex='mathjax')

Learning Objectives
===================

After completing this chapter readers will be able to:

- Define a reference frame with associated unit vectors.
- Define a direction cosine matrix between two oriented reference frames.
- Derive direction cosine matrices for simple rotations.
- Derive direction cosine matrices for successive rotations.
- Manage orientation and direction cosine matrices with SymPy.
- Rotate reference frames using Euler Angles.

Reference Frames
================

In the study of multibody dynamics, we are interested in observing motion of
connected and interacting objects in three dimensional space. This observation
necessitates the concept of a *frame of reference*, or reference frame. A
reference frame is an abstraction which we define as the set of all points in
`Euclidean space`_ that are carried by and fixed to the observer of any given
state of motion. Practically speaking, it is useful to imagine your eye as an
observer of motion. Your eye can orient itself in 3D space to view the motion
of objects from any direction and the motion of objects will appear differently
in the set of points associated with the reference frame attached to your eye
depending on your eye's orientation.

.. _Euclidean space: https://en.wikipedia.org/wiki/Euclidean_space

It is important to note that a reference frame is not equivalent to a
coordinate system. Any number of coordinate systems (e.g., Cartesian or
spherical) can be used to describe the motion of points or objects in a
reference frame. The coordinate system offers a system of measurement in a
reference frame. We will characterize a reference frame by a right-handed_ set
of mutually perpendicular unit vectors that can be used to describe its
orientation relative to other reference frames and we will align a Cartesian
coordinate system with the unit vectors to allow for easy measurement of points
fixed or moving in the reference frame.

.. _right-handed: https://en.wikipedia.org/wiki/Right-hand_rule

Unit Vectors
============

Vectors have a magnitude, direction, and sense (:math:`\pm`) but notably not a
position. Unit vectors have a magnitude of 1. Unit vectors can be fixed,
orientation-wise, to a reference frame. For a reference frame named :math:`N`
we will define the three mutually perpendicular unit vectors as
:math:`\hat{n}_x, \hat{n}_y, \hat{n}_z` where these right-handed `cross
products`_ hold:

.. _cross products: https://en.wikipedia.org/wiki/Cross_product

.. math::
   :label: right-hand-rule

   \hat{n}_x \times \hat{n}_y & = \hat{n}_z \\
   \hat{n}_y \times \hat{n}_z & = \hat{n}_x \\
   \hat{n}_z \times \hat{n}_x & = \hat{n}_y

.. note::

   Unit vectors will be designated using the "hat", e.g. :math:`\hat{v}`.

These unit vectors are fixed in the reference frame :math:`N`. If a second
reference frame :math:`A` is defined, also with its set of right-handed
mutually perpendicular unit vectors :math:`\hat{a}_x, \hat{a}_y, \hat{a}_z`
then we can establish the relative orientation of these two reference frames
based on the angles among the two frames' unit vectors.

.. _orientation-vector-position:

.. figure:: figures/orientation-vector-position.svg
   :align: center

   The image on the left and right represent the same set of right-handed
   mutually perpendicular unit vectors. Vectors, in general, do not have a
   position and can be drawn anywhere in the reference frame. Drawing them with
   their tails coincident is simply done for convenience.

Simple Orientations
===================

Starting with two reference frames :math:`N` and :math:`A` in which their sets
of unit vectors are initially aligned, the :math:`A` frame can then be simply
oriented about the common parallel :math:`z` unit vectors of the two frames. We
then say "reference frame :math:`A` is oriented with respect to reference frame
:math:`N` about the shared :math:`z` unit vectors through an angle
:math:`\theta`. A visual representation of this orientation looking from the
direction of the positive :math:`z` unit vector is:

.. _orientation-simple:

.. figure:: figures/orientation-simple.svg
   :align: center

   View of the parallel :math:`xy` planes of the simply oriented reference
   frames.

From the above figure these relationships between the :math:`\hat{a}` and
:math:`\hat{n}` unit vectors can be deduced:

.. math::
   :label: simple-orientation-unit-vector-relation

   \hat{a}_x & = \cos{\theta} \hat{n}_x + \sin{\theta} \hat{n}_y + 0 \hat{n}_z \\
   \hat{a}_y & = -\sin{\theta} \hat{n}_x + \cos{\theta} \hat{n}_y + 0 \hat{n}_z \\
   \hat{a}_z & = 0 \hat{n}_x + 0 \hat{n}_y + 1 \hat{n}_z

These equations can also be written in a matrix form:

.. math::
   :label: simple-orientation-unit-vector-relation-matrix

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
and so we give it its own variable:

.. math::
   :label: simple-orient-dcm

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

This matrix :math:`{}^A\mathbf{C}^N` maps vectors expressed in the :math:`N`
frame to vectors expressed in the :math:`A` frame. This matrix has an important
property, which we will demonstrate with SymPy.  Start by creating the matrix:

.. jupyter-execute::

   theta = sm.symbols('theta')

   A_C_N = sm.Matrix([[sm.cos(theta), sm.sin(theta), 0],
                      [-sm.sin(theta), sm.cos(theta), 0],
                      [0, 0, 1]])
   A_C_N

If we'd like the inverse relationship between the two sets of unit vectors and
:math:`{}^A\mathbf{C}^N` is invertible, then:

.. math::
   :label: dcm-inverse

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

Notably, the inverse and the transpose are the same here. This indicates that
this matrix is a special `orthogonal matrix`_. All matrices that describe the
orientation between reference frames are orthogonal matrices. Following the
notation convention, this holds:

.. math::
   :label: dcm-inverse-transpose

   {}^N\mathbf{C}^A = \left({}^A\mathbf{C}^N\right)^{-1} = \left({}^A\mathbf{C}^N\right)^T

.. _orthogonal matrix: https://en.wikipedia.org/wiki/Orthogonal_matrix

.. admonition:: Exercise

   Write :math:`{}^A\mathbf{C}^N` for simple rotations about both the shared
   :math:`\hat{n}_x` and :math:`\hat{a}_x` and shared :math:`\hat{n}_y` and
   :math:`\hat{a}_y` axes, rotating :math:`A` with respect to :math:`N` through
   angle :math:`\theta`.

.. admonition:: Solution
   :class: dropdown

   For a :math:`x` orientation:

   .. math::

      \begin{bmatrix}
        1 &  0  & 0 \\
        0 & \cos{\theta} & \sin{\theta} \\
        0 & -\sin{\theta} & \cos{\theta}
      \end{bmatrix}

   For a :math:`y` orientation:

   .. math::

      \begin{bmatrix}
        \cos{\theta} & 0 & -\sin{\theta} \\
        0 &  1  & 0 \\
        \sin{\theta} & 0 & \cos{\theta}
      \end{bmatrix}

Direction Cosine Matrices
=========================

If now :math:`A` is oriented relative to :math:`N` and the pairwise angles
between each :math:`\hat{a}` and :math:`\hat{n}` mutually perpendicular unit
vectors are measured, a matrix for an arbitrary orientation can be defined.
For example, the figure below shows the three angles
:math:`\alpha_{xx},\alpha_{xy},\alpha_{xz}` relating :math:`\hat{a}_x` to each
:math:`\hat{n}` unit vector.

.. _orientation-three-angles:

.. figure:: figures/orientation-three-angles.svg
   :align: center

   Three angles relating :math:`\hat{a}_x` to the unit vectors of :math:`N`.

.. todo:: Use double sided arrows on the above figure so a sign of the angle is
   not implied.

Similar to the simple example above, we can write these equations if the
:math:`\alpha_y` and :math:`\alpha_z` angles relate the :math:`\hat{a}_y` and
:math:`\hat{a}_z` unit vectors to those of :math:`N`:

.. math::
   :label: direction-cosine-unit-vectors

   \hat{a}_x & = \cos\alpha_{xx} \hat{n}_x +\cos\alpha_{xy} \hat{n}_y + \cos\alpha_{xz} \hat{n}_z \\
   \hat{a}_y & = \cos\alpha_{yx} \hat{n}_x +\cos\alpha_{yy} \hat{n}_y + \cos\alpha_{yz} \hat{n}_z \\
   \hat{a}_z & = \cos\alpha_{zx} \hat{n}_x +\cos\alpha_{zy} \hat{n}_y + \cos\alpha_{zz} \hat{n}_z

Since we are working with unit vectors the cosine of the angle between each
pair of vectors is equivalent to the dot product between the two vectors, so
this also holds:

.. math::
   :label:

   \hat{a}_x = (\hat{a}_x \cdot \hat{n}_x) \hat{n}_x + (\hat{a}_x \cdot \hat{n}_y) \hat{n}_y + (\hat{a}_x \cdot \hat{n}_z) \hat{n}_z \\
   \hat{a}_y = (\hat{a}_y \cdot \hat{n}_x) \hat{n}_x + (\hat{a}_y \cdot \hat{n}_y) \hat{n}_y + (\hat{a}_y \cdot \hat{n}_z) \hat{n}_z \\
   \hat{a}_x = (\hat{a}_z \cdot \hat{n}_x) \hat{n}_x + (\hat{a}_z \cdot \hat{n}_y) \hat{n}_y + (\hat{a}_z \cdot \hat{n}_z) \hat{n}_z \\

Now the matrix relating the orientation of :math:`A` with respect to :math:`N`
can be formed:

.. math::
   :label: dcm-dot-full-eq

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

where

.. math::
   :label: dcm-dot-products

   {}^A\mathbf{C}^N
   =
   \begin{bmatrix}
     \hat{a}_x \cdot \hat{n}_x &\hat{a}_x \cdot \hat{n}_y & \hat{a}_x \cdot \hat{n}_z \\
     \hat{a}_y \cdot \hat{n}_x &\hat{a}_y \cdot \hat{n}_y & \hat{a}_y \cdot \hat{n}_z \\
     \hat{a}_z \cdot \hat{n}_x &\hat{a}_z \cdot \hat{n}_y & \hat{a}_z \cdot \hat{n}_z
   \end{bmatrix}

We call :math:`{}^A\mathbf{C}^N` the "direction cosine matrix" as a general
description of the relative orientation of two reference frames. This matrix
uniquely defines the relative orientation between reference frames :math:`N`
and :math:`A`, it is invertible, and its inverse is equal to the transpose, as
shown above in the simple example. The determinant of the matrix is also always 1, to ensure both associated frames are right-handed. The direction cosine matrix found in the
prior section for a simple orientation is a specific case of this more general
definition. The direction cosine matrix is also referred to as a "rotation
matrix" or "orientation matrix" in some texts.

.. todo:: Create a simple exercise involving applying the definition.

Successive Orientations
=======================

Successive orientations of a series of reference frames provides a convenient
way to manage orientation among more than a single pair. Below, an additional
auxiliary reference frame :math:`B` is shown that is simply oriented with
respect to :math:`A` in the same way that :math:`A` is from :math:`N` above in
the prior section.

.. _orientation-simple-successive:

.. figure:: figures/orientation-simple-successive.svg
   :align: center

   Two successive simple orientations through angles :math:`\theta` and then
   :math:`\alpha` for frames :math:`A` and :math:`B`, respectively.

We know from the prior sections that we can define these two relationships
between each pair of reference frames as follows:

.. math::
   :label: dcm-suc-01

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
   :label: dcm-suc-02

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

Now, substitute :math:numref:`dcm-suc-01` into :math:numref:`dcm-suc-02` to
get:

.. math::
   :label: dcm-multiply-eq

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
   :label: dcm-simple-relationship

   {}^B\mathbf{C}^N
   =
   {}^B\mathbf{C}^A
   {}^A\mathbf{C}^N

This holds for any series of general three dimensional successive orientations
and the relation is shown in the following theorem:

.. math::
   :label: dcm-theorem

   {}^Z\mathbf{C}^A
   =
   {}^Z\mathbf{C}^Y
   {}^Y\mathbf{C}^X
   \ldots
   {}^C\mathbf{C}^B
   {}^B\mathbf{C}^A

where frames :math:`A` through :math:`Z` are succesively oriented.

Using :numref:`orientation-simple-successive` as an explicit example of this
property, we start with the already defined :math:`{}^A\mathbf{C}^N`:

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

Simplifying these trigonometric expressions shows the expected result:

.. jupyter-execute::

   sm.trigsimp(B_C_N)

.. admonition:: Exercise

   If you are given :math:`{}^B\mathbf{C}^N` and :math:`{}^A\mathbf{C}^N` from
   the prior example, how would you find :math:`{}^A\mathbf{C}^B`?

.. admonition:: Solution
   :class: dropdown

   .. math::

      {}^B\mathbf{C}^N
      &=
      {}^B\mathbf{C}^A
      {}^A\mathbf{C}^N \\
      {}^A\mathbf{C}^B
      &=
      \left({}^B\mathbf{C}^N
      \left({}^A\mathbf{C}^N\right)^T\right)^T
      =
      {}^A\mathbf{C}^N\left({}^B\mathbf{C}^N\right)^T

SymPy Mechanics
===============

As shown above, SymPy nicely handles the formulation of direction cosine
matrices, but SymPy also offers a more useful tool for tracking orientation
among reference frames. The :external:py:mod:`sympy.physics.mechanics`
:term:`module` includes numerous objects and functions that ease the
bookkeeping and mental models needed to manage various aspects of multibody
dynamics. We will import the module as in this text:

.. jupyter-execute::

   import sympy.physics.mechanics as me

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

:external:py:mod:`sympy.physics.mechanics` includes a way to define and orient
reference frames. To create a reference frame, use
:external:py:class:`~sympy.physics.vector.frame.ReferenceFrame` and provide a
name for your frame as a string.

.. jupyter-execute::

   N = me.ReferenceFrame('N')

The right-handed mutually perpendicular unit vectors associated with a
reference frame are accessed with the :term:`attributes <attribute>` ``.x``,
``.y``, and ``.z``, like so:

.. jupyter-execute::

   N.x, N.y, N.z

Using :numref:`orientation-simple-successive` again as an example, we can
define all three reference frames by additionally creating :math:`A` and
:math:`B`:

.. jupyter-execute::

   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')

   N, A, B

We have already defined the direction cosine matrices for these two successive
orientations. For example:

.. jupyter-execute::

   A_C_N

relates :math:`A` and :math:`N`. ``ReferenceFrame`` objects can be oriented
with respect to one another. The
:external:py:meth:`~sympy.physics.vector.frame.ReferenceFrame.orient_explicit`
:term:`method` allows you to set the direction cosine matrix between two frames
explicitly:

.. jupyter-execute::

   N.orient_explicit(A, A_C_N)

.. warning::

   Note very carefully what version of the direction cosine matrix you pass to
   ``.orient_explicit()``. Check its docstring with ``N.orient_explicit?``.

Now you can ask for the direction cosine matrix of :math:`A` with respect to
:math:`N`, i.e. :math:`{}^A\mathbf{C}^N`, using the
:external:py:meth:`~sympy.physics.vector.frame.ReferenceFrame.dcm` method:

.. jupyter-execute::

   A.dcm(N)

The direction cosine matrix of :math:`N` with respect to :math:`A` is found by
reversing the order of the arguments:

.. jupyter-execute::

   N.dcm(A)

.. admonition:: Exercise

   Orient reference frame :math:`D` with respect to :math:`F` with a simple
   rotation about :math:`y` through angle :math:`\beta` and set this
   orientation with
   :external:py:meth:`~sympy.physics.vector.frame.ReferenceFrame.orient_explicit`.

.. admonition:: Solution
   :class: dropdown

   .. jupyter-execute::

      beta = sm.symbols('beta')

      D = me.ReferenceFrame('D')
      F = me.ReferenceFrame('F')

      F_C_D = sm.Matrix([[sm.cos(beta), 0, -sm.sin(beta)],
                         [0, 1, 0],
                         [sm.sin(beta), 0, sm.cos(beta)]])

      F.orient_explicit(D, F_C_D.transpose())

      F.dcm(D)

:external:py:meth:`~sympy.physics.vector.frame.ReferenceFrame.orient_explicit`
requires you to form the direction cosine matrix yourself, but there are also
methods that relieve you of that necessity. For example,
:external:py:meth:`~sympy.physics.vector.frame.ReferenceFrame.orient_axis`
allows you to define simple orientations between reference frames more
naturally. You provide the frame to orient from, the angle to orient through,
and the vector to orient about and the correct direction cosine matrix will be
formed. As an example, orient :math:`B` with respect to :math:`A` through
:math:`\alpha` about :math:`\hat{a}_z` by:

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
partners).

Now that we've established orientations between :math:`N` and :math:`A` and
:math:`A` and :math:`B`, we might want to know the relationships between
:math:`B` and :math:`N`. Remember that matrix multiplication of the two
successive direction cosine matrices provides the answer:

.. jupyter-execute::

   sm.trigsimp(B.dcm(A)*A.dcm(N))

But, the answer can also be found by calling
:external:py:meth:`~sympy.physics.vector.frame.ReferenceFrame.dcm` with just
the two reference frames in question, :math:`B` and :math:`N`. As long as there
is a successive path of intermediate, or auxiliary, orientations between the
two reference frames, this is sufficient for obtaining the desired direction
cosine matrix and the matrix multiplication is handled internally for you:

.. jupyter-execute::

   sm.trigsimp(B.dcm(N))

Lastly, recall the general definition of the direction cosine matrix. We showed
that the dot product of pairs of unit vectors give the entries to the direction
cosine matrix. ``mechanics`` has a
:external:py:func:`~sympy.physics.vector.functions.dot` function that can
calculate the dot product of two vectors. Using it on two of the unit vector
pairs returns the expected direction cosine matrix entry:

.. jupyter-execute::

   sm.trigsimp(me.dot(B.x, N.x))

.. admonition:: Exercise

   Orient reference frame :math:`D` with respect to :math:`C` with a simple
   rotation through angle :math:`\beta` about the shared :math:`-y` axis.  Use
   the direction cosine matrix from this first orientation to set the
   orientation of reference frame :math:`E` with respect to :math:`D`. Show
   that both pairs of reference frames have the same relative orientations.

.. admonition:: Solution
   :class: dropdown

   .. jupyter-execute::

      beta = sm.symbols('beta')

      C = me.ReferenceFrame('C')
      D = me.ReferenceFrame('D')
      E = me.ReferenceFrame('E')

      D.orient_axis(C, beta, -C.y)

      D.dcm(C)

   .. jupyter-execute::

      E.orient_explicit(D, C.dcm(D))
      E.dcm(D)

Euler Angles
============

The camera stabilization gimbal_ shown in :numref:`camera-gimbal`  has three
`revolute joints`_ that orient the camera :math:`D` relative to the handgrip
frame :math:`A`.

.. _camera-gimbal:

.. figure:: https://moorepants.info/mechmotum-bucket/orientation-camera-gimbal.png
   :align: center

   Four reference frames labeled on the Turnigy Pro Steady Hand Camera Gimbal.
   *Image copyright HobbyKing, used under fair use for educational purposes.*

If we introduce two additional auxiliary reference frames, :math:`B` and
:math:`C`, attached to the intermediate camera frame members, we can use three
successive simple orientations to go from :math:`A` to :math:`D`. We can
formulate the direction cosine matrices for the reference frames using the same
technique for the successive simple orientations shown in :ref:`Successive
Orientations`, but now our sequence of three orientations will enable us to
orient :math:`D` in any way possible relative to :math:`A` in three dimensional
space.

.. _gimbal: https://en.wikipedia.org/wiki/Gimbal
.. _revolute joints: https://en.wikipedia.org/wiki/Revolute_joint

Watch this video to get a sense of the orientation axes for each intermediate
auxiliary reference frame:

.. raw:: html

   <center>
      <iframe
        width="560"
        height="315"
        src="https://www.youtube.com/embed/xQMBIXqWcjI?start=177"
        title="YouTube video player"
        frameborder="0"
        allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture"
        allowfullscreen>
      </iframe>
   </center>

We first orient :math:`B` with respect to :math:`A` about the shared :math:`z`
unit vector through the angle :math:`\psi`, as shown below:

.. _orientation-gimbal-psi:

.. figure:: figures/orientation-gimbal-psi.svg
   :width: 200px
   :align: center

   View of the :math:`A` and :math:`B` :math:`x\textrm{-}y` plane showing the
   orientation of :math:`B` relative to :math:`A` about :math:`z` through angle
   :math:`\psi`.

In SymPy, use :external:py:class:`~sympy.physics.vector.frame.ReferenceFrame`
to establish the relative orientation:

.. jupyter-execute::

   psi = sm.symbols('psi')

   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')

   B.orient_axis(A, psi, A.z)

   B.dcm(A)

Now orient :math:`C` with respect to :math:`B` about their shared :math:`x`
unit vector through angle :math:`\theta`.

.. _orientation-gimbal-theta:

.. figure:: figures/orientation-gimbal-theta.svg
   :width: 200px
   :align: center

   View of the :math:`B` and :math:`C` :math:`y\textrm{-}z` plane showing the
   orientation of :math:`C` relative to :math:`B` about :math:`x` through angle
   :math:`\theta`.

.. jupyter-execute::

   theta = sm.symbols('theta')

   C = me.ReferenceFrame('C')

   C.orient_axis(B, theta, B.x)

   C.dcm(B)

Finally, orient the camera :math:`D` with respect to :math:`C` about their
shared :math:`y` unit vector through the angle :math:`\phi`.

.. figure:: figures/orientation-gimbal-phi.svg
   :width: 200px
   :align: center

   View of the :math:`C` and :math:`D` :math:`x\textrm{-}z` plane showing the
   orientation of :math:`D` relative to :math:`C` about :math:`y` through angle
   :math:`\varphi`.

.. jupyter-execute::

   phi = sm.symbols('varphi')

   D = me.ReferenceFrame('D')

   D.orient_axis(C, phi, C.y)

   D.dcm(C)

With all of the intermediate orientations defined, when can now ask for the
relationship :math:`{}^D\mathbf{C}^A` of the camera :math:`D` relative to the
handgrip frame :math:`A`:

.. jupyter-execute::

   D.dcm(A)

With these three successive orientations the camera can be rotated arbitrarily
relative to the handgrip frame. These successive
:math:`z\textrm{-}x\textrm{-}y` orientations are a standard way of describing
the orientation of two reference frames and are referred to as `Euler Angles`_
[#]_.

.. _Euler Angles: https://en.wikipedia.org/wiki/Euler_angles

There are 12 valid sets of successive orientations that can arbitrarily orient
one reference frame with respect to another. These are the six "Proper Euler
Angles":

.. math::

   z\textrm{-}x\textrm{-}z, x\textrm{-}y\textrm{-}x, y\textrm{-}z\textrm{-}y,
   z\textrm{-}y\textrm{-}z, x\textrm{-}z\textrm{-}x, y\textrm{-}x\textrm{-}y

and the six "Tait-Bryan Angles":

.. math::

   x\textrm{-}y\textrm{-}z, y\textrm{-}z\textrm{-}x, z\textrm{-}x\textrm{-}y,
   x\textrm{-}z\textrm{-}y, z\textrm{-}y\textrm{-}x, y\textrm{-}x\textrm{-}z

Different sets can be more or less suitable for the kinematic nature of the
system you are describing. We will also refer to these 12 possible orientation
sets as "body fixed orientations". As we will soon see, a rigid body and a
reference frame are synonymous from an orientation perspective and each
successive orientation rotates about a shared unit vector fixed in both of the
reference frames (or bodies), thus "body fixed orientations". The method
:external:py:meth:`~sympy.physics.vector.frame.ReferenceFrame.orient_body_fixed`
can be used to establish the relationship between :math:`A` and :math:`D`
without the need to create auxiliary reference frames :math:`B` and :math:`C`:

.. jupyter-execute::

   A = me.ReferenceFrame('A')
   D = me.ReferenceFrame('D')

   D.orient_body_fixed(A, (psi, theta, phi), 'zxy')

   D.dcm(A)

.. todo:: The wikipedia animation is not correct. The lower yellow arrow should
   be colored green. This needs to be replaced with a corrected animation.

.. admonition:: Exercise

   Euler_ discovered 6 of the 12 orientation sets. One of these sets is shown
   in this figure:

   .. _orientation-euler-animation:

   .. figure:: https://objects-us-east-1.dream.io/mechmotum/euler-angle-animation.gif
      :align: center

      An orientation through Euler angles with frame :math:`A` (yellow),
      :math:`B` (red), :math:`C` (green), and :math:`D` (blue). The rightward
      yellow arrow is the :math:`x` direction, leftward yellow arrow is the
      :math:`y` direction, and upward yellow arrow is the :math:`z` direction.
      All frames' unit vectors are aligned before being oriented.

   Take the acute angles between :math:`A` and :math:`B` to be :math:`\psi`,
   :math:`B` and :math:`C` to be :math:`\theta`, and :math:`C` and :math:`D` to
   be :math:`\varphi`. Determine what Euler angle set this is and then
   calculate :math:`{}^D\mathbf{C}^A` using
   :external:py:meth:`~sympy.physics.vector.frame.ReferenceFrame.orient_axis`
   and then with
   :external:py:meth:`~sympy.physics.vector.frame.ReferenceFrame.orient_body_fixed`
   showing that you get the same result.

   .. _Euler: https://en.wikipedia.org/wiki/Leonhard_Euler

.. admonition:: Solution
   :class: dropdown

   The Euler angle set is :math:`z\textrm{-}x\textrm{-}z`.

   .. jupyter-execute::

      psi, theta, phi = sm.symbols('psi, theta, varphi')

   With :external:py:meth:`~sympy.physics.vector.frame.ReferenceFrame.orient_axis`:

   .. jupyter-execute::

      A = me.ReferenceFrame('A')
      B = me.ReferenceFrame('B')
      C = me.ReferenceFrame('C')
      D = me.ReferenceFrame('D')

      B.orient_axis(A, psi, A.z)
      C.orient_axis(B, theta, B.x)
      D.orient_axis(C, phi, C.z)

      D.dcm(A)

   With :external:py:meth:`~sympy.physics.vector.frame.ReferenceFrame.orient_body_fixed`:

   .. jupyter-execute::

      A = me.ReferenceFrame('A')
      D = me.ReferenceFrame('D')

      D.orient_body_fixed(A, (psi, theta, phi), 'zxz')

      D.dcm(A)

Alternatives for Representing Orientation
==========================================

In the previous section, Euler-angles were used to encode the orientation of a
frame or body.
There are `many alternative approaches to representing orientations
<https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions>`_.
Three such representations, which will be used throughout this book, were
already introduced:

* **Euler-angles** themselves, which provides a minimal representation (only 3
  numbers), and a relatively straightforward way to compute the change in
  orientation from the angular velocity (see :ref:`Angular Kinematics`).
* the **direction cosine matrix**, which allow easy rotations or vectors and
  consecutive rotations, both via matrix multiplication,
* the **axis-angle representation** (used in the
  :external:py:meth:`~sympy.physics.vector.frame.ReferenceFrame.orient_axis`
  method), which is often an intuitive way to describe the orientation for
  manual input, and is useful when the axis of rotation is fixed.

Each representation also has downsides.
For example, the direction cosine matrix consists of nine elements; more to
keep track of than three Euler angles.
Furthermore, not all combinations of nine elements form a valid direction
cosine matrix, so we have to be careful to check and enforce validity when
writing code.

Learn more
==========

One more frequently used approach to representing orientations is based on so
called `quaternions`_.
Quaternions are like imaginary numbers, but with three imaginary constants:
:math:`i`, :math:`j` and :math:`k`.
These act as described by the rule

.. math::

   i^2 = j^2 = k^2 = ijk = -1.

A general quaternion :math:`q` can thus be written in terms of its components
:math:`q_0`, :math:`q_i` :math:`q_j`, :math:`q_k` which are real numbers:

.. math::

   q = q_0 + q_ii + q_jj + q_kk

.. _quaternions: https://en.wikipedia.org/wiki/Quaternion

The
:external:py:meth:`~sympy.physics.vector.frame.ReferenceFrame.orient_quaternion`
method enables orienting a reference frame using a quaternion in sympy:

.. jupyter-execute::

   N = me.ReferenceFrame('N')
   A = me.ReferenceFrame('A')

   q_0, qi, qj, qk = sm.symbols('q_0 q_i q_j q_k')
   q = (q_0, qi, qj, qk)
   A.orient_quaternion(N, q)
   A.dcm(N)

A rotation of an angle :math:`\theta` around a unit vector :math:`\hat{e}` can
be converted to a quaternion representation by having :math:`q_0 =
\cos\left(\frac{\theta}{2}\right)`, and the other components equal to a factor
:math:`\sin\left(\frac{\theta}{2}\right)` times the components of the axis of
rotation :math:`\hat{e}`.
For example, if the rotation axis is :math:`\hat{n}_x`, we get:

.. jupyter-execute::

   q = (sm.cos(theta/2), sm.sin(theta/2), 0, 0)
   A.orient_quaternion(N, q)
   sm.trigsimp(A.dcm(N))

The length of a quaternion is the square root of the sum of the squares of its
components.
For a quaternion representing an orientation, this length must always be 1.

It turns out that the multiplication rules for (unit) quaternions provide an
efficient way to compose multiple rotations, and to numerically integrate the
orientation when given an angular velocity.
Due to the interpretation related to the angle and axis representation, it is
also a somewhat intuitive representation.
However, the integration algorithm needs to take an additional step to ensure
the quaternion always has unit length.

The representation of orientations in general, turns out to be related to an
area of mathematics called Lie-groups.
The theory of Lie-groups has further applications to the mechanics and control
of multibody systems.
An example application is finding a general method for simplifying the
equations for symmetric systems, so this can be done more easily and to more
systems. The Lie-group theory is not used in this book.
Instead, the interested reader can look up the `3D rotation group`_ as a
starting point for further study.

.. _3D rotation group: https://en.wikipedia.org/wiki/3D_rotation_group

.. rubric:: Footnotes
.. [#] Technically, this set of angles for the gimbal are one of the 6 Tait-Bryan angles,
   but “Euler Angles” is used as a general term to describe both Tait-Bryan angles
   and “proper Euler angles”.
