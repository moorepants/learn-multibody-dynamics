=======
Vectors
=======

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`vectors` or Jupyter Notebook:
   :jupyter-download:notebook:`vectors`.

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   sm.init_printing(use_latex='mathjax')

What is a vector?
=================

Vectors have three characteristics:

1. magnitude
2. orientation
3. sense

The direction the vector points is derived from both the orientation and the
sense. Vectors are equal when all three characteristics are the same.

.. note::

   In this text we will distinguish scalar variables, e.g. :math:`v`, from
   vectors by including a bar over the top of the symbol, e.g. :math:`\bar{v}`.

Vectors have these mathematical properites:

- scalar mutiplicative: :math:`\lambda\bar{a}` where :math:`\lambda` can only
  change the magnitude and the sense of the vector
- commutative: :math:`\bar{a} + \bar{b} = \bar{b} + \bar{a}`
- distributive: :math:`\lambda(\bar{a} + \bar{b}) = \lambda\bar{a} +
  \lambda\bar{b}`
- assocative: :math:`(\bar{a} + \bar{b}) + \bar{c} = \bar{a} + (\bar{b} +
  \bar{c})`
- common orientation: :math:`\bar{b} = k\bar{a}` where :math:`\bar{b}` and
  :math:`\bar{b}` have the same orientation

Unit vectors are vectors with a magnitude of :math:`1`. If the magnitude of
:math:`\bar{v}` is 1, then we indicate this with :math:`\hat{v}`. Any vector
has an assocated unit vector with the same orientation and sense, found by:

.. math::

   \hat{u} = \frac{\bar{u}}{||\bar{u}||}

where :math:`||\bar{u}||` is the `Euclidean norm`_ (2-norm), or magnitude, of
the vector :math:`\bar{u}`.

.. _Euclidean norm: https://en.wikipedia.org/wiki/Norm_(mathematics)#Euclidean_norm

Vector Functions
================

Vectors can be functions of scalar variables. :math:`\bar{v}` is a vector
function of scalar variable :math:`q` in reference frame :math:`A` if
:math:`\bar{v}` changes when viewed in :math:`A` when :math:`q` changes. Note
that this implies that :math:`\bar{v}` may not be a function of scalar variable
:math:`q` in a different reference frame.

Staring with reference frame :math:`A` and vector :math:`\bar{v}` which is a
function of :math:`n` scalars :math:`q_1,q_2,\ldots,q_n` in :math:`A`, let
:math:`\hat{a}_x,\hat{a}_y,\hat{a}_z` be a set of mutually perpendicular unit
vectors fixed in :math:`A`, i.e. :math:`\hat{a}_x,\hat{a}_y,\hat{a}_z` are
constant when observed from :math:`A`. Then there are three unique scalar
functions :math:`v_x,v_y,v_z` of :math:`q_1,q_2,\ldots,q_n` such that:

.. math::

   \bar{v} = v_x \hat{a}_x + v_y \hat{a}_y + v_y \hat{a}_y

:math:`v_x \hat{a}_x` is called the :math:`\hat{a}_x` component of
:math:`\bar{v}` and :math:`v_x` is called measure number of :math:`\bar{v}`.
Since the components are mutually perpendicular the measure number can also be
found from the dot product of :math:`\bar{v}` and the respective unit vector:

.. math::

   \bar{v} = (\bar{v} \cdot \hat{a}_x) \hat{a}_x +
             (\bar{v} \cdot \hat{a}_y) \hat{a}_y +
             (\bar{v} \cdot \hat{a}_z) \hat{a}_z

Vector with SymPy Mechanics
===========================

Addition
--------

When we add vector :math:`\bar{b}` to vector :math:`\bar{a}`, the result is
a vector that starts at the tail of :math:`\bar{a}` and ends at the tip of
:math:`\bar{b}`:

.. figure:: vector_addition.svg
   :alt: Vector addition
   :align: center

   Vector addition

Vectors in SymPy Mechanics are created by first introducing a reference frame
and using its associated unit vectors to construct vectors of arbritrary
magnitude and direction.

.. jupyter-execute::

   N = me.ReferenceFrame('N')

Now introdcue some scalar variables:

.. jupyter-execute::

   a, b, c, d, e, f = sm.symbols('a, b, c, d, e, f')

The simplest 2D non-unit vector is made up of a single component:

.. jupyter-execute::

   v = a*N.x
   v

A, possible more familiar, column matrix form of a vector is accessed with the
:external:py:meth:`~sympy.physics.vector.vector.Vector.to_matrix`.

.. jupyter-execute::

   v.to_matrix(N)

Fully 3D and arbitray vectors can be created by providing a measure number for
each unit vector of :math:`N`:

.. jupyter-execute::

   w = a*N.x + b*N.y + c*N.z
   w

And the associated column matrix form:

.. jupyter-execute::

   w.to_matrix(N)

Vector addition works by adding the measure numbers of each common component:

.. math::

   \bar{w} = & a \hat{n}_x + b \hat{n}_y + c \hat{n}_z \\
   \bar{x} = & d \hat{n}_x + e \hat{n}_y + f \hat{n}_z \\
   \bar{w} + \bar{x} = & (a + d) \hat{n}_x + (b + e) \hat{n}_y + (c + f) \hat{n}_z

SymPy Mechanics vectors work as expected:

.. jupyter-execute::

   x = d*N.x + e*N.y + f*N.z
   x

.. jupyter-execute::

   w + x

Scaling
-------

Multiplying a vector by a scalar changes its magnitude, but not its
orientation. Scaling by a negative number changes a vector's magnitude and
reverses its sense (rotates it by :math:`\pi` radians).

.. figure:: vector_scaling.svg
   :alt: Vector scaling

   Vector scaling

.. jupyter-execute::

   y = 2*w
   y

.. jupyter-execute::

   z = -w
   z

.. admonition:: Exercise

   Create three vectors that lie in the :math:`xy` plane of reference frame
   :math:`N` where each vector is:

   1. of length :math:`l` that is at an angle of :math:`\frac{\pi}{4}`
      degrees from the :math:`\hat{n}_x` unit vector.
   2. of length :math:`10` and is in the :math:`-\hat{n}_y` direction
   3. of length :math:`l` and is :math:`\theta` radians from the
      :math:`\hat{n}_y` unit vector.

   Finally, add vectors from 1 and 2 and substract :math:`5` times the third
   vector.

   Hint: SymPy has fundamental constants and trigonometic functions, for
   example ``sm.tan, sm.pi``.

.. admonition:: Solution
   :class: dropdown

   .. jupyter-execute::

      N = me.ReferenceFrame('N')
      l, theta = sm.symbols('l, theta')

   .. jupyter-execute::

      v1 = l*sm.cos(sm.pi/4)*N.x + l*sm.sin(sm.pi/4)*N.y
      v1

   .. jupyter-execute::

      v2 = -10*N.y
      v2

   .. jupyter-execute::

      v3 = -l*sm.sin(theta)*N.x + l*sm.cos(theta)*N.y
      v3

   .. jupyter-execute::

      v1 + v2 - 5*v3

Dot Product
-----------

The dot product, which yields a scalar quantity, is defined as:

.. math::

   \bar{v} = & v_x \hat{n}_x + v_y \hat{n}_y + v_z \hat{n}_z \\
   \bar{w} = & w_x \hat{n}_x + w_y \hat{n}_y + w_z \hat{n}_z \\
   \bar{v} \cdot \bar{w} = & v_x w_x + v_v w_y + v_z w_z

and is also equivlant to:

.. math::

   \bar{v} \cdot \bar{w} = ||\bar{v}|| ||\bar{w}|| \cos{\theta}

where :math:`\theta` is the angle between the two vectors.

.. figure:: vector_dot.svg
   :alt: Vector dot product

   Vector dot product

The dot product is often used to determine:

-  the angle between two vectors:
   :math:`\theta = \arccos\frac{\bar{a} \cdot \bar{b}}{|\bar{a}||\bar{b}|}`

-  a vector’s magnitude: :math:`||\bar{v}|| = \sqrt{\bar{v} \cdot \bar{v}}`

-  the length of a vector along a direction of another vector :math:`\hat{u}`
   (called the projection):
   :math:`\mbox{proj}_\hat{u} \bar{v} = \bar{v} \cdot \hat{u}`

-  if two vectors are perpendicular: :math:`\bar{v} \cdot \bar{w} = 0 \mbox{ if }\bar{v} \perp \bar{w}`

Also, dot products are used to convert a vector equation into a scalar equation
by "dotting" an entire equation with a vector.

.. jupyter-execute::

    N = me.ReferenceFrame('N')
    w = a*N.x + b*N.y + c*N.z
    x = d*N.x + e*N.y + f*N.z

The :external:py:func:`~sympy.physics.vector.functions.dot` function
calculates the dot product:

.. jupyter-execute::

    me.dot(w, x)

The :external:py:meth:`~sympy.physics.vector.vector.Vector.normalize`

You can compute a unit vector in the same direction as :math:`\bar{w}` like so:

.. jupyter-execute::

   w.normalize()

.. admonition:: Exercise

   Write your own function that normalizes an arbitrary vector and show that it
   gives the same result as ``w.normalize()``.

.. admonition:: Solution
   :class: dropdown

   The do

   .. jupyter-execute::

      def normalize(vector):
          return vector/sm.sqrt(me.dot(vector, vector))

      normalize(w)

SymPy Mechanics vectors also have a method
:external:py:meth:`~sympy.physics.vector.vector.Vector.magnitude` which is
helpful:

.. jupyter-execute::

   w.magnitude()

.. jupyter-execute::

   w/w.magnitude()

.. admonition:: Exercise

   Given the vectors
   :math:`\bar{v}_1 = a \hat{\mathbf{n}}_x + b\hat{\mathbf{n}}_y + a \hat{\mathbf{n}}_z`
   and
   :math:`\bar{v}_2=b \hat{\mathbf{n}}_x + a\hat{\mathbf{n}}_y + b \hat{\mathbf{n}}_z`
   find the angle between the two vectors using the dot product.

.. admonition:: Solution
   :class: dropdown

   .. jupyter-execute::

      N = me.ReferenceFrame('N')
      v1 = a * N.x + b * N.y + a * N.z
      v2 = b * N.x + a * N.y + b * N.z

   .. jupyter-execute::

      sm.acos(v1.dot(v2) / (v1.magnitude()*v2.magnitude()))


Cross Product
-------------

The `cross product`_, which yields a vector quantity, is defined as:

.. math::  \bar{v} \times \bar{w} = |\bar{v}||\bar{w}| \sin\theta \hat{u}

where :math:`\theta` is the angle between the two vectors, and :math:`\hat{u}`
is the unit vector perpendicular to both :math:`\bar{v}` and :math:`\bar{w}`
whose sense is given by the right-hand rule. For aribtrary measure numbers this
results in the following:

.. math::

   \bar{v} = & v_x \hat{n}_x + v_y \hat{n}_y + v_z \hat{n}_z \\
   \bar{w} = & w_x \hat{n}_x + w_y \hat{n}_y + w_z \hat{n}_z \\
   \bar{v} \times \bar{w} = &
   v_y w_z - v_z w_y  \hat{n}_x +
   v_z w_x - v_x w_z \hat{n}_y +
   v_x w_y - v_y w_x \hat{n}_z

.. _cross product: https://en.wikipedia.org/wiki/Cross_product

The cross product is used to:

-  obtain a vector/direction perpendicular to two other vectors
-  determine if two vectors are parallel:
   :math:`\bar{a} \times \bar{b} = \bar{0} \mbox{ if } \bar{a} \parallel \bar{b}`
-  compute moments: :math:`\bar{r} \times \bar{F}`
-  compute the area of a triangle

.. figure:: vector_cross.svg
   :alt: Vector cross product

   Vector cross product

.. jupyter-execute::

    N = me.ReferenceFrame('N')
    w = a*N.x + b*N.y + c*N.z
    x = d*N.x + e*N.y + f*N.z

.. jupyter-execute::

    me.cross(w, x)

.. note:: Exercise

   Given three points located in reference frame :math:`N` by:

   .. math::

      \bar{p}_1 = 23 \hat{n}_x - 12 \hat{n}_y \\
      \bar{p}_2 = 16 \hat{n}_x + 2 \hat{n}_y - 4 \hat{n}_z \\
      \bar{p}_3 = \hat{n}_x + 14 \hat{n}_z

Find the area of the triangle bounded by these three points using the
cross product.

.. note::

   Hint: Search online for the relationship of the cross product to triangle
   area.

Some vector properties
----------------------

- The order in which you add them does not matter:
  :math:`\bar{a} + \bar{b} = \bar{b} + \bar{a}`
- You can distrubute a scalar among vectors:
  :math:`s (\bar{a} + \bar{b}) = s\bar{a} + s\bar{b}`

**Dot product**

-  You can pull out scalars: $ c :math:`\bar{a}` :math:`\times `d
   :math:`\bar{b}` = cd (:math:`\bar{a}`
   :math:`\times `:math:`\bar{b}`)$
-  Order does not matter: :math:`\bar{a} \cdot \bar{b} = \bar{b} \cdot \bar{a}`
-  You can distribute:
   :math:`\bar{a} \cdot (\bar{b} + \bar{c}) = \bar{a} \cdot \bar{b} + \bar{a} \cdot \bar{c}`

**Cross product**

-  Crossing a vector with itself “cancels” it:
   :math:`\bar{a} \times \bar{b} = \vec{0}`
-  You can pull out scalars: $ c :math:`\bar{a}`
   :math:`\times `d :math:`\bar{b}` = cd (:math:`\bar{a}`
   :math:`\times `:math:`\bar{b}`)$
-  Order DOES matter (because of the right-hand rule):
   :math:`\bar{a} \times \bar{b} = -\bar{b} \times \bar{a}`
-  You can distribute:
   :math:`\bar{a} \times (\bar{b} + \bar{c}) = \bar{a} \times \bar{b} + \bar{a} \times \bar{c}`
-  They are NOT associative:
   :math:`\bar{a} \times ({\bar{b} \times \bar{c}) \neq {(\bar{a} \times \bar{b}) \times \bar{c}`
