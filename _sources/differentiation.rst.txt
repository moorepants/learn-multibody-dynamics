======================
Vector Differentiation
======================

.. note::

   You can download this example as a Python script:
   :jupyter-download-script:`differentiation` or Jupyter Notebook:
   :jupyter-download-notebook:`differentiation`.

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   sm.init_printing(use_latex='mathjax')

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

Learning Objectives
===================

After completing this chapter readers will be able to:

- Calculate the partial derivative of a vector with respect to any variable
  when viewed from any reference frame.
- Use the product rule to find the relationship of changing measure numbers and
  changing unit vectors.
- Explain the difference in expressing a vector in a reference frame and taking
  the derivative of the vector when observed from the reference frame.
- Calculate second partial derivatives.
- Calculate time derivatives of vector functions.

Partial Derivatives
===================

If a vector :math:`\bar{v}` is a function of :math:`n` scalar variables
:math:`q_1,q_2,\ldots,q_n` in reference frame :math:`A` then the first partial
derivatives of :math:`\bar{v}` in :math:`A` with respect to :math:`q_r` where
:math:`r=1\ldots n` can be formed by applying the product rule of
differentation and taking into account that the mutually perpendicular unit
vectors fixed in :math:`A` do not change when observed from :math:`A`. The
partial derivatives are then:

.. math::
   :label: partial-deriv-def

   \frac{{}^A\partial \bar{v}}{\partial q_r} = \sum_{i=1}^3 \frac{\partial
   v_i}{\partial q_r} \hat{a}_i \textrm{ for } r=1\ldots n

where :math:`v_i` are the measure numbers of :math:`\bar{v}` expressed in
:math:`A` associated with the mutually perpendicular unit vectors
:math:`\hat{a}_1,\hat{a}_2,\hat{a}_3`.

If :math:`\bar{v}=v_x\hat{a}_x+v_y\hat{a}_y+v_z\hat{a}_z` the above definition
expands to:

.. math::
   :label: partial-measure

   \frac{{}^A\partial \bar{v}}{\partial q_r} =
   \frac{\partial v_x}{\partial q_r} \hat{a}_x +
   \frac{\partial v_y}{\partial q_r} \hat{a}_y +
   \frac{\partial v_z}{\partial q_r} \hat{a}_z
   \textrm{ for } r=1\ldots n

Many of the vectors we will work with in multibody dynamics will be a function
of a single variable, most often time :math:`t`. If that is the case, the
partial derivative reduces to a single variate derivative:

.. math::
   :label: single-var-deriv-def

   \frac{{}^A d \bar{v}}{dt} := \sum_{i=1}^3 \frac{d v_i}{dt} \hat{a}_i

.. warning::

   A derivative written as :math:`\frac{\partial \bar{v}}{\partial q_r}` is
   meaningless because no reference frame is indicated. The derivative is
   dependent on which reference frame the change is observed from, so without a
   reference frame, the derivative cannot be calculated. This is not the case
   for partial derivatives of scalar expressions, as no reference frame is
   involved.

The above definition implies that a vector must be expressed in the reference
frame one is observing the change from before calculating the partial
derivatives of the scalar measure numbers. For example, here is a vector that
is expressed with unit vectors from three different reference frames:

.. jupyter-execute::

   alpha, beta = sm.symbols('alpha, beta')
   a, b, c, d, e, f = sm.symbols('a, b, c, d, e, f')

   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')
   C = me.ReferenceFrame('C')

   B.orient_axis(A, alpha, A.x)
   C.orient_axis(B, beta, B.y)

   v = a*A.x + b*A.y + c*B.x + d*B.y + e*C.x + f*C.y
   v

To calculate :math:`\frac{{}^A\partial\bar{v}}{\partial \alpha}` we first need
to project the vector :math:`\bar{v}` onto the unit vectors of :math:`A` and
take the partial derivative of those measure numbers with respect to
:math:`\alpha`. The dot product provides the projection and the resulting
scalar is differentiated:

.. jupyter-execute::

   dvdalphaAx = v.dot(A.x).diff(alpha)
   dvdalphaAx

.. jupyter-execute::

   dvdalphaAy = v.dot(A.y).diff(alpha)
   dvdalphaAy

.. jupyter-execute::

   dvdalphaAz = v.dot(A.z).diff(alpha)
   dvdalphaAz

We can then construct the vector :math:`\frac{{}^A\partial \bar{v}}{\partial
\alpha}` from the new measure numbers know that the :math:`A` unit vectors are
fixed:

.. jupyter-execute::

   dvdalphaA = dvdalphaAx*A.x + dvdalphaAy*A.y + dvdalphaAz*A.z
   dvdalphaA

SymPy Mechanics vectors have a special
:external:py:meth:`~sympy.physics.vector.vector.Vector.diff` method that
manages taking partial derivatives from different reference frames. For the
vector ``.diff()`` method you provide first the variable :math:`\alpha`
followed by the reference frame you are observing from:

.. jupyter-execute::

   v.diff(alpha, A)

This gives the identical result as our manually constructed partial derivative
above.

.. admonition:: Exercise

   Calculate :math:`\frac{{}^B\partial \bar{v}}{\partial e}` manually and with
   :external:py:meth:`~sympy.physics.vector.vector.Vector.diff` and show the
   results are the same.

.. admonition:: Solution
   :class: dropdown

   .. jupyter-execute::

      dvdeBx = v.dot(B.x).diff(e)
      dvdeBy = v.dot(B.y).diff(e)
      dvdeBz = v.dot(B.z).diff(e)
      dvdeBx*B.x + dvdeBy*B.y + dvdeBz*B.z

   .. jupyter-execute::

      v.diff(e, B).express(B)

.. warning:: What's the difference in `.express()` and `.diff()`?

   Any vector can be "expressed" in any reference frame. To express a vector in
   a reference frame means to project it onto the three mutually perpendicular
   unit vectors fixed in the reference frame and then to rewrite the vector in
   terms of measure numbers associated with those three unit vectors using the
   relevant direction cosine matrix entries. This has nothing to do with
   differentiation.

   We can also take the derivative of a vector when viewed from a specific
   reference frame. To do so, we observe how the vector changes when viewed
   from the reference frame and formulate that derivative. Once the derivative
   is taken, we can express the new vector in any reference frame we desire.

   Expressing a vector in a reference frame and taking a derivative of a vector
   when observered from a reference frame are two different things! Try not to
   get tripped up by this important distinction.

Product Rule
============

Consider again vector :math:`\bar{v}=v_x\hat{a}_x+v_y\hat{a}_y+v_z\hat{a}_z`.
Previously, only the measure numbers of this vector were scalar functions of
:math:`q_r`. Now consider a reference frame :math:`N` that is oriented relative
to :math:`A` such that the relative orientation also depends on :math:`q_r`.
This means, that when observed from :math:`N`, the unit vectors
:math:`\hat{a}_x,\hat{a}_y,\hat{a}_z` may be a function of :math:`q_r`. With
both the measure numbers and unit vectors dependent on :math:`q_r` the
derivative of :math:`\bar{v}` in :math:`N` requires the use of the product rule
when taking the partial derivative. For example:

.. math::
   :label: product-rule-big

   \frac{{}^N\partial \bar{v}}{\partial q_r} =
   \frac{{}^N\partial v_x}{\partial q_r}\hat{a}_x + v_x \frac{{}^N\partial \hat{a}_x}{\partial q_r} +
   \frac{{}^N\partial v_y}{\partial q_r}\hat{a}_y + v_y \frac{{}^N\partial \hat{a}_y}{\partial q_r} +
   \frac{{}^N\partial v_z}{\partial q_r}\hat{a}_z + v_z \frac{{}^N\partial \hat{a}_z}{\partial q_r}

The three similar terms with scalar derivatives have the same interpretation of
the ones in the prior section.

.. math::
   :label: product-rule-part-01

   \frac{{}^N\partial v_x}{\partial q_r}\hat{a}_x,
   \frac{{}^N\partial v_y}{\partial q_r}\hat{a}_y,
   \frac{{}^N\partial v_z}{\partial q_r}\hat{a}_z

But the part with unit vector derivatives is more interesting. The partial
derivative of a unit vector depends on how it changes. But unit vectors do not
change in length, only in orientation.

.. math::
   :label: product-rule-part-02

   v_x \frac{{}^N\partial \hat{a}_x}{\partial q_r},
   v_y \frac{{}^N\partial \hat{a}_y}{\partial q_r},
   v_z \frac{{}^N\partial \hat{a}_z}{\partial q_r}

You will learn in the next chapter how to interpret and use these terms to
simplify the calculations of common derivatives. But for now, just be aware of
the nature of this partial derivative in :math:`N`.

The product rule also applies to the dot and cross products:

.. math::
   :label: product-dot-cross

   \frac{\partial}{\partial q_r}(\bar{v} \cdot \bar{w}) = &
   \frac{\partial \bar{v}}{\partial q_r} \cdot \bar{w} +
   \bar{v} \cdot \frac{\partial \bar{w}}{\partial q_r}

   \frac{\partial}{\partial q_r}(\bar{v} \times \bar{w}) = &
   \frac{\partial \bar{v}}{\partial q_r} \times \bar{w} +
   \bar{v} \times \frac{\partial \bar{w}}{\partial q_r}

and generalizes to any series of products. Let :math:`G=f_1 \cdots f_n` be a
series of products, then:

.. math::
   :label: product-rule-gen

   \frac{\partial G}{\partial q_r} =
   \frac{\partial f_1}{\partial q_r}\cdot f_2 \cdots f_n +
   f_1 \cdot\frac{\partial f_2}{\partial q_r}\cdot f_3 \cdots f_n +
   \dots +
   f_1 \cdots f_{n-1} \cdot \frac{\partial f_n}{\partial q_r}

Second Derivatives
==================

:math:`\frac{{}^A\partial \bar{v}}{\partial q_r}` is also a vector and, just
like :math:`\bar{v}`, may be a vector function. We can thus calculate the
second partial derivative with respect to :math:`q_s` where :math:`s=1\ldots
n`. This second partial derivative need not be taken with respect to the same
reference frame as the first partial derivative. If we first differentiate with
when viewed from :math:`A` and then when viewed from :math:`B`, the second
partial derivative is:

.. math::
   :label: second-derivative

   \frac{{}^B\partial}{\partial q_s} \left(\frac{{}^A\partial\bar{v}}{\partial
   q_r}\right)

Second partials in different reference frames do not necessarily commute:

.. math::
   :label: no-commute-second-deriv

   \frac{{}^B\partial}{\partial q_s} \left(\frac{{}^A\partial\bar{v}}{\partial
   q_r}\right)
   \neq
   \frac{{}^A\partial}{\partial q_r} \left(\frac{{}^B\partial\bar{v}}{\partial
   q_s}\right)

If the reference frames of each partial derivative are the same, then `mixed
partials do commute`_.

.. _mixed partials do commute: https://en.wikipedia.org/wiki/Symmetry_of_second_derivatives

.. todo:: Make an example of second derivatives not commuting.

Vector Functions of Time
========================

In multibody dynamics we are primarily concerned with how motion changes with
respect to time :math:`t` and our vectors and measure numbers will often be
implicit functions of time, i.e. :math:`q_r(t)`. When that is the case the
chain rule can be used to take total derivatives:

.. math::
   :label: time-deriv

   \frac{{}^A d\bar{v}}{dt} =
   \sum_{i=1}^n \frac{{}^A\partial \bar{v}}{\partial q_r(t)} \frac{d q_r(t)}{dt} +
   \frac{{}^A \partial \bar{v}}{\partial t}
   \textrm{ where } r=1,\ldots,n

.. note::

   We will typically use the "dot" notation for time derivatives, i.e.
   :math:`\frac{dq}{dt}` as :math:`\dot{q}` and :math:`\frac{d^2q}{dt^2}` as
   :math:`\ddot{q}` and so on.

In SymPy Mechanics, scalar functions of time can be created like so:

.. jupyter-execute::

   t = sm.symbols('t')
   q_of = sm.Function('q')

   q = q_of(t)
   q

And these scalar functions can be differentiated:

.. jupyter-execute::

   q.diff(t)

SymPy Mechanics provides the convenience function
:external:py:func:`~sympy.physics.vector.dynamicsymbols` to create scalar
functions of time just like ``symbols()``:

.. jupyter-execute::

   q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')
   q1, q2, q3

The time variable used in ``q1,q2,q3`` can be accessed like so:

.. jupyter-execute::

   t = me.dynamicsymbols._t

SymPy Mechanics also provide a special printing function
:external:py:func:`~sympy.physics.vector.printing.init_vprinting` which enables
the dot notation on functions of time:

.. jupyter-execute::

   me.init_vprinting(use_latex='mathjax')
   q1.diff(t), q2.diff(t, 2), q3.diff(t, 3)

Now these scalar functions of time can be used to formulate vectors:

.. jupyter-execute::

   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')
   B.orient_body_fixed(A, (q1, q2, q3), 'ZXZ')
   v = q1*A.x + q2*A.y + t**2*A.z
   v

And the time derivative can be found with:

.. jupyter-execute::

   v.diff(t, A)

Lastly, vectors have a
:external:py:meth:`~sympy.physics.vector.vector.Vector.dt` method that
calculates time derivatives when viewed from a reference frame, saving a few
characters of typing:

.. jupyter-execute::

   v.dt(A)

We will use time derivatives in the next chapters to formulate velocity and
acceleration.
