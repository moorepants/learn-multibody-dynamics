======================
Vector Differentiation
======================

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`differentiation` or Jupyter Notebook:
   :jupyter-download:notebook:`differentiation`.

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   sm.init_printing(use_latex='mathjax')

Partial Derivatives
===================

If a vector :math:`\bar{v}` is a function of :math:`n` scalar variables
:math:`q_1,q_2,\ldots,q_n` in reference frame :math:`A` then the first partial
derivatives of :math:`\bar{v}` in :math:`A` with respect to :math:`q_r` where
:math:`r=1\ldots n` can be formed with the following definition:

.. math::
   :label: partial-deriv-def

   \frac{{}^A\partial \bar{v}}{\partial q_r} :=
   \sum_{i=1}^3 \frac{\partial v_i}{\partial q_r} \hat{a}_i

where :math:`v_i` are the measure numbers of :math:`\bar{v}` expressed in
:math:`A` with mutually perpendicular unit vectors
:math:`\hat{a}_1,\hat{a}_2,\hat{a}_3`. This definition relies on the fact that
the unit vectors are fixed in :math:`A` and thus do not change when observed
from :math:`A`.

Given :math:`\bar{v}=v_x\hat{a}_x+v_y\hat{a}_y+v_z\hat{a}_z` the above
definition expands to:

.. math::
   :label: partial-measure

   \frac{{}^A\partial \bar{v}}{\partial q_r} =
   \frac{\partial v_x}{\partial q_r} \hat{a}_x +
   \frac{\partial v_y}{\partial q_r} \hat{a}_y +
   \frac{\partial v_z}{\partial q_r} \hat{a}_z

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
   reference frame, the derivative cannot be calculated.

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

   dvdax = v.dot(A.x).diff(alpha)
   dvdax

.. jupyter-execute::

   dvday = v.dot(A.y).diff(alpha)
   dvday

.. jupyter-execute::

   dvdaz = v.dot(A.z).diff(alpha)
   dvdaz

We can then construct the vector :math:`\frac{{}^A\partial \bar{v}}{\partial
\alpha}` from the new measure numbers know that the :math:`A` unit vectors are
fixed:

.. jupyter-execute::

   dvda = dvdax*A.x + dvday*A.y + dvdaz*A.z
   dvda

SymPy Mechanics vectors have a special
:external:py:meth:`~sympy.physics.vector.vector.Vector.diff` method that
manages taking partial derivatives from different reference frames. For the
vector ``.diff()`` method you provide first the variable :math:`\alpha`
followed by the reference frame you are observing from:

.. jupyter-execute::

   dvdalpha = v.diff(alpha, A)
   dvdalpha

The result is not so simplified because SymPy attempts to express the
derivative in the same components as the vector was, so you can use the vector
:external:py:meth:`~sympy.physics.vector.vector.Vector.simplify` method, which
applies :external:py:func:`~sympy.simplify.trigsimp.trigsimp` to each measure
number:

.. jupyter-execute::

   v.diff(alpha, A).simplify()

This multi reference frame form can be shown to be the same as we calculated
above by expressing it fully in :math:`A` and simplifying:

.. jupyter-execute::

   v.diff(alpha, A).express(A).simplify()

.. _product-rule:

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
respect to :math:`A` and then with respect to :math:`B`, the second partial
derivative is:

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

If the reference frames of each partial derivative are the same, then mixed
partials do commute.

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

SymPy Mechanics provides the convince function
:external:py:func:`~sympy.physics.vector.functions.dynamicsymbols` to create
scalar functions of time just like ``symbols()``:

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
calculates time derivatives, saving a few characters of typing:

.. jupyter-execute::

   v.dt(A)
