======================
Vector Differentiation
======================

.. warning::

   This page is a draft until February 19, 2022. Report issues at:
   https://github.com/moorepants/learn-multibody-dynamics/issues

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
derivative of :math:`\bar{v}` in :math:`A` with respect to :math:`q_r` where
:math:`r=1\ldots n` can be formed with the following definition:

.. math::

   \frac{{}^A\partial \bar{v}}{\partial q_r} :=
   \sum_{i=1}^3 \frac{\partial v_i}{\partial q_r} \hat{a}_i

where :math:`v_i` are the measure numbers of :math:`\bar{v}` expressed in
:math:`A` with mutually perpendicular unit vectors
:math:`\hat{a}_1,\hat{a}_2,\hat{a}_3`. This definition relies on the fact that
the unit vectors are fixed in and do not change when observed from :math:`A`.

Given :math:`\bar{v}=v_x\hat{a}_x+v_y\hat{a}_y+v_z\hat{a}_z` the above
definition expands to:

.. math::

   \frac{{}^A\partial \bar{v}}{\partial q_r} =
   \frac{\partial v_x}{\partial q_r} \hat{a}_x +
   \frac{\partial v_y}{\partial q_r} \hat{a}_y +
   \frac{\partial v_z}{\partial q_r} \hat{a}_z

Many of the vectors we will work with in multibody dynamics will be a function
of a single variable, most often time :math:`t`. If that is the case, the
partial derivative reduces to a single variate derivative:

.. math::

   \frac{{}^A d \bar{v}}{dt} := \sum_{i=1}^3 \frac{d v_i}{dt} \hat{a}_i

.. warning::

   A derivative written as :math:`\frac{\partial \bar{v}}{\partial q_r}` is
   meaningless because no reference frame is indicated. The derivative is
   dependent on which reference frame the change is observed from, so without a
   reference frame, the derivative cannot be calculated.

The above definition implies that a vector must be expressed in the reference
frame one is observing the change from before calculating the partial
derivatives of the scalar measure numbers. For example here is a vector that is
expressed in three different reference frames:

.. jupyter-execute::

   a, b, c, d, e, f = sm.symbols('a, b, c, d, e, f')
   alpha, beta = sm.symbols('alpha, beta')

   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')
   C = me.ReferenceFrame('C')

   B.orient_axis(A, alpha, A.x)
   C.orient_axis(B, beta, B.y)

   v = a*A.x + b*A.y + c*B.x + d*B.y + e*C.x + f*C.y
   v

To calculate

.. math::

   \frac{{}^A\partial\bar{v}}{\partial \alpha}

we first need to project the vector onto the unit vectors of :math:`A` and take
the partial derivative of those measure numbers with respect to :math:`\alpha`.


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
\alpha}` from the new measure numbers:

.. jupyter-execute::

   dvda = dvdax*A.x + dvday*A.y + dvdaz*A.z
   dvda

..
   .. todo:: Open an issue on SymPy about Vector.diff() producing unecessarily
      complex results (seemingly). Here if v.diff() is called it is a mess. If
      v.diff(alpha, A).express(A) it's even more of a mess.

SymPy Mechanics vectors have a special
:external:py:meth:`~sympy.physics.vector.vector.Vector.diff` method that
manages taking partial derivatives from different reference frames. For the
vector ``.diff()`` you provide first the variable :math:`\alpha` followed by
the reference frame:

.. jupyter-execute::

   dvdalpha = v.diff(alpha, A)
   dvdalpha

The result is not so simplified so you can use the vector
:external:py:meth:`~sympy.physics.vector.vector.Vector.simplify` method, which
applies ``trigsimp()`` to each measure number:

.. jupyter-execute::

   v.diff(alpha, A).simplify()

This multi reference frame form can be shown to be the same as we calculated
above by expressing it fully in :math:`A` and simplifying:

.. jupyter-execute::

   v.diff(alpha, A).express(A).simplify()

Product Rule
============

Consider again vector :math:`\bar{v}=v_x\hat{a}_x+v_y\hat{a}_y+v_z\hat{a}_z`.
If reference frames :math:`N` and :math:`A` are oriented relative to each other
and the orientation is also a function of :math:`q_r` then we must use the
product rule when taking the partial derivative. For example:

.. math::

   \frac{{}^N\partial \bar{v}}{\partial q_r} =
   \frac{{}^N\partial v_x}{\partial q_r}\hat{a}_x + v_x \frac{{}^N\partial \hat{a}_x}{\partial q_r} +
   \frac{{}^N\partial v_y}{\partial q_r}\hat{a}_y + v_y \frac{{}^N\partial \hat{a}_y}{\partial q_r} +
   \frac{{}^N\partial v_z}{\partial q_r}\hat{a}_z + v_z \frac{{}^N\partial \hat{a}_z}{\partial q_r}

The three similar terms with scalar derivatives have the same interpretation of
the ones in the prior section.

.. math::

   \frac{{}^N\partial v_x}{\partial q_r}\hat{a}_x,
   \frac{{}^N\partial v_y}{\partial q_r}\hat{a}_y,
   \frac{{}^N\partial v_z}{\partial q_r}\hat{a}_z

But the part with unit vector derivatives is more interesting. The partial
derivative of a unit vector depends on how it changes. But unit vectors do not
change in length, only in orientation.

.. math::

   v_x \frac{{}^N\partial \hat{a}_x}{\partial q_r},
   v_y \frac{{}^N\partial \hat{a}_y}{\partial q_r},
   v_z \frac{{}^N\partial \hat{a}_z}{\partial q_r}

You will learn in the next chapter how to interpret and use these terms to
simplify the calculations of common derivatives.

The product rule also applies to the dot and cross products:

.. math::

   \frac{\partial}{\partial q_r}(\bar{v} \cdot \bar{w}) = &
   \frac{\partial \bar{v}}{\partial q_r} \cdot \bar{w} +
   \bar{v} \cdot \frac{\partial \bar{w}}{\partial q_r}

   \frac{\partial}{\partial q_r}(\bar{v} \times \bar{w}) = &
   \frac{\partial \bar{v}}{\partial q_r} \times \bar{w} +
   \bar{v} \times \frac{\partial \bar{w}}{\partial q_r}

This generalizes to any series of products. Let :math:`G=f_1\cdot\ldots\cdot
f_n` be a series of products. Then:

.. math::

   \frac{\partial G}{\partial q_r} =
   \frac{\partial f_1}{\partial q_r}\cdot f_2 \cdot\ldots\cdot f_n +
   f_1 \cdot\frac{\partial f_2}{\partial q_r}\cdot f_3 \cdot\ldots\cdot f_n +
   \ldots +
   f_1 \cdot \ldots \cdot f_{n-1}\cdot  \frac{\partial f_n}{\partial q_r}

Second Derivatives
==================

If :math:`\frac{{}^A\partial \bar{v}}{\partial q_r}` is a vector function both
in :math:`A` and also any other reference frame it can then change and thus be
differentiated with respect to any variable :math:`q_r` in any reference frame,
e.g. :math:`B`. We then arrive at the second partial derivative:

.. math::

   \frac{{}^B\partial}{\partial q_s} \left(\frac{{}^A\partial\bar{v}}{\partial
   q_r}\right)

Second partials in different reference frames do not necessarily commute:

.. math::

   \frac{{}^B\partial}{\partial q_s} \left(\frac{{}^A\partial\bar{v}}{\partial
   q_r}\right)
   \neq
   \frac{{}^A\partial}{\partial q_r} \left(\frac{{}^B\partial\bar{v}}{\partial
   q_s}\right)

If the reference frames of each partial derivative is the same, then mixed
partials do commute.

Vector Functions of Time
========================

In multibody dynamics we are primarily concerned with how motion changes with
respect to time :math:`t` and our vectors and measure numbers will often be
implicit functions of time, i.e. :math:`q_r(t)`. When that is the case the
chain rule can be used to take total derivatives:

.. math::

   \frac{{}^A d\bar{v}}{dt} = \sum_{i=1}^n \frac{{}^A\partial v_i}{\partial
   q_r(t)} \frac{d q_r(t)}{dt} + \frac{{}^A \partial \bar{v}}{\partial t}

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
