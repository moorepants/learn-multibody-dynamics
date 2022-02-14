======================
Vector Differentiation
======================

.. warning::

   This page is a draft until February 19, 2022. Report issues at:
   https://github.com/moorepants/learn-multibody-dynamics/issues

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`differentiation` or Jupyter Notebook:
   :jupyter-download:notebook:`differentation`.

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   sm.init_printing(use_latex='mathjax')

If a vector :math:`\bar{v}` is a function of scalar variables
:math:`q_1,q_2,\ldots,q_n` in reference frame :math:`A` then the first partial
derivatives of :math:`\bar{v}` in :math:`A` can be formed. We define the rth
partial derivative as:

.. math::

   \frac{{}^A\partial \bar{v}}{\partial q_r} :=
   \sum_{i=1}^3 \frac{\partial v_i}{\partial q_r} \hat{a}_i =
   \frac{\partial v_x}{\partial q_r} \hat{a}_x +
   \frac{\partial v_y}{\partial q_r} \hat{a}_y +
   \frac{\partial v_z}{\partial q_r} \hat{a}_z


where :math:`v_i` are the measure numbers of :math:`\bar{v}` expressed in
:math:`A` with mutually perpendicular unit vectors
:math:`\hat{a}_1,\hat{a}_2,\hat{a}_3`. This definition relies on the fact that
the unit vectors are fixed in :math:`A`.

Many of the vectors we will work with in multibody dynamics will be a function
of a single variable, time :math:`t`. If that is the case then the partial
derivatives reduce to a single variate derivative:

.. math::

   \frac{{}^A d \bar{v}}{dt} :=
   \sum_{i=1}^3 \frac{d v_i}{dt} \hat{a}_i =
   \sum_{i=1}^3 \dot{v}_i \hat{a}_i

.. warning::

   A dertiave written as :math:`\frac{\partial \bar{v}}{\partial q_r}` is
   meaningless because no reference frame is indicated. The derivative is
   dependent on which reference frame the change is viewed from.

This implies that a vector must be expressed in the reference frame one is
observing the change from before calculating the partial derivatives of the
measure numbers. For example:

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

.. jupyter-execute::

   dvdax = v.express(A).dot(A.x).diff(alpha)
   dvdax

.. jupyter-execute::

   dvday = v.express(A).dot(A.y).diff(alpha)
   dvday

.. jupyter-execute::

   dvdaz = v.express(A).dot(A.z).diff(alpha)
   dvdaz

So :math:`\frac{{}^A\partial \bar{v}}{\partial \alpha}` is:

.. jupyter-execute::

   dvda = dvdax*A.x + dvday*A.y + dvdaz*A.z
   dvda

:external:py:meth:`~sympy.physics.vector.vector.Vector.diff`

.. todo:: Open an issue on SymPy about Vector.diff() producing unecessarily
   complex results (seemingly). Here if v.diff() is called it is a mess. If
   v.diff(alpha, A).express(A) it's even more of a mess.

:external:py:meth:`~sympy.physics.vector.vector.Vector.simplify`

.. jupyter-execute::

   v.diff(alpha, A)

.. jupyter-execute::

   v.diff(alpha, A).simplify()

.. jupyter-execute::

   v.diff(alpha, A).express(A).simplify()

Product Rule
============

If you consider taking the derivative of the vector :math:`\bar{v}` in a
reference frame that it is not expressed in you must use the product rule. For
example, for vector :math:`\bar{v}` expressed in :math:`A` taking the
derivative in :math:`N` gives:

.. math::

   \frac{{}^N\partial \bar{v}}{\partial q_r} =
   \frac{{}^N\partial v_x}{\partial q_r}\hat{a}_x + v_x \frac{{}^N\partial \hat{a}_x}{\partial q_r} +
   \frac{{}^N\partial v_y}{\partial q_r}\hat{a}_y + v_y \frac{{}^N\partial \hat{a}_y}{\partial q_r} +
   \frac{{}^N\partial v_z}{\partial q_r}\hat{a}_z + v_z \frac{{}^N\partial \hat{a}_z}{\partial q_r}

The three similar terms with scalar derivatives have the same interpretation of
the ones in the prior section.

.. math::

   \frac{{}^N\partial v_x}{\partial q_r}\hat{a}_x +
   \frac{{}^N\partial v_y}{\partial q_r}\hat{a}_y +
   \frac{{}^N\partial v_z}{\partial q_r}\hat{a}_z

This part is more interesting. The partial derivative of a unit vector depends
on how it changes. But unit vectors don't change in length, only in
orientation.

.. math::

   v_x \frac{{}^N\partial \hat{a}_x}{\partial q_r} +
   v_y \frac{{}^N\partial \hat{a}_y}{\partial q_r} +
   v_z \frac{{}^N\partial \hat{a}_z}{\partial q_r}

We will see in the next chapter how to interpret and use these terms to
simplify the calculations of common derivatives.

The product rule also applies to the dot and cross products:

.. math::

   \frac{\partial}{\partial q_r}(\bar{v} \cdot \bar{w}) = &
   \frac{\partial \bar{v}}{\partial q_r} \cdot \bar{w} +
   \bar{v} \cdot \frac{\partial \bar{w}}{\partial q_r}

   \frac{\partial}{\partial q_r}(\bar{v} \times \bar{w}) = &
   \frac{\partial \bar{v}}{\partial q_r} \times \bar{w} +
   \bar{v} \times \frac{\partial \bar{w}}{\partial q_r}

This generalizes to any series of products. Let
:math:`P=F_1\cdot\ldots\cdot F_n` be a series of products. Then

.. math::

   \frac{\partial P}{\partial q_r} =
   \frac{\partial F_1}{\partial q_r}\cdot F_2 \cdot\ldots\cdot F_n +
   F_1 \cdot\frac{\partial F_2}{\partial q_r}\cdot F_3 \cdot\ldots\cdot F_n +
   \ldots +
   F_1 \cdot \ldots \cdot F_{n-1}\cdot  \frac{\partial F_n}{\partial q_r}

Second Derivatives
==================

If :math:`\frac{\partial \bar{v}}{\partial q_r}` is a vector function both in A
and any other reference frame it can change and be differentiated with respect
to any variable :math:`q_i` in any reference frame. We then arrive at the
second parital derivative

.. math::

   \frac{{}^B\partial}{\partial q_s} \left(\frac{{}^A\partial\bar{v}}{\partial
   q_r}\right)

Second partials in different reference frames do not commute:

.. math::

   \frac{{}^B\partial}{\partial q_s} \left(\frac{{}^A\partial\bar{v}}{\partial
   q_r}\right)
   \neq
   \frac{{}^A\partial}{\partial q_r} \left(\frac{{}^B\partial\bar{v}}{\partial
   q_s}\right)

If the reference frames of each diervative is the same, then mixed partials do
commute.

Vector Functions of Time
========================

In multibody dynamics we are primarily concern with how motion changes with
respect to time :math:`t` and our vectors and measure numbers will be implicit
functions of time, i.e. :math:`q_r(t)`. When that is the case the chain rule
can be used to take total derivatives.

.. math::

   \frac{{}^A d\bar{v}}{dt} = \sum_{i=1}^n \frac{{}^A\partial \bar{v}}{\partial
   q_r} \frac{d q_r}{dt} + \frac{{}^A \partial \bar{v}}{\partial t}

.. note::

   We will typically write :math:`\frac{dq}{dt}` as :math:`\dot{q}` and
   :math:`\frac{d^2q}{dt^2}` as :math:`\ddot{q}` and so on.

In SymPy Mechanics, scalar functions of time can be created like so:

.. jupyter-execute::

   t = sm.symbols('t')
   q_of = sm.Function('q')
   q = q_of(t)
   q

.. jupyter-execute::

   q.diff(t)

:external:py:func:`~sympy.physics.vector.dynamicsymbols`

.. jupyter-execute::

   q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')
   q1, q2, q3

.. jupyter-execute::

   t = me.dynamicsymbols._t

:external:py:func:`~sympy.physics.vector.printing.init_vprinting`

.. jupyter-execute::

   me.init_vprinting(use_latex='mathjax')
   q1.diff(t), q2.diff(t), q3.diff(t)
