==================
Angular Kinematics
==================

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`angular` or Jupyter Notebook:
   :jupyter-download:notebook:`angular`.

Angular Velocity
================

In Ch. :ref:`orientation` we learned that reference frames can be oriented
relative to each other. If the relative orientation of two reference frames
change with respect to time, then we can calculate the angular velocity
:math:`{}^A\bar{\omega}^B` of reference frame :math:`B` when observed from
reference frame :math:`A`. If :math:`\hat{b}_x,\hat{b}_y,\hat{b}_z` are right
handed mutually perpendicular unit vectors fixed in :math:`B` then the angular
velocity of :math:`B` in :math:`A` is defined as ([Kane1985]_, pg. 16):

.. math::
   :label: angular-velocity-definition

   {}^A\bar{\omega}^B :=
   \left(\frac{{}^A d\hat{b}_y}{dt} \cdot \hat{b}_z\right) \hat{b}_x +
   \left(\frac{{}^A d\hat{b}_z}{dt} \cdot \hat{b}_x\right) \hat{b}_y +
   \left(\frac{{}^A d\hat{b}_x}{dt} \cdot \hat{b}_y\right) \hat{b}_z

If :math:`B` is oriented with respect to :math:`A` and mutually perpendicular
unit vectors :math:`\hat{a}_x,\hat{a}_y,\hat{a}_z` are fixed in :math:`A` then
there are these general relationships among the unit vectors of each frame (see
:ref:`direction-cosine-matrix`):

.. math::
   :label: unit-vector-general-relation

   \hat{b}_x & = c_{xx} \hat{a}_x + c_{xy} \hat{a}_y + c_{xz} \hat{a}_z \\
   \hat{b}_y & = c_{yx} \hat{a}_x + c_{yy} \hat{a}_y + c_{yz} \hat{a}_z \\
   \hat{b}_z & = c_{zx} \hat{a}_x + c_{zy} \hat{a}_y + c_{zz} \hat{a}_z

We can create these equations in SymPy to demonstrate how to work with the
definition of angular velocity.

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   me.init_vprinting(use_latex='mathjax')

We first create the direction cosine matrix with time varying elements:

.. jupyter-execute::

   cxx, cyy, czz = me.dynamicsymbols('c_{xx}, c_{yy}, c_{zz}')
   cxy, cxz, cyx = me.dynamicsymbols('c_{xy}, c_{xz}, c_{yx}')
   cyz, czx, czy = me.dynamicsymbols('c_{yz}, c_{zx}, c_{zy}')

   B_C_A = sm.Matrix([[cxx, cxy, cxz],
                      [cyx, cyy, cyz],
                      [czx, czy, czz]])

and establish the orientation using
:external:py:meth:`~sympy.physics.vector.frame.ReferenceFrame.orient_explicit`:

.. warning::

   Remember this method takes the transpose of the direction cosine matrix.

.. jupyter-execute::

   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')
   B.orient_explicit(A, B_C_A.transpose())
   B.dcm(A)

This now let's us write the :math:`B` unit vectors in terms of the :math:`A`
unit vectors:

.. jupyter-execute::

   B.x.express(A)

.. jupyter-execute::

   B.y.express(A)

.. jupyter-execute::

   B.z.express(A)

Recalling the definition of angular velocity above, each of the measure numbers
of the angular velocity is calculated by dotting the derivative of a :math:`B`
unit vector in :math:`A` with a unit vector in :math:`B`. :math:`\frac{{}^A
\hat{b}_y}{dt}` is for example:

.. jupyter-execute::

   B.y.express(A).dt(A)

Each of the measure numbers of :math:`{}^A\bar{\omega}^B` are then:

.. jupyter-execute::

   mnx = me.dot(B.y.express(A).dt(A), B.z)
   mnx

.. jupyter-execute::

   mny = me.dot(B.z.express(A).dt(A), B.x)
   mny

.. jupyter-execute::

   mnz = me.dot(B.x.express(A).dt(A), B.y)
   mnz

The angular velocity vector is then:

.. jupyter-execute::

   A_w_B = mnx*B.x + mny*B.y + mnz*B.z
   A_w_B

Simple Orientations
===================

For a simple orientation about the :math:`z` axis through :math:`\theta` the
direction cosine matrix is:

.. jupyter-execute::

   theta = me.dynamicsymbols('theta')

   B_C_A = sm.Matrix([[sm.cos(theta), sm.sin(theta), 0],
                      [-sm.sin(theta), sm.cos(theta), 0],
                      [0, 0, 1]])
   B_C_A

Following the same pattern as before the angular velocity of :math:`B` in
:math:`A` can be formed:

.. jupyter-execute::

   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')
   A.orient_explicit(B, B_C_A)

   mnx = me.dot(B.y.express(A).dt(A), B.z)
   mny = me.dot(B.z.express(A).dt(A), B.x)
   mnz = me.dot(B.x.express(A).dt(A), B.y)

   A_w_B = mnx*B.x + mny*B.y + mnz*B.z
   A_w_B.simplify()

.. note::

   Don't confuse the left and right superscripts on direction cosine matrices
   and angular velocities. :math:`{}^B\mathbf{C}^A` describes the orientation
   of :math:`B` rotated with respect to :math:`A` and the mapping of vectors in
   :math:`A` to vectors expressed in :math:`B`. Whereas
   :math:`{}^A\bar{\omega}^B` describes the angular velocity of :math:`B` when
   observed from :math:`A`.

The angular velocity of a simple orientation is simply the time rate of change of
:math:`\theta` about :math:`\hat{b}_z`, the axis of the simple orientation. SymPy
Mechanics offers the
:external:py:meth:`~sympy.physics.vector.frame.ReferenceFrame.ang_vel_in`
method for automatically calculating the angular velocity if a direction cosine
matrix exists between the two reference frames:

.. jupyter-execute::

   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')
   B.orient_axis(A, theta, A.z)
   B.ang_vel_in(A)

.. todo:: Should this return the angular velocity expressed in the body fixed
   frame?

A simple orientation and associated simple angular velocity can be formulated for
any arbitrary orientation axis vector, not just one of the three mutually
perpendicular unit vectors as shown above. There is a simple angular velocity
between two reference frames :math:`A` and :math:`B` if there exists a single
unit vector :math:`\hat{k}` which is fixed in both :math:`A` and :math:`B` for
some finite time. If this is the case, then :math:`{}^A\bar{\omega}^B = \omega
\hat{k}` where :math:`\omega` is the time rate of change of the angle
:math:`\theta` between a line fixed in :math:`A` and another line fixed in
:math:`B` both of which are perpendicular to the orientation axis :math:`\hat{k}`.
We call :math:`\omega=\dot{\theta}` the angular speed of :math:`B` in
:math:`A`.

:external:py:meth:`~sympy.physics.vector.frame.ReferenceFrame.orient_axis` can
take any arbitrary vector fixed in :math:`A` and :math:`B` to establish the
orientation:

.. jupyter-execute::

   theta = me.dynamicsymbols('theta')

   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')
   B.orient_axis(A, theta, A.x + A.y)
   B.ang_vel_in(A)

The angular speed is then:

.. jupyter-execute::

   B.ang_vel_in(A).magnitude()

.. note:: This result should be :math:`|\dot{\theta}|`. This is a bug in SymPy,
   see https://github.com/sympy/sympy/issues/23173 for more info.

.. todo:: Why doesn't this simplify to theta dot? I tried ``real=True`` on
   theta.

Body Fixed Orientations
=======================

If you establish a Euler :math:`z\textrm{-}x\textrm{-}z` orientation with
angles :math:`\psi,\theta,\varphi` respectively, then the angular velocity
vector is:

.. jupyter-execute::

   psi, theta, phi = me.dynamicsymbols('psi, theta, varphi')

   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')
   B.orient_body_fixed(A, (psi, theta, phi), 'ZXZ')

   mnx = me.dot(B.y.express(A).dt(A), B.z)
   mny = me.dot(B.z.express(A).dt(A), B.x)
   mnz = me.dot(B.x.express(A).dt(A), B.y)

   A_w_B = mnx*B.x + mny*B.y + mnz*B.z
   A_w_B.simplify()

.. todo::

   ``simplify()`` shouldn't be needed here: https://github.com/sympy/sympy/issues/23130

:external:py:meth:`~sympy.physics.vector.frame.ReferenceFrame.ang_vel_in` gives
the same result:

.. jupyter-execute::

   B.ang_vel_in(A).simplify()

Time Derivatives of Vectors
===========================

Using the definition of angular velocity one can show ([Kane1985]_, pg. 17)
that the time derivative of a unit vector **fixed in** :math:`B` is related to
:math:`B`'s angular velocity as so:

.. math::
   :label: time-derivative-fixed-unit-vector

   \frac{{}^Ad\hat{b}_x}{dt} = {}^A\bar{\omega}^B \times \hat{b}_x

This shows that the derivative is always normal to the rotating unit vector,
because the magnitude of the unit vector is constant, and the derivative scales
with the magnitude of the angular velocity:

.. math::
   :label: time-derivative-unit-vector-scalar-mag

   \frac{{}^Ad\hat{b}_x}{dt} = \left| {}^A\bar{\omega}^B \right| \left( {}^A\hat{\omega}^B \times \hat{b}_x \right)

Now if vector :math:`\bar{v} = v\hat{b}_x` and :math:`v` is constant with
respect to time then:

.. math::
   :label: time-derivative-fixed-vector

   \frac{{}^A d\bar{v}}{dt} =
   v({}^A\bar{\omega}^B \times \hat{b}_x) =
   {}^A\bar{\omega}^B \times v\hat{b}_x =
   {}^A\bar{\omega}^B \times \bar{v}

This extends to any vector **fixed in** :math:`B` and observed from :math:`A`,
making the time derivative equal to the cross product of the angular velocity
of :math:`B` in :math:`A` with the vector.

Now, if :math:`\bar{u}` is a vector that is **not fixed in** :math:`B` we
return to the product rule in Section :ref:`product-rule`. First expressed
:math:`\bar{u}` in :math:`B`:

.. math::
   :label: time-varying-vector

   \bar{u} = u_1\hat{b}_x + u_2\hat{b}_y + u_3\hat{b}_z

The derivative in another reference frame :math:`A` is then:

.. math::
   :label: deriv-arb-vector

   \frac{{}^Ad\bar{u}}{dt} &=
   \dot{u}_1\hat{b}_x + \dot{u}_2\hat{b}_y + \dot{u}_3\hat{b}_z +
   u_1\frac{{}^Ad\hat{b}_x}{dt} + u_2\frac{{}^Ad\hat{b}_y}{dt} + u_3\frac{{}^Ad\hat{b}_z}{dt} \\
   &=
   \frac{{}^Bd\bar{u}}{dt} +
   u_1{}^A\bar{\omega}^B\times\hat{b}_x + u_2{}^A\bar{\omega}^B\times\hat{b}_y + u_3{}^A\bar{\omega}^B\times\hat{b}_z \\
   &=
   \frac{{}^Bd\bar{u}}{dt} +
   {}^A\bar{\omega}^B\times\bar{u}

We can show that Eq. :math:numref:`deriv-arb-vector` holds with an example.
Take a :math:`z\textrm{-}x` orientation and an arbitrary vector that is not fixed
in :math:`B`:

.. jupyter-execute::

   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')
   B.orient_body_fixed(A, (psi, theta, 0), 'ZXZ')

   u1, u2, u3 = me.dynamicsymbols('u1, u2, u3')

   u = u1*B.x + u2*B.y + u3*B.z
   u

As we learned in the last chapter we can express the vector in :math:`A` and
then take the time derivative of the measure numbers to arrive at
:math:`\frac{{}^Ad\bar{u}}{dt}`:

.. jupyter-execute::

   u.express(A)

.. jupyter-execute::

   u.express(A).dt(A)

But applying the theorem above we can find the derivative with a cross product.
The nice aspect of this formulation is there is no need to express the vector
in :math:`A`. First :math:`\frac{{}^Bd\bar{u}}{dt}`:

.. jupyter-execute::

   u.dt(B)

and then :math:`{}^A\bar{\omega}^B\times\bar{u}`:

.. jupyter-execute::

   A_w_B = B.ang_vel_in(A)
   A_w_B

:math:`\frac{{}^Ad\bar{u}}{dt}` is then:

.. jupyter-execute::

   u.dt(B) + me.cross(A_w_B, u)

We can show that the first result is equivalent by expressing in :math:`B` and
simplifying:

.. jupyter-execute::

   u.express(A).dt(A).express(B).simplify()

.. _addition-of-angular-velocity:

Addition of Angular Velocity
============================

Similar to the relationship in direction cosine matrices of successive
orientations (Sec. :ref:`successive-orientations`), there is a relationship
among the angular velocities of successively oriented reference frames
([Kane1985]_, pg. 24) but it is tied to the addition of vectors instead of
multiplication of matrices.

.. math::
   :label: addition-angular-velocity

   {}^A\bar{\omega}^Z =
   {}^A\bar{\omega}^B +
   {}^B\bar{\omega}^C +
   \ldots +
   {}^Y\bar{\omega}^Z

We can demonstrate this by creating three simple orientations for a Euler
:math:`y\textrm{-}x\textrm{-}y` orientation:

.. jupyter-execute::

   psi, theta, phi = me.dynamicsymbols('psi, theta, varphi')

   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')
   C = me.ReferenceFrame('C')
   D = me.ReferenceFrame('D')

   B.orient_axis(A, psi, A.y)
   C.orient_axis(B, theta, B.x)
   D.orient_axis(C, phi, C.y)

The simple angular velocity of each successive orientation is shown:

.. jupyter-execute::

   A_w_B = B.ang_vel_in(A)
   A_w_B

.. jupyter-execute::

   B_w_C = C.ang_vel_in(B)
   B_w_C

.. jupyter-execute::

   C_w_D = D.ang_vel_in(C)
   C_w_D

Summing the successive angular velocities gives the compact result:

.. jupyter-execute::

   A_w_D = A_w_B + B_w_C + C_w_D
   A_w_D

Similarly, we can skip the auxiliary frames and form the relationship between
:math:`A` and :math:`D` directly and calculate :math:`{}^A\bar{\omega}^D`:

.. jupyter-execute::

   A2 = me.ReferenceFrame('A')
   D2 = me.ReferenceFrame('D')
   D2.orient_body_fixed(A2, (psi, theta, phi), 'YXY')
   D2.ang_vel_in(A2).simplify()

If we express our prior result in :math:`D` we see the results are the same:

.. jupyter-execute::

   A_w_D.express(D)

.. todo:: I could show with three generic direction cosine matrices that the
   angular velocities add up, but that would be a bit messy presentation.

Angular Acceleration
====================

The angular acceleration of :math:`B` when observed from :math:`A` is defined
as:

.. math::
   :label: angular-acceleration-definition

   {}^A\bar{\alpha}^B := \frac{{}^Ad}{dt} {}^A\bar{\omega}^B

:math:`{}^A\bar{\omega}^B` is simply a vector so we can time differentiate it
with respect to frame :math:`A`. Using Eq. :math:numref:`deriv-arb-vector` we
can write:

.. math::
   :label: angular-acceleration-cross

   \frac{{}^Ad}{dt} {}^A\bar{\omega}^B & =
   \frac{{}^Bd}{dt} {}^A\bar{\omega}^B + {}^A\bar{\omega}^B \times {}^A\bar{\omega}^B \\

and since :math:`{}^A\bar{\omega}^B \times {}^A\bar{\omega}^B=0`:

.. math::
   :label: ang-acc-frame

   \frac{{}^Ad}{dt} {}^A\bar{\omega}^B = \frac{{}^Bd}{dt} {}^A\bar{\omega}^B

which is rather convenient.

With SymPy Mechanics :math:`{}^A\bar{\alpha}^B` is found automatically with
:external:py:meth:`~sympy.physics.vector.frame.ReferenceFrame.ang_acc_in` if
the orientations are established. For a simple orientation:

.. jupyter-execute::

   theta = me.dynamicsymbols('theta')

   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')
   B.orient_axis(A, theta, A.z)
   B.ang_acc_in(A)

Similarly we can calcualte the derivative manually:

.. jupyter-execute::

   B.ang_vel_in(A).dt(A)

and see that that Eq. :math:numref:`ang-acc-frame` holds:

.. jupyter-execute::

   B.ang_vel_in(A).dt(B)

For a body fixed orientation we get:

.. jupyter-execute::

   psi, theta, phi = me.dynamicsymbols('psi, theta, varphi')

   A = me.ReferenceFrame('A')
   D = me.ReferenceFrame('D')
   D.orient_body_fixed(A, (psi, theta, phi), 'YXY')

   D.ang_acc_in(A).simplify()

and with manual derivatives:

.. jupyter-execute::

   D.ang_vel_in(A).dt(A).simplify()

.. jupyter-execute::

   D.ang_vel_in(A).dt(D).simplify()

Addition of Angular Acceleration
================================

The calculation of angular acceleration is relatively simple, but the addition
of angular velocities explained in Sec. :ref:`addition-of-angular-velocity`
does not extend to angular accelerations.

.. math::
   :label: addition-angular-acceleration

   {}^A\bar{\alpha}^Z \neq
   {}^A\bar{\alpha}^B +
   {}^B\bar{\alpha}^C +
   \ldots +
   {}^Y\bar{\alpha}^Z

Coming back to the successive orientations that form a
:math:`y\textrm{-}x\textrm{-}y` Euler rotation, we can see that the result is
not the same as above:

.. jupyter-execute::

   psi, theta, phi = me.dynamicsymbols('psi, theta, varphi')

   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')
   C = me.ReferenceFrame('C')
   D = me.ReferenceFrame('D')

   B.orient_axis(A, psi, A.y)
   C.orient_axis(B, theta, B.x)
   D.orient_axis(C, phi, C.y)

The simple angular acceleration of each successive orientations is shown:

.. jupyter-execute::

   A_alp_B = B.ang_acc_in(A)
   A_alp_B

.. jupyter-execute::

   B_alp_C = C.ang_acc_in(B)
   B_alp_C

.. jupyter-execute::

   C_alp_D = D.ang_acc_in(C)
   C_alp_D

Summing the successive angular accelerations gives this result:

.. jupyter-execute::

   A_alp_D = A_alp_B + B_alp_C + C_alp_D
   A_alp_D.express(D).simplify()

which is not equal to the correct, more complex, result:

.. jupyter-execute::

   D.ang_acc_in(A).express(D).simplify()

.. warning:: Do not sum successive angular accelerations!
