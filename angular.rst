==================
Angular Kinematics
==================

Angular Velocity
================

If the relative orientation of two reference frames change with respect to
time, then we can calculate the angular velocity :math:`{}^A\bar{\omega}^B` of
reference frame :math:`B` in refrence frame :math:`A`. If
:math:`\hat{b}_x,\hat{b}_y,\hat{b}_z` are right handed mutally perpendicular
unit vectors fixed in :math:`B` then the angular velocity of :math:`B` in
:math:`A` is defined as ([Kane1985]_, pg. 16):

.. math::

   {}^A\bar{\omega}^B :=
   \left(\frac{{}^A d\hat{b}_y}{dt} \cdot \hat{b}_z\right) \hat{b}_x +
   \left(\frac{{}^A d\hat{b}_z}{dt} \cdot \hat{b}_x\right) \hat{b}_y +
   \left(\frac{{}^A d\hat{b}_x}{dt} \cdot \hat{b}_y\right) \hat{b}_z

If :math:`B` is oriented with respect to :math:`A` and mutually perpendicular
unit vectors :math:`\hat{a}_x,\hat{a}_y,\hat{a}_z` are fixed in :math:`A` then
there are these general relationships among the unit vectors of each frame:

.. math::

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
   B.orient_explicit(A, B_C_A.T)
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

Each of the measure numbers of :math:`{}^B\omega^A` are then:

.. jupyter-execute::

   mnx = me.dot(B.y.express(A).dt(A), B.z)
   mnx

.. jupyter-execute::

   mny = me.dot(B.z.express(A).dt(A), B.x)

.. jupyter-execute::

   mnz = me.dot(B.x.express(A).dt(A), B.y)
   mnz

The angular velocity vector is then:

.. jupyter-execute::

   A_w_B = mnx*B.x + mny*B.y + mnz*B.z
   A_w_B

For a simple rotation about the :math:`z` axis through :math:`\theta` the
direction cosine matrix is:

.. jupyter-execute::

   theta = me.dynamicsymbols('theta')

   B_C_A = sm.Matrix([[sm.cos(theta), sm.sin(theta), 0],
                      [-sm.sin(theta), sm.cos(theta), 0],
                      [0, 0, 1]])
   B_C_A

Following the same pattern as before the angular velocity of :math:`B` in
:math:`A` can be formed:

.. note::

   Don't confuse the left and rigth superscripts on direction cosine matrices
   and angular velocities. :math:`{}^B\mathbf{C}^A` describes the orientation
   of :math:`B` rotated with respect to :math:`A` and the mapping of vectors in
   :math:`A` to vectors expressed in :math:`B`. Whereas :math:`{}^A\omega^B`
   describes the angular velocity of :math:`B` when observed from :math:`A`.

.. jupyter-execute::

   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')
   A.orient_explicit(B, B_C_A)

   mnx = me.dot(B.y.express(A).dt(A), B.z)
   mny = me.dot(B.z.express(A).dt(A), B.x)
   mnz = me.dot(B.x.express(A).dt(A), B.y)

   A_w_B = mnx*B.x + mny*B.y + mnz*B.z
   A_w_B.simplify()

The angular velocity of a simple rotation is simply the time rate of change of
:math:`\theta` about :math:`\hat{b}_z`, the axis of rotation.

Similarly, if you establish a Euler :math:`z\textrm{-}x\textrm{-}z` orientation
with angles :math:`\psi,\theta,\varphi` respectively, then the angular velocity
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
