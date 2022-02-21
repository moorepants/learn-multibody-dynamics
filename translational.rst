.. |_| unicode:: 0xA0
   :trim:

=======================
Translational Kinematics
========================

.. warning::

   This page is a draft until February 25, 2022. Report issues at:
   https://github.com/moorepants/learn-multibody-dynamics/issues

In multibody dynamics, we are going to need to calculate the velocities and
accelerations of points. We will learn that the acceleration of the mass
centers of the bodies in a multibody system will be a primary ingredient in
forming Newton's Second Law of motion :math:`\bar{F} = m\bar{a}`. This chapter
will equip you to calculate the relative translational velocities and
acclerations of points in a system.

Translational Velocity
======================

If a point :math:`P` is moving with respect to a point :math:`O` that is fixed
in reference frame :math:`A` the translational velocity vector of point
:math:`P` is defined as:

.. math::

   {}^A\bar{v}^P := \frac{{}^Ad\bar{r}^{P/O}}{dt}

We also know from Eq. :math:numref:`deriv-arb-vector` that the time derivative
of any vector can be written in terms of the angular velocity of the associated
reference frames. So:

.. math::
   :label: point-velocity-two-frames

   {}^A\bar{v}^P =
   \frac{{}^Ad\bar{r}^{P/O}}{dt} =
   \frac{{}^Bd\bar{r}^{P/O}}{dt} +
   {}^A\bar{\omega}^B\times\bar{r}^{P/O} =
   {}^B\bar{v}^P + {}^A\bar{\omega}^B\times\bar{r}^{P/O}

This formulation allows us to utlize reference frames to simplify velocity
calculations. Take for example this piece of kinetic art that now stands in
Rotterdam:

.. figure:: https://upload.wikimedia.org/wikipedia/commons/thumb/0/03/Rickey_Rotterdam_04.JPG/360px-Rickey_Rotterdam_04.JPG
   :align: center

   Kinetic sculpture "Two Turning Vertical Rectangles"(1971) in Rotterdam/The
   Netherlands (FOP) by George Rickey.
   https://nl.wikipedia.org/wiki/Two_Turning_Vertical_Rectangles

   K.Siereveld, Public domain, via Wikimedia Commons

.. todo:: Need a non-breaking space between K. and Si...

User https://www.reddit.com/user/stravalnak posted this video of the sculpture
to Reddit during the 2022 storm Eunice:

.. raw:: html

   <center>
   <video width="320" controls>
     <source src="https://v.redd.it/egrhlwlbrmi81/DASH_1080.mp4"
     type="video/mp4">
   </video>
   <p> From: https://www.reddit.com/r/Rotterdam/comments/svo3cs/hij_maakt_overuren/</p>
   </center>

and it looks very dangerous. It would be interesting to know the velocity and
acceleration of various points on this sculpture. So we will use this as an
example.

.. figure:: figures/translational-kinetic-sculpture.svg

   Sketch of one of the two plates. Note the pigeon trying to walk across one
   edge of the plate at point :math:`R`.

   Pigeon SVG from https://freesvg.org/vector-clip-art-of-homing-pigeon Public DOmain

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   me.init_vprinting(use_latex='mathjax')

.. jupyter-execute::

   alpha, beta = me.dynamicsymbols('alpha, beta')

   N = me.ReferenceFrame('N')
   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')

   A.orient_axis(N, alpha, N.z)
   B.orient_axis(A, beta, A.x)

.. jupyter-execute::

   h, d, w, c, l = sm.symbols('h, d, w, c, l')

   r_O_P = h*N.z
   r_P_S = -d*A.x
   r_S_Q = -w*B.x - (c + l/2)*B.z

   r_O_P, r_P_S

.. math::

   {}^N\bar{v}^S = {}^A\bar{v}^S + {}^N\bar{\omega}^A\times\bar{r}^{S/O}

.. jupyter-execute::

   (r_O_P + r_P_S).dt(A)

.. jupyter-execute::

   A.ang_vel_in(N)

.. jupyter-execute::

   me.cross(A.ang_vel_in(N), r_O_P + r_P_S)

.. jupyter-execute::

   (r_O_P + r_P_S + r_S_Q).dt(B)

.. jupyter-execute::

   me.cross(B.ang_vel_in(N), r_O_P + r_P_S + r_S_Q)

SymPy Mechanics provites the
:external:py:class:`~sympy.physics.vector.point.Point` object that simplifies
working with position vectors.

.. jupyter-execute::

   O = me.Point('O')
   P = me.Point('P')
   S = me.Point('S')
   Q = me.Point('Q')

   P.set_pos(O, h*N.z)
   S.set_pos(P, -d*A.x)
   Q.set_pos(S, -w*B.x - (c + l/2)*B.z)

   Q.pos_from(O)

.. jupyter-execute::

   O.set_vel(N, 0)
   Q.vel(N)

Velocity Two Point Thereom
==========================

If there are two points :math:`P` and :math:`S` fixed in a reference frame
:math:`A` and you know the angular velocity :math:`{}^N\bar{\omega}^A` and the
velocity :math:`{}^N\bar{v}^P` then :math:`{}^N\bar{v}^S` can be caluculated if
the vector :math:`\bar{r}^{S/P}`, which is fixed in :math:`A`, is known. The
following thereom provides a conveniet formulation:

.. math::

   {}^N\bar{v}^S &=  \frac{{}^N d\bar{r}^{S/O} }{dt} \\
   &= \frac{{}^N d\left(\bar{r}^{P/O} + \bar{r}^{S/P}\right)}{dt} \\
   &= {}^A\bar{v}^P + \frac{{}^A d\bar{r}^{S/P} }{dt} \\
   &= {}^A\bar{v}^P + {}^N\bar{\omega}^A \times \bar{r}^{S/P}

Both :math:`O` and :math:`P` are fixed in :math:`N`, so
:math:`{}^N\bar{v}^P=0`.

.. jupyter-execute::

   N_v_P = 0*N.z

.. jupyter-execute::

   N_v_S = N_v_P +  me.cross(A.ang_vel_in(N), S.pos_from(P))
   N_v_S

Point objects have the
:external:py:meth:`~sympy.physics.vector.point.Point.v2pt_theory` method for
calculating the above equation given the other point fixed in the same frame,
the frame you want the velocity in, and the frame both points are fixed in. The
velocity of :math:`P` is set to zero using
:external:py:meth:`~sympy.physics.vector.point.Point.set_vel` first to ensure
we start with a known velocity.

.. jupyter-execute::

   P.set_vel(N, 0)
   S.v2pt_theory(P, N, A)

Note that when you
call:external:py:meth:`~sympy.physics.vector.point.Point.v2pt_theory` it also
sets the velocity of point :math:`S` to this value:

.. jupyter-execute::

   S.vel(N)

Both points :math:`S` and :math:`Q` are fixed in reference frame :math:`B` and
we just caclulated :math:`{}^N\bar{v}^S`, so we can use the two point theorem
to find the velocity of :math:`Q` in a similar fashion. First, manually:

.. math::

   {}^N\bar{v}^Q = {}^N\bar{v}^S + {}^N\bar{\omega}^B \times \bar{r}^{Q/S}

.. jupyter-execute::

   N_v_Q = N_v_S +  me.cross(B.ang_vel_in(N), Q.pos_from(S))
   N_v_Q

And with the :external:py:meth:`~sympy.physics.vector.point.Point.v2pt_theory`:

.. jupyter-execute::

   Q.v2pt_theory(S, N, B)

.. admonition:: Exercise

   Calculate the velocity of the center of mass of the plate using the two
   point theorem.

Velocity One Point Thereom
==========================

If you are interested in the velocity of a point :math:`R` that is moving in a reference
frame :math:`B` and you know the velocity of a point in that reference frame
then

.. math::

   {}^N\bar{v}^R = {}^B\bar{v}^R + {}^N\bar{v}^T

where point :math:`T` is a point that conicides with :math:`R` at that instant.

Combined with the the two point theorem for :math:`T` you get:

.. math::

   {}^N\bar{v}^R = {}^B\bar{v}^R + {}^N\bar{v}^S + {}^N\bar{\omega}^B \times \bar{r}^{T/S}

If the pigeon :math:`R` is walking at a constant rate :math:`s` in the
:math:`\hat{b}_x` direction, then we can calculate the velocity of the pigeon
when observed from the :math:`N` reference frame.

.. jupyter-execute::

   s = me.dynamicsymbols('s')
   t = me.dynamicsymbols._t

   R = me.Point('R')
   R.set_pos(Q, l*B.z + s*B.x)

   B_v_R = s.diff(t)*B.x
   B_v_R

.. jupyter-execute::

   r_T_S = R.pos_from(S)
   r_T_S

.. jupyter-execute::

   N_v_T = N_v_S + me.cross(B.ang_vel_in(N), r_T_S)
   N_v_T

.. jupyter-execute::

   N_v_R = B_v_R + N_v_T
   N_v_R

And with the :external:py:meth:`~sympy.physics.vector.point.Point.v1pt_theory`:

.. jupyter-execute::

   S.set_vel(B, 0)
   R.v1pt_theory(S, N, B)

Translational Acceleration
==========================

The acceleration of point :math:`P` in reference frame :math:`A` is defined as

.. math::

   {}^A\bar{a}^P := \frac{{}^A d {}^A\bar{v}^P}{dt}

:external:py:meth:`~sympy.physics.vector.point.Point.acc`:

.. jupyter-execute::

   S.acc(N)

Acceleration Two Point Thereom
==============================

.. math::

   {}^N\bar{v}^S = {}^A\bar{v}^P + {}^N\bar{\omega}^A \times \bar{r}^{S/P}

.. math::

   {}^N\bar{a}^S
   & = \frac{{}^N d\left({}^A\bar{v}^P\right)}{dt} +\frac{{}^N d \left( {}^N\bar{\omega}^A \times \bar{r}^{S/P}\right)}{dt} \\
   & = {}^N\bar{a}^P +
   \frac{{}^N d \left( {}^N\bar{\omega}^A \right)}{dt} \times \bar{r}^{S/P} +
   {}^N\bar{\omega}^A \times \frac{{}^N d  \left(\bar{r}^{S/P}\right)}{dt} \\
   & =
   {}^N\bar{a}^P +
   {}^N\bar{\alpha}^A \times\bar{r}^{S/P}
   {}^N\bar{\omega}^A\times\left({}^N\bar{\omega}^A \times\bar{r}^{S/P}\right)

Here we see clear the tangential component of acceleration:

.. math::

   {}^N\bar{\alpha}^A \times\bar{r}^{S/P}

.. jupyter-execute::

   me.cross(A.ang_acc_in(N), S.pos_from(P))

And the radial component:

.. math::

   {}^N\bar{\omega}^A\times\left({}^N\bar{\omega}^A \times\bar{r}^{S/P}\right)

.. jupyter-execute::

   me.cross(A.ang_vel_in(N), me.cross(A.ang_vel_in(N), S.pos_from(P)))

.. jupyter-execute::

   S.a2pt_theory(P, N, A)

Acceleration One Point Thereom
==============================

Starting with the expanded one point theorem for velocity:

.. math::

   {}^N\bar{v}^R = {}^B\bar{v}^R + {}^N\bar{v}^S + {}^N\bar{\omega}^B \times \bar{r}^{T/S}

and taking the time derivative in the frame :math:`N` the corollary formula for
acceleration can be derived:

.. math::

   {}^N\bar{a}^R
   & =
   \frac{{}^Nd {}^B\bar{v}^R}{dt} +
   \frac{{}^Nd {}^N\bar{v}^S}{dt} +
   \frac{{}^Nd {}^N\bar{\omega}^B \times \bar{r}^{T/S}}{dt} \\
   & =
   \frac{{}^Nd {}^N\bar{v}^R }{dt} +
   {}^N\bar{\omega}^B \times {}^N\bar{v}^R +
   {}^N\bar{a}^S +
   \frac{{}^Nd {}^N\bar{\omega}^B}{dt} \times \bar{r}^{T/S} +
   {}^N\bar{\omega}^B \times \frac{{}^Nd \bar{r}^{T/S}}{dt} \\
   & =
   {}^B\bar{a}^R +
   {}^N\bar{\omega}^B \times {}^B\bar{v}^R +
   {}^N\bar{a}^S +
   {}^N\bar{\alpha}^B \times \bar{r}^{T/S} +
   {}^N\bar{\omega}^B \times \left( {}^B\bar{v}^T +
   {}^N\bar{\omega}^B \times \bar{r}^{T/S} \right) \\
   & =
   {}^B\bar{a}^R +
   2{}^N\bar{\omega}^B \times {}^B\bar{v}^R +
   {}^N\bar{a}^S +
   {}^N\bar{\alpha}^B \times \bar{r}^{T/S} +
   {}^N\bar{\omega}^B \times \left(
   {}^N\bar{\omega}^B \times \bar{r}^{T/S} \right)

This is equivalent to

.. math::

   {}^N\bar{a}^R
   =
   {}^B\bar{a}^R +
   {}^N\bar{a}^T +
   2{}^N\bar{\omega}^B \times {}^B\bar{v}^R

:math:`2{}^N\bar{\omega}^B \times {}^N\bar{v}^R` is the Coriolis acceleration
that arises from :math:`R` moving in the rotating frame :math:`B`.

.. jupyter-execute::

   B_a_R = R.acc(B)
   B_a_R

.. jupyter-execute::

   N_a_T = R.a2pt_theory(S, N, B)
   N_a_T

.. jupyter-execute::

   2*me.cross(B.ang_vel_in(N), R.vel(B))

And with the :external:py:meth:`~sympy.physics.vector.point.Point.a1pt_theory`:

.. jupyter-execute::

   R.a1pt_theory(S, N, B)
