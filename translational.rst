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

Two Point Thereom
=================

One Point Thereom
=================

Translational Acceleration
==========================

Two Point Thereom
=================

One Point Thereom
=================
