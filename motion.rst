========================
Nonholonomic Constraints
========================

.. note::

   You can download this example as a Python script:
   :jupyter-download-script:`motion` or Jupyter Notebook:
   :jupyter-download-notebook:`motion`.

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   me.init_vprinting(use_latex='mathjax')

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

- Formulate nonholonomic constraints to constrain motion.
- Determine if a nonholonomic constraint is essential and not simply the time
  derivative of a  holonomic constraint.
- Formulate a motion constraint for rolling without slip.
- Define kinematical differential equations and solve them to put in first
  order form.
- Select different choices of generalized speeds.
- Solve for the dependent generalized speeds in terms of the independent
  generalized speeds.
- Calculate the degrees of freedom of a multibody system.

Motion Constraints
==================

In :ref:`Holonomic Constraints`, we discussed constraints on the configuration
of a system. When defining configuration we are only concerned with the
locations of points and how reference frames are oriented. In this chapter, we
will consider constraints on the motion of a system. Motion concerns how points
and reference frames move in time. Take parallel parking a car as a motivating
example, :numref:`fig-motion-parallel`.

.. _fig-motion-parallel:
.. figure:: figures/motion-parallel.svg
   :align: center

   a) two positions (or configurations) of car 2 relative to cars 1 and 3, b)
   simplest motion to move car 2 into an empy spot between cars 1 and 3, c)
   actual motion to move car 2 into the empty spot

We know that car 2 can be in either the left or right location in a), i.e. the
car's configuration permits either location. But the motion scenario in b) is
not possible. A car cannot move from the left configuration to the right
configuration by simply sliding directly to the right (see the note below if
you question this). Although, this surely would be nice if we could. A car has
wheels and only the front wheels can be steered, so the scenario in c) is a
viable motion for the car to end up in the correct final configuration. The car
has to *move* in a specific way to get from one configuration to another. This
implies that we have some kind of constraint on the motion but not the
configuration. Constraints such as these are called *nonholonomic constraints*
and they take the form:

.. math::
   :label: eq-nonholonomic-qdot

   \bar{f}_n(\dot{\bar{q}}, \bar{q}, t) = 0 \\
   \textrm{ where } \\
   \bar{f}_n \in \mathbb{R}^m \\
   \bar{q} = \left[ q_1, \ldots, q_n\right]^T \in \mathbb{R}^n

The :math:`m` constraints involve the time derivatives of the generalized
coordinates and arise from scalar equations derived from velocities.

.. note::

   We could find a very strong person to push the car sideways, overcoming the
   very high resisting friction force. It is important to note that any
   constraint is just a model of a physical phenomena. We know that if we push
   hard enough and low enough that the car's lateral motion is not constrained.
   Also, if the car were on ice, then the nonholomonic constraint may be a poor
   modeling decision.

Chaplygin Sleigh
================

Take the simple example of the `Chaplygin Sleigh`_, sketched out in
:numref:`fig-motion-sleigh`. A sleigh can slide along a flat plane, but can
only move in the direction it is pointing, much like the wheels of the car
above. This system is described by three generalized coordinates
:math:`x,y,\theta`. For the motion to only occur along its body fixed
:math:`\hat{a}_x` direction, the component of velocity in the body fixed
:math:`\hat{a}_y` direction must equal zero at all times.

.. _Chaplygin Sleigh: https://en.wikipedia.org/wiki/Chaplygin_sleigh

.. _fig-motion-sleigh:
.. figure:: figures/motion-sleigh.svg
   :align: center

   Configuration diagram of a Chaplygin Sleigh. The rectangle :math:`A`
   represents a sleigh moving on a plane. Point :math:`P` represents the center
   of the sleigh.

Using SymPy Mechanics we can find the velocity of :math:`P` and express it in
the :math:`A` reference frame:

.. jupyter-execute::

   x, y, theta = me.dynamicsymbols('x, y, theta')

   N = me.ReferenceFrame('N')
   A = me.ReferenceFrame('A')

   A.orient_axis(N, theta, N.z)

   O = me.Point('O')
   P = me.Point('P')

   P.set_pos(O, x*N.x + y*N.y)

   O.set_vel(N, 0)

   P.vel(N).express(A)

The single scalar nonholonomic constraint then takes this form:

.. math::
   :label: eq-chaplygin-sleigh-constraint

   {}^N\bar{v}^P \cdot \hat{a}_y = 0

because there can be no velocity component in the :math:`\hat{a}_y` direction.
With SymPy, this is:

.. jupyter-execute::

   fn = P.vel(N).dot(A.y)
   fn

How do we know that this is, in fact, a nonholonomic constraint and not simply
the time derivative of a holonomic constraint?

Recall one of the four-bar linkage holonomic constraints arising from Eq.
:math:numref:`constraint-expression` and time differentiate it:

.. jupyter-execute::

   t = me.dynamicsymbols._t

   q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')
   la, lb, lc, ln = sm.symbols('l_a, l_b, l_c, l_n')

   fhx = la*sm.cos(q1) + lb*sm.cos(q1 + q2) + lc*sm.cos(q1 + q2 + q3) - ln
   sm.trigsimp(fhx.diff(t))

This looks like a nonholonomic constraint, i.e. it has time derivatives of the
coordinates, but we know that if we integrate this equation with respect to
time we can retrieve the original holonomic constraint, so it really isn't a
nonholonomic constraint even though it looks like one.

So if we can integrate :math:`f_n` with respect to time and we arrive at a
function of only the generalized coordinates and time, then we do not have a
nonholonomic constraint, but a holonomic constraint in disguise. Unfortunately,
it is not generally possible to integrate :math:`f_n` so we must check the
integrability of :math:`f_n` indirectly.

If :math:`f_n` of the sleigh was the time derivative of a holonomic constraint
:math:`f_h` then it must be able to be expressed in this form:

.. math::
   :label: eq-diff-holonomic

   f_n = \frac{d f_h}{dt} =
   \frac{\partial f_h}{\partial x} \frac{dx}{dt} +
   \frac{\partial f_h}{\partial y} \frac{dy}{dt} +
   \frac{\partial f_h}{\partial \theta} \frac{d\theta}{dt} +
   \frac{\partial f_h}{\partial t}

and a `condition of integrability is that the mixed partial derivatives must
commute <https://en.wikipedia.org/wiki/Symmetry_of_second_derivatives>`_. By
inspection of :math:`f_n` we see that we can extract the partial derivatives by
collecting the coefficients. SymPy's
:external:py:meth:`~sympy.core.expr.Expr.coeff` can extract the linear
coefficients for us:

.. jupyter-execute::

   dfdx = fn.coeff(x.diff(t))
   dfdy = fn.coeff(y.diff(t))
   dfdth = fn.coeff(theta.diff(t))

   dfdx, dfdy, dfdth

Each pair of mixed partials can be calculated. For example
:math:`\frac{\partial^2 f_h}{\partial y \partial x}` and
:math:`\frac{\partial^2 f_h}{\partial x \partial y}`:

.. jupyter-execute::

   dfdx.diff(y), dfdy.diff(x)

and the other two pairs:

.. jupyter-execute::

   dfdx.diff(theta), dfdth.diff(x)

.. jupyter-execute::

   dfdy.diff(theta), dfdth.diff(y)

We see that to for the last two pairs, the mixed partials do not commute. This
proves that :math:`f_n` is not integrable and is thus an essential nonholonomic
constraint that is not a holonomic constraint in disguise.

.. admonition:: Exercise

   Check whether the mixed partials of the time derivative of the four-bar
   linkage constraints commute.

.. admonition:: Solution
   :class: dropdown

   .. jupyter-execute::

      fnx = fhx.diff(t)
      dfdq1 = fnx.diff(q1)
      dfdq2 = fnx.diff(q2)
      dfdq3 = fnx.diff(q3)

   All of the mixed partials are the same:

   .. jupyter-execute::

      dfdq1.diff(q2) -  dfdq2.diff(q1)

   .. jupyter-execute::

      dfdq2.diff(q3) - dfdq3.diff(q2)

   .. jupyter-execute::

      dfdq3.diff(q1) - dfdq1.diff(q3)

   All of the mixed partials are the same so this is a holonomic constraint in
   disguise.

Rolling Without Slip
====================

It is quite common to make the modeling assumption that a wheel rolls without
slip. A wheel best provides its beneficial properties of rolling and propulsion
by ensuring that the friction between the wheel and the surface it rolls on is
sufficiently high. This avoids relative motion between a point fixed on the
wheel and a point fixed on the surface located at the wheel-surface contact
location at any given time. This nature can be modeled by a motion constraint.
The key to developing the constraint to ensure there is no relative slip
velocity is to identify the correct two points, calculate the velocity of those
points, and specify that the relative velocity is zero.

.. _fig-motion-wheel:
.. figure:: figures/motion-wheel.svg
   :align: center

   A 2D disc :math:`B` rolling on a motionless plane :math:`N`.

For example, when a 2D disc :math:`B` rolls without slip over a motionless
plane :math:`N` (:numref:`fig-motion-wheel`), the velocity of a point :math:`C`
fixed in :math:`B` at the contact point with the plane must be zero to ensure
no slip when observed from the plane's reference frame. We can state this
mathematically as:

.. math::

   {}^N\bar{v}^{C} = 0

One must be careful about calculating this velocity and recognizing that there
are numerous points of possible interest at the same wheel-plane contact
location. You may consider these points, for example:

- A point :math:`B_C` that moves in the plane :math:`N` which is *always*
  located at the wheel-plane contact location. The coordinate :math:`q_1`
  tracks this point in the figure.
- A point :math:`G_C` that is fixed in the wheel which follows a cycloid_ curve
  as it rolls along.
- A point :math:`G` that is fixed in the plane which is located at the
  wheel-plane contact point at any given instance of time.
- A point :math:`C` that is fixed in the wheel which is located at the
  wheel-plane contact point at any given instance of time.

.. _cycloid: https://en.wikipedia.org/wiki/Cycloid

A motion constraint that ensures rolling without slip, can only be formed by
considering the last two points. The vector constraint equation is:

.. math::

   {}^N\bar{v}^{C} - {}^N\bar{v}^{G} = 0

Point :math:`G` is fixed in :math:`N` so it has no velocity in :math:`N`:

.. math::

   {}^N\bar{v}^{G} = 0

Point :math:`C` is fixed in :math:`B`. To determine its velocity, take
:math:`B_o` to be the wheel center which is also fixed in :math:`B`. Since both
points are fixed in :math:`B` we can apply the two point velocity theorem.

.. math::

   {}^N\bar{v}^{C} = {}^N\bar{v}^{B_o} + {}^N\bar{\omega}^B \times \bar{r}^{C/B_o}

We can then use two generalized coordinates to describe the position
:math:`q_1` (from :math:`O` fixed in :math:`N`) and rotation :math:`q_2` of the
wheel. The velocity of the wheel center is then:

.. math::

   {}^N\bar{v}^{B_o} = \dot{q}_1\hat{n}_x

The cross product terms are found with the radius of the wheel with :math:`r`
and the angular velocity to give the velocity of :math:`C`:

.. math::
   :label: eq-roll-slip

   {}^N\bar{v}^{C} = & \dot{q}_1\hat{n}_x + \dot{q}_2 \hat{n}_z \times -r\hat{n}_y \\
   {}^N\bar{v}^{C} = & \dot{q}_1\hat{n}_x + \dot{q}_2 r \hat{n}_x

Applying the motion constraint and knowing that :math:`{}^N\bar{v}^{G} = 0`
gives us this scalar constraint equation directly from
:math:numref:`eq-roll-slip`:

.. math::

   \dot{q}_1 + \dot{q}_2 r = 0

This is a scalar constraint equation that ensures rolling without slip and
involves the time derivatives of the coordinates. It is integrable and thus
actually a holonomic constraint, i.e. :math:`q_1 + q_2 r = 0`. General rolling
without slip in three dimensions will be nonholonomic. Take care to calculate
the relative velocities of the two points fixed in each of the bodies in
rolling contact that are located at the contact point at that *instance of
time*.

.. todo:: Exercise to calculate the constraint if the plane has a horizontal
   velocity.

Kinematical Differential Equations
==================================

In Eq. :math:numref:`eq-nonholonomic-qdot` we show the form of the nonholonomic
constraints in terms of :math:`\dot{\bar{q}}`. Newton's and Euler's Second Laws
of motion will require calculation of acceleration and angular acceleration
respectively. These laws of motion are second order differential equations
because it involves second time derivatives of distances and angles. Any second
order differential equation can be equivalently represented by two first order
differential equations by introducing a new variable for any first derivative
terms. We are working towards writing the equations of motion of a multibody
system, which will be differential equations that are most useful for
simulation when in a first order form. To do this, we now introduce the
variables :math:`\bar{u} = \left[u_1, \ldots, u_n\right]^T` and define them as
linear functions of the time derivatives of the generalized coordinates
:math:`\dot{q}_1, \ldots, \dot{q}_n`. These variables are called *generalized
speeds*. They take the form:

.. math::
   :label: eq-generalized-speeds

   \bar{u} := \mathbf{Y}_k(\bar{q}, t) \dot{\bar{q}} + \bar{z}_k(\bar{q}, t)

:math:`\bar{u}` must be chosen such that :math:`\mathbf{Y}_k` is invertible. If
it is, then we solve for :math:`\dot{\bar{q}}` we can write these first order
differential equations as such:

.. math::
   :label: eq-kinematical-diff-eq

   \dot{\bar{q}} = \mathbf{Y}_k^{-1}\left(\bar{u} - \bar{z}_k\right)

Eq. :math:numref:`eq-kinematical-diff-eq` are called the *kinematical
differential equations*.

The most common, and always valid, choice of generalized speeds is:

.. math::
   :label: eq-simplest-generalized-speeds

   \bar{u} = \mathbf{I} \dot{\bar{q}}

where :math:`\mathbf{I}` is the identity matrix. This results in :math:`u_i =
\dot{q}_i` for :math:`i=1,\ldots,n`.

Now that we have introduced generalized speeds, the nonholonomic constraints
can then be written as:

.. math::
   :label: nonholonomic-constraints-u

   \bar{f}_n(\bar{u}, \bar{q}, t) = 0 \\
   \textrm{ where } \\
   \bar{f}_n \in \mathbb{R}^m \\
   \bar{u} = \left[ u_1, \ldots, u_n\right]^T \in \mathbb{R}^n\\
   \bar{q} = \left[ q_1, \ldots, q_n\right]^T \in \mathbb{R}^n

Choosing Generalized Speeds
===========================

There are many possible choices for generalized speed and you are free to
select them as you please, as long as they fit the form of equation
:math:numref:`eq-generalized-speeds` and :math:`\mathbf{Y}_k` is invertible.
Some selections of generalized speeds can reduce the complexity of important
velocity expressions and if selected carefully may reduce the complexity of the
equations of motion we will derive in a later chapters (see [Mitiguy1996]_ for
examples). To see some examples of selecting generalized speeds, take for
example the angular velocity of a reference frame which is oriented with a
:math:`z\textrm{-}x\textrm{-}y` body fixed orientation:

.. jupyter-execute::

   q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')

   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')

   B.orient_body_fixed(A, (q1, q2, q3), 'ZXY')

   A_w_B = B.ang_vel_in(A).simplify()
   A_w_B

Choice 1
--------

If we choose the simplest definition for the :math:`u`'s, i.e.
:math:`u_1=\dot{q}_1`, :math:`u_2=\dot{q}_2`, and :math:`u_3=\dot{q}_3`, the
angular velocity takes this form:

.. jupyter-execute::

   u1, u2, u3 = me.dynamicsymbols('u1, u2, u3')

   t = me.dynamicsymbols._t
   qdot = sm.Matrix([q1.diff(t), q2.diff(t), q3.diff(t)])
   u = sm.Matrix([u1, u2, u3])

   A_w_B = A_w_B.xreplace(dict(zip(qdot, u)))
   A_w_B

.. jupyter-execute::

   Yk_plus_zk = qdot
   Yk_plus_zk

Recall from :ref:`Solving Linear Systems` that the Jacobian is a simple way
to extract the coefficients of linear terms into a coefficient matrix for a
system of linear equations. In this case, we see that this results in the
identity matrix.

.. jupyter-execute::

   Yk = Yk_plus_zk.jacobian(qdot)
   Yk

Now find :math:`\bar{z}_k` by setting the time derivatives of the generalized
coordinates to zero:

.. jupyter-execute::

   qd_zero_repl = dict(zip(qdot, sm.zeros(3, 1)))
   qd_zero_repl

.. jupyter-execute::

   zk = Yk_plus_zk.xreplace(qd_zero_repl)
   zk

The linear equation can be solved for the :math:`\dot{q}`'s, (Eq.
:math:numref:`eq-kinematical-diff-eq`):

.. jupyter-execute::

   sm.Eq(qdot, Yk.LUsolve(u - zk))

Choice 2
--------

Another valid choice is to set the :math:`u`'s equal to each measure number of
the angular velocity expressed in :math:`B`:

.. math::
   :label: u-choice-2

   u_1 = {}^A\bar{\omega}^B \cdot \hat{b}_x \\
   u_2 = {}^A\bar{\omega}^B \cdot \hat{b}_y \\
   u_3 = {}^A\bar{\omega}^B \cdot \hat{b}_z

so that:

.. math::
   :label: omega-choice-2

   {}^A\bar{\omega}^B = u_1\hat{b}_x + u_2\hat{b}_y + u_3\hat{b}_z

.. jupyter-execute::

   A_w_B = B.ang_vel_in(A).simplify()
   A_w_B

.. jupyter-execute::

   u1_expr = A_w_B.dot(B.x)
   u2_expr = A_w_B.dot(B.y)
   u3_expr = A_w_B.dot(B.z)

   Yk_plus_zk = sm.Matrix([u1_expr, u2_expr, u3_expr])
   Yk_plus_zk

.. jupyter-execute::

   Yk = Yk_plus_zk.jacobian(qdot)
   Yk

.. jupyter-execute::

   zk = Yk_plus_zk.xreplace(qd_zero_repl)
   zk

Now we form:

.. jupyter-execute::

   sm.Eq(qdot, sm.trigsimp(Yk.LUsolve(u - zk)))

.. note::

   Notice how the kinematical differential equations are not valid when
   :math:`q_2` or :math:`q_3` are even multiples of :math:`\pi/2`. If your
   system must orient through these values, you'll need to select a different
   body fixed rotation or an orientation method that isn't suseptible to these
   issues.

Choice 3
--------

Another valid choice is to set the :math:`u`'s equal to each measure number of
the angular velocity expressed in :math:`A`:

.. math::
   :label: u-choice-3

   u_1 = {}^A\bar{\omega}^B \cdot \hat{a}_x \\
   u_2 = {}^A\bar{\omega}^B \cdot \hat{a}_y \\
   u_3 = {}^A\bar{\omega}^B \cdot \hat{a}_z

so that:

.. math::
   :label: omega-choice-3

   {}^A\bar{\omega}^B = u_1\hat{a}_x + u_2\hat{a}_y + u_3\hat{a}_z

.. jupyter-execute::

   A_w_B = B.ang_vel_in(A).express(A).simplify()
   A_w_B

.. jupyter-execute::

   u1_expr = A_w_B.dot(A.x)
   u2_expr = A_w_B.dot(A.y)
   u3_expr = A_w_B.dot(A.z)

   Yk_plus_zk = sm.Matrix([u1_expr, u2_expr, u3_expr])
   Yk_plus_zk

.. jupyter-execute::

   Yk = Yk_plus_zk.jacobian(qdot)
   Yk

.. jupyter-execute::

   zk = Yk_plus_zk.xreplace(qd_zero_repl)
   zk

.. jupyter-execute::

   sm.Eq(qdot, sm.trigsimp(Yk.LUsolve(u - zk)))

Snakeboard
==========

A snakeboard_ is a variation on a skateboard that can be propelled via
nonholonomic locomotion [Ostrowski1994]_. Similar to the Chaplygin Sleigh, the
wheels can generally only travel in the direction they are pointed. This
classic video from 1993 shows how to propel the board:

.. raw:: html

   <center>
   <iframe width="560" height="315"
   src="https://www.youtube.com/embed/yxlC95YjmEs" title="YouTube video player"
   frameborder="0" allow="accelerometer; autoplay; clipboard-write;
   encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
   </center>

.. _snakeboard: https://en.wikipedia.org/wiki/Snakeboard

:numref:`fig-snakeboard` shows what a real Snakeboard looks like and
:numref:`fig-snakeboard-configuration` shows a configuration diagram.

.. _fig-snakeboard:
.. figure:: https://upload.wikimedia.org/wikipedia/commons/thumb/6/61/Snakeboard_down.jpg/640px-Snakeboard_down.jpg
   :align: center

   Example of a snakeboard that shows the two footpads each with attached truck
   and pair of wheels that are connected by the coupler.

   Николайков Вячеслав, `CC BY-SA 3.0
   <https://creativecommons.org/licenses/by-sa/3.0>`_, via Wikimedia Commons

.. _fig-snakeboard-configuration:
.. figure:: figures/motion-snakeboard.svg
   :align: center

   Configuration diagram of a planar Snakeboard model.

Start by defining the time varying variables and constants:

.. jupyter-execute::

   q1, q2, q3, q4, q5 = me.dynamicsymbols('q1, q2, q3, q4, q5')
   l = sm.symbols('l')

The reference frames are all simple rotations about the axis normal to the
plane:

.. jupyter-execute::

   N = me.ReferenceFrame('N')
   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')
   C = me.ReferenceFrame('C')

   A.orient_axis(N, q3, N.z)
   B.orient_axis(A, q4, A.z)
   C.orient_axis(A, q5, A.z)

The angular velocities of each reference frame are then:

.. jupyter-execute::

   A.ang_vel_in(N)

.. jupyter-execute::

   B.ang_vel_in(N)

.. jupyter-execute::

   C.ang_vel_in(N)

Establish the position vectors among the points:

.. jupyter-execute::

   O = me.Point('O')
   Ao = me.Point('A_o')
   Bo = me.Point('B_o')
   Co = me.Point('C_o')

   Ao.set_pos(O, q1*N.x + q2*N.y)
   Bo.set_pos(Ao, l/2*A.x)
   Co.set_pos(Ao, -l/2*A.x)

The velocity of :math:`A_o` in :math:`N` is a simple time derivative:

.. jupyter-execute::

   O.set_vel(N, 0)
   Ao.vel(N)

The two point theorem is handy for computing the other two velocities:

.. jupyter-execute::

   Bo.v2pt_theory(Ao, N, A)

.. jupyter-execute::

   Co.v2pt_theory(Ao, N, A)

The unit vectors of :math:`B` and :math:`C` are aligned with the wheels of the
Snakeboard. This lets us impose that there is no velocity in the direction
normal to the wheel's rolling direction by taking dot products with the
respectively reference frames' :math:`y` direction unit vector to form the two
nonholonomic constraints:

.. math::
   :label: eq-snakeboard-constraints

   {}^N\bar{v}^{Bo} \cdot \hat{b}_y = 0 \\
   {}^N\bar{v}^{Co} \cdot \hat{c}_y = 0

.. jupyter-execute::

   fn = sm.Matrix([Bo.vel(N).dot(B.y),
                   Co.vel(N).dot(C.y)])
   fn = sm.trigsimp(fn)
   fn

Now we introduce some generalized speeds. By inspection of the above constraint
equations, we can see that defining a generalized speed equal to
:math:`\frac{l\dot{q}_3}{2}` can simplify the equations a bit. So define these
generalized speeds:

.. math::
   :label: eq-snakeboard-gen-speeds

   u_i = \dot{q}_i \textrm{ for } i=1,2,4,5 \\
   u_3 = \frac{l\dot{q}_3}{2}

Now replace all of the time derivatives of the generalized coordinates with the
generalized speeds. We use :external:py:meth:`~sympy.core.basic.Basic.subs`
here because the replacement isn't an exact replacement (in the sense of
:external:py:meth:`~sympy.core.basic.Basic.xreplace`).

.. jupyter-execute::

   u1, u2, u3, u4, u5 = me.dynamicsymbols('u1, u2, u3, u4, u5')

   u_repl = {
       q1.diff(): u1,
       q2.diff(): u2,
       l*q3.diff()/2: u3,
       q4.diff(): u4,
       q5.diff(): u5
   }

   fn = fn.subs(u_repl)
   fn

These nonholonomic constraints take this form:

.. math::
   :label: eq-general-con

   \bar{f}_n(u_1, u_2, u_3, q_3, q_4, q_5) = 0 \textrm{ where } \bar{f}_n \in \mathbb{R}^2

We now have two equations with three unknown generalized speeds. Note that all
of the generalized coordinates are not present in these constraints which is
common. We can solve for two of the generalized speeds in terms of the third.
So we select two as dependent generalized speeds and one as an independent
generalized speed.  Because nonholonomic constraints are derived from measure
numbers of velocity vectors, the nonholonomic constraints are always linear in
the generalized speeds. If we introduce :math:`\bar{u}_s` as a vector of
independent generalized speeds and :math:`\bar{u}_r` as a vector of dependent
generalized speeds, the nonholonomic constraints can be written as:

.. math::
   :label: eq-constraint-linear-form

   \bar{f}_n(\bar{u}_s, \bar{u}_r, \bar{q}, t) =
   \mathbf{A}_r(\bar{q}, t) \bar{u}_r +
   \mathbf{A}_s(\bar{q}, t) \bar{u}_s +
   \bar{b}_{rs}(\bar{q}, t) = 0

or

.. math::
   :label: eq-constraint-linear-form-solve

   \bar{u}_r = \mathbf{A}_r^{-1}\left(-\mathbf{A}_s \bar{u}_s - \bar{b}_{rs}\right) \\
   \bar{u}_r = \mathbf{A}_n \bar{u}_s + \bar{b}_n

For the Snakeboard let's choose :math:`\bar{u}_s = [u_3, u_4, u_5]^T` as the
independent generalized speeds and :math:`\bar{u}_r = [u_1, u_2]^T` as the
dependent generalized speeds.

.. jupyter-execute::

   us = sm.Matrix([u3, u4, u5])
   ur = sm.Matrix([u1, u2])

:math:`\mathbf{A}_r` are the linear coefficients of :math:`\bar{u}_r` so:

.. jupyter-execute::

   Ar = fn.jacobian(ur)
   Ar

:math:`\mathbf{A}_s` are the linear coefficients of :math:`\bar{u}_s` so:

.. jupyter-execute::

   As = fn.jacobian(us)
   As

:math:`\bar{b}_{rs}` remains when :math:`\bar{u}=0`:

.. jupyter-execute::

   brs = fn.xreplace(dict(zip([u1, u2, u3, u4, u5], [0, 0, 0, 0, 0])))
   brs

:math:`\mathbf{A}_n` and :math:`\bar{b}_n` are formed by solving the linear
system:

.. jupyter-execute::

   An = Ar.LUsolve(-As)
   An = sm.simplify(An)
   An

.. jupyter-execute::

   bn = Ar.LUsolve(-brs)
   bn

We now have the :math:`m=2` dependent generalized speeds
:math:`\bar{u}_r=\left[u_1,u_2\right]^T` written as functions of the `n=1`
independent generalized speeds :math:`\bar{u}_s=\left[u_3\right]`:

.. jupyter-execute::

   sm.Eq(ur, An*us + bn)

Degrees of Freedom
==================

For simple nonholonomic systems observed in a reference frame :math:`A`, such
as the Chaplygin Sleigh or the Snakeboard, the *degrees of freedom* in
:math:`A` are equal to the number of independent generalized speeds. The number
of degrees of freedom :math:`p` is defined as:

.. math::
   :label: eq-dof-definition

   p := n - m

where :math:`n` is the number of generalized coordinates and :math:`m` are the
number of nonholonomic constraints (and thus dependent generalized speeds). If
there are no nonholonomic constraints, the system is a holonomic system in
:math:`A` and :math:`p=n` making the number of degrees of freedom equal to the
number of generalized coordinates.

.. admonition:: Exercise

   What are the number of degrees of freedom for the Chaplygin Sleigh,
   Snakeboard, and Four-bar linkage?

.. admonition:: Solution
   :class: dropdown

   The Chapylgin Sleigh has :math:`n=3` generalized coordinates
   :math:`x,y,\theta` and :math:`m=1` nonholonomic constraints. The degrees of
   freedom are then :math:`p = 3 - 1 = 2`.

   The Snakeboard has :math:`n=5` generalized coordinates and :math:`m=2`
   nonholonomic constraints. The degrees of freedom are then :math:`p = 5 - 2 =
   3`.

   We described the four-bar linkage with :math:`N=3` coordinates and there
   were :math:`M=2` holonomic constraints leaving us with :math:`n=N-M=3-2=1`
   generalized coordinates. There are no nonholonomic constraints so
   :math:`m=0`. This means that there :math:`p=n-m=1-0=1` degrees of freedom.

It is not always easy to visualize the degrees of freedom of a nonholonomic
system when thinking of its motion, but for holonomic systems thought
experiments where you vary one or two generalized coordinates at a time can
help you visualize the motion.

If you have a holonomic system (no nonholonomic constraints) the degrees of
freedom are equal to the number of generalized coordinates. Nonholonomic
systems (those with non-integrable motion constraints) have fewer degrees of
freedom than the number of generalized coordinates.

