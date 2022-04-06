==================
Forces and Torques
==================

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`loads` or Jupyter Notebook:
   :jupyter-download:notebook:`loads`.

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   me.init_vprinting(use_latex='mathjax')

.. todo:: Define force and torque.

Bound and Free Vectors
======================

Vectors may have a *line of action* or not. A line of action is parallel to the
vector and passes through a particular point. If a vector has a line of action,
it is said to be *bound* to its line of action. If a vector is not bound to a
line of action it is said to be *free*.

Angular velocity is an example of a free vector. It has a direction and
magnitude, but is not associated with any line of action. A force vector, on
the other hand, is bound. If a force is applied to a rigid body, we must know
where on the body it is applied to resolve the force's effect. A force vector
acting on rigid body :math:`B` at point :math:`P` has a line of action through
:math:`P` and parallel to the force vector.

Moment
======

If a vector is a bound vector, then we can define its *moment* about a point.
The moment math:`\bar{M}` of bound vector :math:`\bar{v}` about point :math:`P`
is a vector itself and is defined as ([Kane1985]_, pg XX):

.. math::
   :label: eq-moment-definition

   \bar{M} := \bar{r}^{L/P} \times \bar{v}

:math:`\bar{r}^{L/P}` is a position vector from :math:`P` to any point
:math:`L` on the line of action of :math:`\bar{v}`.

A moment can be the result of a set of vectors. The *resultant* of a set
:math:`S` of vectors :math:`\bar{v}_1,\ldots,\bar{v}_\nu` is defined as:

.. math::
   :label: eq-resultant-definition

   \bar{R}^{S} := \sum{i=1}^{\nu} \bar{v}_i

If each vector in the resultant is bound, the sum of the moments due to each
vector about :math:`P` is call the moment of :math:`\bar{R}^{S}` about
:math:`P`.  This summation can be written as:

.. math::
   :label: eq-sum-moments

   \bar{M}^{S/P} = \sum{i=1}^{\nu} \bar{r}^{L_i/P} \times \bar{v}_i

The moment of bound vectors :math:`S` about one point :math:`P` is related to
the moment about another point :math:`Q` by ([Kane1985_], pg. 91):

.. math::
   :label: eq-moment-another-point

   \bar{M}^{S/P} = \bar{M}^{S/Q} + \bar{r}^{P/Q} \times \bar{R}^{S/Q}

:math:`\bar{R}^{S/Q}` is the resultant of the set :math:`S` bound to a line of
action through :math:`Q`.

.. todo:: Show a figure here to represent the left and right sides of the
   equation.

Couple
======

A set of bound vectors with a resultant equal to zero is called a *couple*. A
couple can have as many vectors as desired or needed with a minimum number
being two, such that :math:`\bar{R}^{S}=0`. A couple composed of two vectors is
called a *simple couple*.

.. todo:: add figure with three sketches of couples

The *torque* of a couple is the moment of the couple about a point. Because the
resultant of a couple is zero, the torque of a couple is the same about all
points. The torque, being a moment, is also a vector.

Equivalence & Replacement
=========================

Two sets of bound vectors are *equivalent* when they have these two properties:

1. equal resultants
2. equal moments about *any* point

If 1. and 2. are true, the sets are said to be *replacements* of each other.
Couples that have equal torques are equivalent, because the resultants are zero
and moments about any point are equal to the torque.

Given a set of bound vectors :math:`S` and a set of bound vectors that consist
of a torque of a couple :math:`\bar{T}` and vector :math:`\bar{v}` bound to an
arbitrary point :math:`P` it is a necessary and sufficient condition that the
second set is a replacement of the first if ([Kane1985]_, pg XX):

.. math::
   :label: eq-couple-torque-repl

   \bar{T} = \bar{M}^{S/P} \\
   \bar{v} = \bar{R}^S

This means that every set of bound vectors can be replaced by an equivalent
torque of a couple and a single bound vector that is the resultant of the
replaced set. This replacement simplifies the description of forces acting on
bodies.

Take for example the birds eye view of a four wheeled car which has front
steering and motors at each wheel allowing for precise control of the
propulsion forces at each wheel. A diagram of the forces acting at each wheel
is shown in Figure.

.. todo:: Creating car figure.

.. jupyter-execute::

   l, w = sm.symbols('l, w')
   Ffl, Ffr, Frl, Frr = me.dynamicsymbols('F_{fl}, F_{fr}, F_{rl}, F_{rr}')
   alphafl, alphafr = me.dynamicsymbols(r'\alpha_{fl}, \alpha_{fr}')
   alpharl, alpharr = me.dynamicsymbols(r'\alpha_{rl}, \alpha_{rr}')
   delta = me.dynamicsymbols('delta')

   B = me.ReferenceFrame('B')
   W = me.ReferenceFrame('W')
   FR = me.ReferenceFrame('FR')
   FL = me.ReferenceFrame('FL')
   RR = me.ReferenceFrame('RR')
   RL = me.ReferenceFrame('RL')

   W.orient_axis(B, delta, B.z)
   FR.orient_axis(W, alphafr, W.z)
   FL.orient_axis(W, alphafl, W.z)
   RR.orient_axis(B, alpharr, B.z)
   RL.orient_axis(B, alpharl, B.z)

The resultant of the forces is:

.. jupyter-execute::

   R = Ffl*FL.y + Ffr*FR.y + Frl*RL.y + Frr*RR.y
   R

This resultant is bound to a line of action through :math:`B_o`. The associated
couple is then:

.. jupyter-execute::

   T = (me.cross(l/2*B.y - w/2*B.x, Ffl*FL.y) +
        me.cross(l/2*B.y + w/2*B.x, Ffr*FR.y) +
        me.cross(-l/2*B.y - w/2*B.x, Frl*RL.y) +
        me.cross(-l/2*B.y + w/2*B.x, Frr*RR.y))
   T = T.simplify()
   T

Since we can always describe the forces acting on a rigid body as a resultant
force and an associate torque of a couple, we will take advantage of this
simpler form.

Specifying Forces and Torques
=============================

Forces are bound vectors that can be considered acting on specific points, thus
we will always need a vector and a point to fully describe the force. Methods
and functions in SymPy Mechanics that make use of forces will typically require
a tuple containing a point and a vector, for example the resultant force
:math:`R^{B/B_o}` acting on the mass center of of the car would be specified
like so:

.. jupyter-execute::

   Bo = me.Point('Bo')
   force = (Bo, R)
   force

Torques of a couple are free vectors (not bound to a line of action) but
represent the couple acting on a rigid body, thus a reference frame associated
with a rigid body and the vector representing the torque will be used to
describe the torque in SymPy Mechanics. For example:

.. jupyter-execute::

   torque = (B, T)
   torque

We will often refer to forces and torques collectively as *loads*.

Equal & Opposite
================

Both forces and torques applied to a multibody system must obey `Newton's Third
Law`_, i.e. that forces and torques act equal and opposite. Take for example a
torque from a motor that causes a pinned lever :math:`B` to rotate relative to
the ground :math:`N`. The motor torque occurs between the ground and the lever
(or more precisely the stator and the rotor which are fixed to the ground and
the lever). A sign convention must be chosen for the equal and opposite torque.

.. _Newton's Third Law: https://en.wikipedia.org/wiki/Newton's_laws_of_motion#Third_law

The motor torque can be specified as a time varying vector:

.. jupyter-execute::

   T, q = me.dynamicsymbols('T, q')

   N = me.ReferenceFrame('N')
   B = me.ReferenceFrame('B')
   B.orient_axis(N, q, N.z)

   Tm = T*N.z

Then the equal and opposite torques are captured by these two tuples:

.. jupyter-execute::

   (B, Tm), (N, -Tm)

with equal and opposite torques applied to each body.

.. warning::

   Careful about your sign convention. It is equally valid to choose `(B, -Tm),
   (N, Tm)`. But it is useful to choose a sign convention such that when the
   signs of angular velocity and torque are the same it corresponds to power
   into the system (from the motor in this case). So `B.orient_axis(N, q, N.z)`
   corresponds to `(T*N.z, B)` to power in with both are positive or both are
   negative. This is just a convention though and the choice of force and
   torque signs can be anything, just make sure you know and understand what it
   is!

Contributing and Non-contributing
=================================

*Contributing forces and torques* are those that do work on the multibody
system. Work is defined as:

.. math::
   :label: eq-work-definition

   W = \int \bar{F} \cdot d\bar{x}

The gravitational force acting on a particle moving through a unidirectional
constant gravitational field does work on the system.

*Non-contributing forces and torques* do no work on the system. For example,
when a force acts between two points that have no relative motion, no work is
done. The contact forces between two rigid bodies if the two bodies are
connected at a single point.

Gravity
=======

We will often be interested in a multibody systems motion when it is subject to
gravitational forces. The simplest case is a constant unidirectional
gravitional field, which is appropriate model for small objects moving on and
near the Earth's surface. The gravitational forces can be applied soley to the
mass centers of each rigid body in a multibody system as a resultant force. The
gravitional torque on the bodies is zero because the force is equal in
magnitude for each particle in the body. See [Kane1985]_ pg. XX for the more
general case of Newton's Law of Gravitation where this is not the case which
often comes into play for modeling spacecraft.

In SymPy Mechanics a gravitational force acting on a particle of mass :math:`m`
with acceleration due to gravity being :math:`g` in the :math:`-\hat{n}_y`
direction would take this form:

.. jupyter-execute::

   m, g = sm.symbols('m, g')
   Fg = -m*g*N.y

Springs & Dampers
=================

Idealized springs and dampers are useful models of elements that have distance
and velocity dependent forces and torques. A spring with free length
:math:`q_0` and where :math:`q_1,q_2` locate the ends of the spring along a
line parallel to the :math:`\hat{n}_x` direction taking a sign convention that
a positive spring force acting on the :math:`q_2` end of the spring is in the
negative :math:`\hat{n}_x` direction . If the spring is linear with stiffness
:math:`k` the spring force vector is then:

.. jupyter-execute::

   q0, k = sm.symbols('q0, k')
   q1, q2 = me.dynamicsymbols('q1, q2')

   displacement = q2 - q1 - q0

   Fs = -k*displacement*N.x
   Fs

.. todo:: Add figure of spring and damper with force directions.

Similarly, a linear damping force with damping coefficient :math:`c` is defined
as:

.. jupyter-execute::

   c = sm.symbols('c')
   t = me.dynamicsymbols._t

   Fc = -c*displacement.diff(t)*N.x
   Fc

Aerodynamic
===========

Aerodynamic drag of a blunt body is dominated by the frontal area drag and the
magnitude of this drag force can be modeled with the following equation:

.. math::
   :label: eq-aerodynamic-drag

   \frac{1}{2}\rhoC_dAv^2

where :math:`\rho` is the density of the air, :math:`C_d` is the drag
coefficient, :math:`A` is the frontal area, and :math:`v` is the air speed
relative to the body.

If a body is moving in still air at an aribtrary velocity and point :math:`P`
is the aerodynamic center of the body then the aerodynamic drag force vector
that opposes the motion can be found with such an equation:

.. jupyter-execute::

   A, Cd, rho = sm.symbols('A, C_d, rho')
   ux, uy, uz = me.dynamicsymbols('u_x, u_y, u_z')

   N_v_P = ux*N.x + uy*N.y + uz*N.z

   Fd = -N_v_P.normalize()*Cd*A*rho/2*N_v_P.dot(N_v_P)
   Fd

If the motion is only along the :math:`\hat{n}_x` direction, for example, the
equation for the drag force vector reduces to:

.. jupyter-execute::

   Fd.xreplace({uy: 0, uz:0})

.. todo:: This may be incorrect should have Abs(ux).

Friction
========

Coulomb's Law is the simplest model of friction constant friction.

.. jupyter-execute::

   mu, m, g = sm.symbols('mu, m, g')

   Fn = m*g

   Ff = sm.Piecewise((mu*Fn, displacement.diff(t) > 0),
                     (-mu*Fn, displacement.diff(t) < 0),
                     (0, True))*N.x
   Ff

.. jupyter-execute::

   Ff = mu*Fn*sm.sign(displacement.diff(t))*N.x
   Ff


Collision
=========

If two points, a point and a surface, or two surfaces collide the impact
behavior depends on the material properties and mass of the colliding bodies. A
simple way to model impact is to create a stiff spring that only engages if one
body pentrates the other body.

.. jupyter-execute::

   x, y, z = me.dynamicsymbols('x, y, z')

   r_O_P = x*N.x + y*N.y + z*N.z

   penetration = r_O_P.dot(N.z)

   Fc = sm.Piecewise((-k*penetration, penetration < 0), (0, True))
   Fc

.. jupyter-execute::

   Fc = sm.Abs(penetration) / penetration
   Fc

.. todo:: Add a model with some damping in the plane direction.
