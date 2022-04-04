==========================
Forces and Torques (Loads)
==========================

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`loads` or Jupyter Notebook:
   :jupyter-download:notebook:`loads`.

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   me.init_vprinting(use_latex='mathjax')

Bound and Free Vectors
======================

Vectors may have a *line of action* or not. If a vector has a line of action,
it is said to be *bound* to its line of action :math:`L`. If a vector is not
bound to a line of action it is said to be *free*.

Angular velocity is an example of a free vector. It has a direction and
magnitude, but is not tied to any line of action. A force vector, on the other
hand, is bound. If a force is applied to a rigid body, we must know where on
the body it is applied to resolve the force's effect. A force vector acting on
rigid body :math:`B` at point :math:`P` has a line of action through :math:`P`
and parallel to the force.

Moment
======

If a vector is a bound vector, then we can define its *moment* about a point.
The moment of :math:`\bar{M}` of bound vector :math:`\bar{v}` about point
:math:`P` is then defined as:

.. math::
   :label: eq-moment-definition

   \bar{M} := \bar{r}^{P_L/P} \times \bar{v}

:math:`\bar{r}^{P_L/P}` is a position vector from :math:`P` to any point on the
line of action of :math:`\bar{v}`.

A moment can be the result of a set of vectors acting about a point. The
*resultant* of a set of vectors :math:`\bar{v}_1,\ldots,\bar{v}_n` is defined
as:

.. math::
   :label: eq-resultant-definition

   \bar{R} := \sum{i=1}^{n} \bar{v}_i

If each vector in the resultant is bound, the sum of the moments about
:math:`P` is call the moment of :math:`\bar{R}` about :math:`P`.

Couple
======

A set of bound vectors with a resultant equal to zero is called a *couple*.  A
couple can have as many vectors as desired or needed with a minimum number
being two, such that :math:`\bar{R}=0`. A couple composed of two vectors is
called a *simple couple*.

.. todo:: add figure with three sketches of couples

The *torque* of a couple is the moment of the couple about a point. Because the
resultant of a couple is zero, the torque of a couple is the same about all
points.

Equivalence Replacement
=======================

Two sets of bound vectors are equivalent when they have these two properties:

1. equal resultants
2. equal moments about any point

If 1 and 2 are true, the sets are said to be replacements of each other.
Couples that have equal torques are equivalent, since the resultants are zero
and moments about any point is equal to the value of the torque.

The moment about one point :math:`P` is related to the moment about another
point :math:`Q` by ([Kane1985_], pg. 91):

.. math::
   :label: eq-moment-another-point

   \bar{M}^{S/P} = \bar{M}^{S/Q} + \bar{r}^{P/Q} \times \bar{R}^{S/Q}

:math:`\bar{R}^{S/Q}` is the resultant of the set of bound vectors about point
:math:`Q`.

.. todo:: Show a figure here to represent the left and right sides of the
   equation.

Given a set of bound vectors and a set of bound vectors that consist of a
torque of a couple :math:`\bar{T}` and vector :math:`\bar{v}` bound to an
arbitrary point :math:`P` it is a necessary and sufficient condition that the
second set is a replacement of the first if:

.. math::
   :label: eq-couple-torque-repl

   \bar{T} = \bar{M}^{S/P} \\
   \bar{v} = \bar{R}^S

This means that every set of bound vectors can be replaced by an equivalent
torque of a couple and a single bound vector that is the resultant of the
replaced set.  This is the simplest replacement and simplifies the descritpion
of forces acting on bodies and particle sets.

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

   R = Ffl*FL.y + Ffr*FR.y + Frl*RL.y + Frr*RR.y
   R

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
we will require a force vector and a point to fully describe these bound
vectors. Methods and functions in SymPy mechanics that make use of forces will
typically require a tuple containing a vector and a point, for example the
resultant force acting on the mass center of of the car :math:`B_o` would be
specified like so:

.. jupyter-execute::

   Bo = me.Point('Bo')

   (R, Bo)

Torques of a couple are free vectors (not bound to a line of action) but
represent the couple acting on a rigid body (or set of particles), thus a
reference frame associated with a rigid body will be used to describe the
torque.

.. jupyter-execute::

   (T, B)

We will refer to both forces and torques as *loads*.

Equal & Opposite
----------------

Both forces and torques applied to a multibody system must obey `Newton's Third
Law`_, i.e. forces and toorques act equal and opposite. Take for example a
torque from a motor that causes a pinned lever :math:`B` to rotate relative to the
ground :math:`N`. 

.. jupyter-execute::

   N = me.ReferenceFrame('N')
   B = me.ReferenceFrame('B')

   T = me.dynamicsymbols('T')

   Tm = T*N.z

   (Tm, B), (-Tm, N)

.. warning::

   Careful about your sign convention. It is equally valid to choose `(-Tm, B),
   (Tm, N)`. But it is useful to choose a sign convention such that when the
   signs of angular velocity and torque are the same it means power into the
   system (from the motor in this case). So `B.orient_axis(N, q, N.z)`
   corresponds to `(T*N.z, B)` to power in with both are positive or both are
   negative. This is just a convention though and the choice of force and
   torque signs can be anything, just make sure you know and understand what it
   is!

.. _Newton's Third Law: https://en.wikipedia.org/wiki/Newton's_laws_of_motion#Third_law

Contributing and Non-contributing
---------------------------------

Contributing forces and torques are those that do work on the multibody system.
Non-contributing forces and torques do no work on the system. An example of
non-contributing forces are the contact forces between two rigid bodies if the
two bodies are connected at a single point.

Gravity
-------

We will often be interested in a multibody systems motion when it is subject to
gravitational forces. The simplest case is a constant unidirectional
gravitional field, which is appropriate model for small objects moving about on
and near the Earth's surface. The gravitational forces can be applied to the
mass center of each rigid body or particle in a multibody system. See
[Kane1985]_ pg. XX for the more general case of Newton's Law of Gravitation
which often comes into play for modeling spacecraft.

.. jupyter-execute::

   m, g = sm.symbols('m, g')

   Fg = -m*g*N.y

Springs & Dampers
-----------------

Idealized springs and dampers are useful models of elements that have distance
and velocity depedent forces and torques. A spring with free length :math:`q_0`
and :math:`q_1,q_2` locate the ends of the spring along a line parallel to the
:math:`\hat{n}_x` direction taking a sign convention that a positive spring
force acting on the :math:`q_2` end of the spring is in the negative
:math:`\hat{n}_x` direction . If the spring is linear with stiffness :math:`k`
the spring force vector is then:

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

Friction
--------

Coulomb's Law is the simplest model of friction constant friction.

.. jupyter-execute::

   mu, m, g = sm.symbols('mu, m, g')

   Fn = m*g

   Ff = sm.Piecewise((mu*Fn, displacement.diff(t) > 0),
                     (0, True),
                     (-mu*Fn, displacement.diff(t) < 0))*N.x
   Ff

.. jupyter-execute::

   Ff = mu*Fn*sm.sign(displacement.diff(t))*N.x
   Ff

Aerodynamic
-----------

Aerodynamic drag of a blunt body is dominated by the frontal area drag. If
:math:`P` is always located at the aerodynamic center of a body then the
aerodynamic drag is:

.. jupyter-execute::

   A, Cd, rho = sm.symbols('A, Cd, rho')
   ux, uy, uz = me.dynamicsymbols('u_x, u_y, u_z')

   N_v_P = ux*N.x + uy*N.y + uz*N.z

   Fd = -N_v_P.normalize()*Cd*A*rho/2*N_v_P.dot(N_v_P)
   Fd

Collision
---------

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
