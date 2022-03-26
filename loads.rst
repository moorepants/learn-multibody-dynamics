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
   sm.init_printing(use_latex='mathjax')

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

Specifying Forces and Torques
=============================

Equal & Opposite

Gravity

Springs & Dampers

Friction

Aerodynamic

Contact
