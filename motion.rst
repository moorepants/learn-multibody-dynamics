========================
Nonholonomic Constraints
========================

In the prior chapter we discussed constraints on configuration of a system. The
configuration considers where where points are and how reference frames are
oriented. In this chapter, we will consider constraints on the motion of a
system.

Take for example parallel parking a car as a motivating example.

.. _motion-parallel:
.. figure:: figures/motion-parallel.svg
   :align: center

We know that car 2 can be in either the left or right location in a), i.e. the
car's configuration permits either location. But the scenario in b) isn't
possible. A car can't move from the left configuration to the right
configuration by simply moving to the right [*]_. Although, this surely would
be nice if we could. A car has wheels and only the front wheels can be steered,
so the scenario in c) is the only way for the car to end up in the right
configuration. The car has to move in a specific way to get from one
configruation to another. This entails that we have some kind of constraint on
the motion but not the configuration. Constraints such as these are called
*nonholonomic constraints* and they take the form:

.. math::
   :label: nonholonomic-constraints

   \bar{f}_n(\bar{u}, \bar{q}, t) = 0 \\
   \textrm{ where } \\
   \bar{f}_n \in \mathbb{R}^m \\
   \bar{u} = \left[ u_1, \ldots, u_n\right]^T \in \mathbb{R}^n\\
   \bar{q} = \left[ q_1, \ldots, q_n\right]^T \in \mathbb{R}^n

Kinematical Differential Equations
==================================

The variables :math:`u_1, \ldots, u_n` are defined as linear functions of the
time derivatives of the generalized coordinates :math:`\dot{q}_1, \ldots,
\dot{q}_n`. These variables are called generalized speeds. They take the form:

.. math::
   :label: generalized-speeds

   \bar{u} := \mathbf{Y}_k \dot{\bar{q}} + \bar{z}_k(\bar{q}, t)

:math:`\bar{u}` must be chosen such that :math:`\mathbf{Y}_k` is invertible.
Eq. :math:numref:`generalized-speeds` are called *kinematical differential
equations*. The most common, and always valid, choice of generalized speeds
is:

.. math::
   :label: generalized-speeds

   \bar{u} := \mathbf{I} \dot{\bar{q}}

where :math:`u_i = \dot{q}_i` for :math:`i=1,\ldots,n`.

.. [*] Well, we could find a very strong person to push th ecar sideways,
   overcoming the very high resisting friction force.
