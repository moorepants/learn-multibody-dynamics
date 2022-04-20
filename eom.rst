===================
Equations of Motion
===================

In the previous chapter, we introduced the generalized active forces and the
generalized inertia forces. Together, these two pieces give us the dynamical
differential equations. The dynamical differential equations are:

.. math::
   :label: eq-kanes-equations

   \bar{F}_r + \bar{F}^*_r = \bar{f}_d(\dot{\bar{u}}, \bar{u}, \bar{q}, t)  = 0

We also call these equations *Kane's Equations* due to the formulation
presented in [Kane1985]_.

:math:`\bar{F}^*_r` is linear in the time derivatives of the generalized speeds
and also contains velocitiy dependent terms such as centripal and Coriolis
forces and rotational velocity couplings. Dynamics texts will often present
it in this form:

.. math::

   \mathbf{M}(q, t) \dot{\bar{u}} + \bar{C}(\bar{u}, \bar{q}, t)

where :math:`\mathbf{M}` is called the *mass matrix* and :math:`\bar{C}` is are
the forces due to velocity effects.

The kinematical and dynamical differential equations constitute the *equations
of motion* for a holonomic multibody system. These equations are ordinary
differntial equations in the generalized speeds and generalized coordinates.
They are also linear in :math:`\dot{\bar{u}}` and :math:`\dot{\bar{q}}`

.. math::
   :label: eq-equations-of-motion

   \bar{f}_d(\dot{\bar{u}}, \bar{u}, \bar{q}, t)  = 0 \\
   \bar{f}_k(\dot{\bar{q}}, \bar{u}, \bar{q}, t)  = 0


