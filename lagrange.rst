============================================
Equations of Motion with the Lagrange Method
============================================

.. note::

   You can download this example as a Python script:
   :jupyter-download-script:`tmt` or Jupyter Notebook:
   :jupyter-download-notebook:`lagrange`.

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   me.init_vprinting(use_latex='mathjax')

This book has already discussed three methods to derive the equations
of motion of multibody systems: Newton-Euler, Kane's method and the TMT
method. This chapter will add a third: the Lagrange method, originally 
developed by Joseph-Louis Lagrange.



Inertial forces from kinetic energy
===================================

When learning about Kane's method, we learned that the negated inertial
forces equal the applied forces, see :ref:`Unconstrained Equations of Motion`.
A large part of Kane's method of deriving the equations of motions for a 
system is involved with finding the inertial forces.

As an alternative, the following equation also get the inertial forces of a
system, now by starting from the `kinetic energy TODO`_ :math:`T(\bar{u}, \bar{q})`
expressed as function of the generalized coordinates :math:`\bar{q}`, and 
their corresponding generalized speeds :math:`\bar{u}`

.. _`kinetic energy TODO`: https://en.wikipedia.org/wiki/Work_in_process


.. math::
   :label: eq-lagrange-inertial

    -\bar{F}^*_r = \frac{\mathrm{d}}{\mathrm{d}t}\left(\frac{\partial T}{\partial u_r}
        \right) - \frac{\partial T}{\partial q_r}

The derivation of this formula is discussed in the section
:ref:`Lagrange equation from virtual work`.
 
.. warning:: Note the two minus signs in the above equation

Example: 3D body in outer space


.. jupyter-execute::
    # Setting up reference frames
    psi,theta, phi, x, y, z = me.dynamicsymbols('psi theta phi x y z')
    N = me.ReferenceFrame('N')
    B = me.ReferenceFrame('B')
    B.orient_body_fixed(N, (psi, theta, phi), 'zxy')

    # Mass and inertia
    m, Ixx, Iyy, Izz = sm.symbols('M, I_{xx}, I_{yy}, I_{zz}')
    I_B = me.inertia(B, Ixx, Iyy, Izz)

    # Kinematics and kinetic energy

    omega_B = B.ang_vel_in(N)
    r_com = x*N.x + y*N.y + z*N.z
    v_com = r_com.dt(N)
    T = omega_B.dot(I_B.dot(omega_B))/2 + m*v_com.dot(v_com)/2

    # Euler-Lagrange equation

    t = me.dynamicsymbols._t
    q = sm.Matrix([psi, theta, phi, x, y, z])
    qdot = q.diff(t)
    qddot = qdot.diff(t)
    p = sm.Matrix([T]).jacobian(qdot).transpose()
    g = -sm.Matrix([T]).jacobian(q).transpose()
    left_hand_side = p.diff(t) + g


This gives the equations of motion, but the terms, particularly the terms
involving :math:`\dot{u}_r` are mangled. It is common to extract the system
mass matrix and velocity forces vector like so:

.. jupyter-execute::

    mass_matrix = left_hand_side.jacobian(qddot)
    dynamic_bias = left_hand_side - mass_matrix*qddot


Conservative Forces
===================

Some applied forces, known as conservative forces `conservative forces`_, can
be expressed using the gradient of a scalar function of the generalized coordinates,
known as the `potential energy`_ :math:`V(\bar{q})`:

.. math::
    :label: eq-potential-energy

    \bar{F}_r = -\frac{\partial V}{\partial q_r}

.. warning:: Note the minus sign in the above equation.

.. _`conservative forces`: https://en.wikipedia.org/wiki/Conservative_force
.. _`potential energy`: https://en.wikipedia.org/wiki/Potential_energy

Some examples of conservative forces are:

* linearized gravity on the surface of the earth, with potential :math:`m g h`,
* gravity from Newton's universal gravitation, with potential :math:`-G \frac{m_1m_2}{r}`,
* a linear spring, with potential :math:`\frac{1}{2}k(l - l_0)`.

For conservative forces, it is often convenient to derive the applied forces via 
the potential energy.


The Lagrange-method
===================

Both the equation for computing the inertial forces from the kinetic energy, and 
the equation for computing the applied forces from a potential energy have a term
in them with the partial derivative with respect to the generalized coordinate. 
Furtermore, the potential energy does not depend on the generalized speeds. 
Therefore, we can derive the resulting (inertial and conservative applied) forces
in one go, by combining the two equations.

Step 1. Compute the so called Lagrangian :math:`L`, the difference between the 
kinetic energy and potential energy:

.. math::
    :label: eq-lagrangian

    L = T - V

Step 2. Use the Euler-Lagrange equations (the name for the equation 
:ref:`eq-lagrange-inertial`) to find the equations of motion:

.. math::
    :label: eq-euler-lagrange

    \frac{\mathrm{d}}{\mathrm{d}t}\left(\frac{\partial L}{\partial u_r}
        \right) - \frac{\partial L}{\partial q_r} = \bar{F}_r,
    
while being careful to include a force either in the applied forces 
:math:`\bar{F}_r`, or in the potential energy :math:`V`, but never
in both.


EXAMPLE: the same example used throughout, with funny double pendulum etc

Note that when we extracted the mass matrix from the left hand side of these
equations, the residual is not just the velocity force vector, but also
includes the conservative forces.



Constrained equations of motion
===============================

Compute the applied forces as before, add the constraint equation as before


Lagrange equation from virtual work
===================================

Show that the Euler-Lagrange equation leads to the same results as virtual work,
for a system consisting of :math:`n` particles.

(Learn more) Generalized momentum
=================================

The partial derivative of the Lagrangian with respect to generalized speed is
called the generalized momentum.

Examples showing that this matches to momentum and angular momentum in relevant 
particle cases.

If the lagrangian does not depend on a generalized coordinates, its associated
generalized momentum is conserved.

Example to show "conservation of angular momentum", and rotating body like falling cat:
 * body 1 can rotate wrt ground around same axis as gravity (say z)
 * (massless) body 2 can rotate wrt body 1 around same axis as gravity
 * body 3 can rotate wrt body 2 around a (body fixed) axis perpendicular to gravity (say x)
This example will also show how to apply motor torques at joints.

Practice problem: add a damping force or a coulomb friction force in the first joint 
(the example and this problem are inspired by a talk by A. Ruina, https://www.youtube.com/watch?v=j-wHI764dWU)




(Learn more) Euler-Lagrange in optimization
===========================================

Euler-Lagrange as key result in variational calculus.

We can optimize more things + references to other advanced concepts.






