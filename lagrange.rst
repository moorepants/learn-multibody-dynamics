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

When learning about Kane's method, we learned that the negated generalized inertial
forces equal the applied forces, see :ref:`Unconstrained Equations of Motion`.
A large part of Kane's method of deriving the equations of motions for a 
system is involved with finding the generalized inertial forces.

As an alternative, the following equation also get the generalized inertial forces of a
system, now by starting from the `kinetic energy TODO`_ :math:`T(\dot{\bar{q}}, \bar{q})`
expressed as function of the generalized coordinates :math:`\bar{q}`, and 
their time derivatives.

.. _`kinetic energy TODO`: https://en.wikipedia.org/wiki/Work_in_process


.. math::
   :label: eq-lagrange-inertial

    -\bar{F}^*_r = \frac{\mathrm{d}}{\mathrm{d}t}\left(\frac{\partial T}{\partial u_r}
        \right) - \frac{\partial T}{\partial q_r}

 
.. warning:: Note the two minus signs in the above equation

.. note::

    In Kane's method, we are free to choose generalized speeds :math:`\bar{u}` that differ from
    the time derivatives of the generalized coordinates, that is :math:`\dot{\bar{q}}`. This
    freedom is not there when using the Lagrange method. Therefore we will consistently use
    :math:`\dot{\bar{q}}` in this chapter.

The generalized inertial forces computed in this manner are the same as when following
Kane's method. This can be shown by carefully matching terms in both formulations, as
is done for a system of point-masses in [Vallery2020]_.

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
involving :math:`\ddot{q}_r` are mangled. It is common to extract the system
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

* linearized gravity on the surface of the earth, with potential :math:`m g h(\bar{q})`,
* gravity from Newton's universal gravitation, with potential :math:`-G \frac{m_1m_2}{r(\bar{q})}`,
* a linear spring, with potential :math:`\frac{1}{2}k(l(\bar{q}) - l_0)`.

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

When using Kane's method, constraints are handled by dividing the generalized speeds into two sets:
the dependent and independent generalized speeds. Then, the dependent generalized speeds are eliminated 
by solving the (time derivative of the) constraint equation.

In the Lagrange method, the generalized speeds should always mach the generalized coordinates.
Therefore, to handle constraints, the generalized coordinates should be likewise eliminated. This
is not possible for non-holonomic constraint (by definition), and requires to solve often difficult
non-linear equations when considering holonomic constraints. This method of elimination is therefore
not useful within the Lagrange method.

Instead, we generalize the approach in :ref:`Exposing Noncontributing Forces`. We will first ommit the
constraint, and add a constraint force, for which we can specify the direction, but not the magnitude. 
The (second) time derivative of the constraint equation is then added to the equations found with the
Euler-Lagrange equation.

For example, a four bar linkage:






It can be tricky to find the direction of the constraint force from the geometric of the system directly.
There is a trick, called the method of the Lagrange multupliers, to quickly find the correct generalized
forces associated with the constraint forces. 

Given a constraint in the general form

.. math::

    \sum_r a_r(\bar{q}) \dot{q}_r = 0

We find the generalized force as:

.. math::

    F_r = \lambda a_r(\bar{q})

Here :math:`\lambda` is a variable encoding the magnitude of the constraint force. It is
called  the Lagrange multiplier. The same :math:`lambda`` is used for each :math:`r`, that is, 
one constraint has one associated Lagrange multiplier.


**Example: turning the freely floating body discussed earlier into a rolling sphere.**

The non-slip condition of the rolling sphere is split into three constraints: the velocity of
the contact point (:math:`G`) is zero in both the :math:`\hat{n}_x`, :math:`\hat{n}_y` and :math:`\hat{n}_z`
direction. These constraints are enforced by contact forces in their respective directions.

The contact point can be found according by :math:`\bar{r}^{G/C} = -r \hat{n}_z`. We therefore get the
constraint:

.. math::

    \begin{align*}
        \bar{n}_x\cdot ({}^N\bar{v}^C + {}^N\bar{\omega}^B \times -r\hat{n}_z) &= 0 \\
        \bar{n}_y\cdot ({}^N\bar{v}^C + {}^N\bar{\omega}^B \times -r\hat{n}_z) &= 0 \\
        \bar{n}_z\cdot ({}^N\bar{v}^C + {}^N\bar{\omega}^B \times -r\hat{n}_z) &= 0 \\
    \end{align*}

These can be used in the Lagrange-multiplier method as follows:

.. jupyter-execute:

    1 + 1



    
The method of the Lagrange multiplier can of course also be used within Kane's method. However,
this results in a larger system of equations, which is why the elimination approach is often
preferred there. An exception being scenarios where the constraint force itself is a useful output,
for instance to check no-slip conditions in case of limited friction.


Lagrange's vs Kane's
====================

Why did we learn both methods?

* Kane i


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


The generalized momenta are an invertable function of the generalized speeds. We can therefore replace the
Lagrangian equation by:

.. math::

    \dot{p_r} = \frac{\partial L}{\partial q_r}

.. math::

    \dot{q_r} = \dot{q_r}(\bar{p})  

which are equivalent to the equations obtained using Hamilton's method. Hamiltonian systems and their
extension Port-Hamiltonian system are often used in physics and control theory respectively.



(Learn more) Euler-Lagrange in optimization
===========================================

There is an 

Euler-Lagrange as key result in variational calculus.

We can optimize more things + references to other advanced concepts.






