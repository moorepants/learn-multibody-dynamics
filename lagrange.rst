============================================
Equations of Motion with the Lagrange Method
============================================

.. note::

   You can download this example as a Python script:
   :jupyter-download-script:`lagrange` or Jupyter Notebook:
   :jupyter-download-notebook:`lagrange`.

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

- Derive the Lagrangian for a system of particles and rigid bodies
- Use the Euler-Lagrange equation to derive equations of motions given a Lagrangian
- Use the method of Lagrange multipliers to add constraints to the equations of motions
- Know how to determine the generalized momenta of a system.

Introduction
============

This book has already discussed two methods to derive the equations
of motion of multibody systems: Newton-Euler and Kane's method. This
chapter will add a third: the `Lagrange method`_, originally 
developed by Joseph-Louis Lagrange. These materials focus on Engineering
applications for multi-body systems, and therefore build the Lagrange method around
the terms found earlier in Kane's equations. In other textbooks, the Lagrange method
is often derived from the `Variational principles`_, such as virtual work or the principle
of least action. A good starting point for studying the physical
and mathematical background of the Lagrange approach is [Lanczos1970]_. 

.. _Variational principles: https://en.wikipedia.org/wiki/Variational_principle
.. _`Lagrange method`: https://en.wikipedia.org/wiki/Lagrangian_mechanics



Inertial forces from kinetic energy
===================================

In Kane's method the negated generalized inertial
forces equal the applied forces, see :ref:`Unconstrained Equations of Motion`.
A large part of Kane's method of deriving the equations of motions for a 
system is involved with finding the generalized inertial forces.

As an alternative, the following equation also calculates the generalized inertial forces of a
system, now by starting from the kinetic energy :math:`K (\dot{\bar{q}}, \bar{q})`
expressed as function of the generalized coordinates :math:`\bar{q}`, and 
their time derivatives. See :ref:`Energy and Power` for the definition of kinetic energy.

.. math::
  :label: eq-lagrange-inertial

   -\bar{F}^*_r = \frac{\mathrm{d}}{\mathrm{d}t}\left(\frac{\partial K}{\partial \dot{q}_r}
        \right) - \frac{\partial K}{\partial q_r}

.. warning:: Note the two minus signs in the above equation

.. note::

   In Kane's method, it is possible to choose generalized speeds :math:`\bar{u}` that differ from
   the time derivatives of the generalized coordinates :math:`\dot{\bar{q}}`. By convention
   :math:`\bar{u} = \dot{\bar{q}}` is assumed when using the Lagrange method. Therefore, throughout
   this chapter :math:`\dot{\bar{q}}` is used.

The generalized inertial forces computed in this manner are the same as when following
Kane's method, or the TMT method used in the next chapter. This can be shown by carefully
matching terms in these formulations, as is done for a a system of point-masses in [Vallery2020]_.

Example: freely moving 3D body
------------------------------

This example is largely the same as the example in :ref:`Body Fixed Newton-Euler Equations`. A key difference
is a difference between the generalized speeds describing the rotation. In the calculation with Kane's method,
they were body-fixed angular velocities, whereas here they are the rates of change of the Euler angles. 

First, set up the generalized coordinates, reference frames and mass properties:

.. jupyter-execute::

   t = me.dynamicsymbols._t
   psi,theta, phi, x, y, z = me.dynamicsymbols('psi theta phi x y z')
   q = sm.Matrix([psi, theta, phi, x, y, z])
   qd = q.diff(t)
   qdd = qd.diff(t)
   N = me.ReferenceFrame('N')
   B = me.ReferenceFrame('B')
   B.orient_body_fixed(N, (psi, theta, phi), 'zxy')
   m, Ixx, Iyy, Izz = sm.symbols('M, I_{xx}, I_{yy}, I_{zz}')
   I_B = me.inertia(B, Ixx, Iyy, Izz)
   q

Then compute the kinetic energy: 

.. jupyter-execute::

   N_w_B = B.ang_vel_in(N)
   r_O_P = x*N.x + y*N.y + z*N.z
   N_v_C = r_O_P.dt(N)
   K = N_w_B.dot(I_B.dot(N_w_B))/2 + m*N_v_C.dot(N_v_C)/2
   K

Use the kinetic energy to find the generalized inertial forces. Here we start with 
the generalized coordinate :math:`\psi`

.. jupyter-execute::

    psid = psi.diff(t)
    F_psi_s = K.diff(psid).diff(t) - K.diff(psi)

A similar derivation should be made for all generalized coordinates. We could write
a loop, but there there is a method to derive all the equations in one go.
The vector of partial derivatives of a function, that is the gradient, can be created
using the :external:py:meth:`~sympy.matrices.matrices.MatrixCalculus.jacobian` method. The generalized inertial forces can then be found like this: 

.. jupyter-execute::

   K_as_matrix = sm.Matrix([K])
   Fs_transposed = K_as_matrix.jacobian(qd).diff(t) - K_as_matrix.jacobian(q)
   Fs = Fs_transposed.transpose()
   Fs


While these are correct equations of motion, the terms, particularly the terms
involving :math:`\ddot{q}_r` are mangled. It is common to extract the system
mass matrix :math:`\mathbf{M}_d` and velocity forces vector :math:`\bar{g}_d` like before:

.. jupyter-execute::

   qdd_zerod = {qddr: 0 for qddr in qdd}
   Md = Fs.jacobian(qdd)
   gd = Fs.xreplace(qdd_zerod)
   Md.simplify()
   gd.simplify()
   Md, gd


Conservative Forces
===================

Recall from :ref:`Energy and Power` that `conservative forces`_, can
be expressed using the gradient of a scalar function of the generalized coordinates,
known as the `potential energy`_ :math:`V(\bar{q})`:

.. math::
   :label: eq-potential-energy

   \bar{F}_r = -\frac{\partial V}{\partial q_r}

.. warning:: Note the minus sign in the above equation.

.. _`conservative forces`: https://en.wikipedia.org/wiki/Conservative_force
.. _`potential energy`: https://en.wikipedia.org/wiki/Potential_energy

Some examples of conservative forces are:

* a uniform gravitational field, for example on the surface of the earth, with potential :math:`V = m g h(\bar{q})`,
* gravity from Newton's universal gravitation, with potential :math:`V = -G \frac{m_1m_2}{r(\bar{q})}`,
* a linear spring, with potential :math:`V = \frac{1}{2}k(l(\bar{q}) - l_0)^2`.

For conservative forces, it is often convenient to derive the applied forces via 
the potential energy.


The Lagrange Method
===================

Both the equation for computing the inertial forces from the kinetic energy, and 
the equation for computing the applied forces from a potential energy have a term
in them with the partial derivative with respect to the generalized coordinate. 
Furtermore, the potential energy does not depend on the generalized speeds. 
Therefore, the resulting (inertial and conservative applied) forces can be derived
in one go, by combining the two equations.

Step 1. Compute the so called `Lagrangian`_ :math:`L`, the difference between the 
kinetic energy and potential energy:

.. math::
   :label: eq-lagrangian

   L = K - V

.. _`Lagrangian`: https://en.wikipedia.org/wiki/Lagrangian

Step 2. Use the Euler-Lagrange equations (the name for the equation 
:math:numref:`eq-lagrange-inertial`) to find the equations of motion:

.. math::
   :label: eq-euler-lagrange

   \frac{\mathrm{d}}{\mathrm{d}t}\left(\frac{\partial L}{\partial \dot{q}_r}
     \right) - \frac{\partial L}{\partial q_r} = F_r,

while being careful to include a force either in the applied forces 
:math:`\bar{F}_r`, or in the potential energy :math:`V`, but never
in both.


Example: Double pendulum with springs and sliding pointmass
-----------------------------------------------------------

This example will use the Lagrange method to derive the equations of motion 
for the system introduced in :ref:`Example of Kane's Equations`. The description
of the system is shown again in :numref:`fig-eom-double-rod-pendulum-repeat`.

.. _fig-eom-double-rod-pendulum-repeat:
.. figure:: figures/eom-double-rod-pendulum.svg
   :align: center
   :width: 600px

   Three dimensional pendulum made up of two pinned rods and a sliding mass on
   rod :math:`B`. Each degree of freedom is resisted by a linear spring. When
   the generalized coordinates are all zero, the two rods are perpendicular to
   each other.

The first step is to define the relevant variables, constants and frames. This step
is the same as for Kane's method.

.. admonition:: Frames and Bodies Setup
   :class: dropdown

    .. jupyter-execute::

        m, g, kt, kl, l = sm.symbols('m, g, k_t, k_l, l')
        q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')

        N = me.ReferenceFrame('N')
        A = me.ReferenceFrame('A')
        B = me.ReferenceFrame('B')

        A.orient_axis(N, q1, N.z)
        B.orient_axis(A, q2, A.x)

        O = me.Point('O')
        Ao = me.Point('A_O')
        Bo = me.Point('B_O')
        Q = me.Point('Q')

        Ao.set_pos(O, l/2*A.x)
        Bo.set_pos(O, l*A.x)
        Q.set_pos(Bo, q3*B.y)

        O.set_vel(N, 0)

        I = m*l**2/12
        I_A_Ao = I*me.outer(A.y, A.y) + I*me.outer(A.z, A.z)
        I_B_Bo = I*me.outer(B.x, B.x) + I*me.outer(B.z, B.z)

Then, set up the Lagrangian:

.. jupyter-execute::

   t = sm.symbols('t')
   q = sm.Matrix([q1, q2, q3])
   qd = q.diff(t)
   qdd = qd.diff(t)

   K = m/2*(Ao.vel(N).dot(Ao.vel(N)) + Bo.vel(N).dot(Bo.vel(N)) + Q.vel(N).dot(Q.vel(N))) + (
       A.ang_vel_in(N).dot(I_A_Ao.dot(A.ang_vel_in(N))) + B.ang_vel_in(N).dot(I_B_Bo.dot(B.ang_vel_in(N)))
   )/2
   V = m*g*(Ao.pos_from(O).dot(-N.x) + Bo.pos_from(O).dot(-N.x)) + kt/2*(q1**2) + kt/2*q2**2 + kl/2*q3**2

   L = sm.Matrix([K - V])
   sm.trigsimp(L)

Finally, derive the equations of motion:

.. jupyter-execute::

   left_hand_side = L.jacobian(qd).diff(t) - L.jacobian(q)
   qdd_zerod = {qddr: 0 for qddr in qdd}
   Md = left_hand_side.jacobian(qdd)
   gd = left_hand_side.xreplace(qdd_zerod)
   me.find_dynamicsymbols(Md), me.find_dynamicsymbols(gd)

The mass matrix :math:`\mathbf{M}_d` only depends on :math:`\bar{q}`, and :math:`\bar{g}_d` depends
on :math:`\dot{\bar{q}}` and :math:`\bar{q}`, just as in Kane's method. Note that :math:`\bar{g}_d` now
combines the effects of the velocity force vector and the conservative forces. In this setting, 
:math:`\bar{g}_d` is often called the dynamic bias. 

It is often useful to use a vector of intermediate variables when finding the Euler-Lagrange equations. The variables
are defined as:

.. math::

    p_r = \frac{\partial L}{\partial \dot{q_r}}

The variables are collected in a vector :math:`\bar{p}`. 

They are called the generalized momenta, 
as they coincide with linear momentum in the case of a Lagrangian describing a particle. 
Similar to the situation in the dynamics of particles, there can 
be conservation of generalized momentum. This is the case for the generalized momentum 
associated with ignorable coordinates, as defined in :ref:`Equations of Motion with Nonholonomic Constraints`. 

For the example pendulum, the generalized momenta are calculated as:

.. jupyter-execute::

    p = L.jacobian(qd).transpose()
    sm.trigsimp(p)


Constrained equations of motion
===============================

When using Kane's method, constraints are handled by dividing the generalized speeds into two sets:
the dependent and independent generalized speeds. Depending on the type of constraints, the 
dependent generalized speeds are eliminated by solving the constraint equation (for non-holonomic 
constraints) or the time derivative of the constraint equation (holonomic constraints). Kane's 
method only gives rise to :math:`p = n - m` dynamical equations, one for each independent generalized
speed.

The Lagrange method gives rise to :math:`n` dynamical equations, one for each generalized coordinate. 
To eliminate the constraints, and end up with the right number of equations (:math:`n - m`, one for
each degree of freedom), both the generalized speeds and the generalized coordinates should be solved 
using the constraint equation. For non-holonomic constraints, this elimination is not possible (by the
definition of non-holonomic), and for holonomic constraints this elimination requires solving often 
difficult non-linear equations for the generalized coordinates. The method of elimination is therefore 
not useful within the Lagrange method.

Instead, constraints are handled using a generalized version of the approach in 
:ref:`Exposing Noncontributing Forces`. First the constraints are omitted. Then a constraint force is added,
with a known direction, but unknown magnitude. Finally, the (second) time derivative of the constraint 
equation is then appended to the equations found with the Euler-Lagrange equation.

For example, consider a particle of mass :math:`m` and position 
:math:`\bar{r}^{P/O} = q_1 \hat{n}_x + q_2 \hat{n}_y + q_3\hat{n}_z`
on a slope :math:`q_1 = q_2`.  The unconstrained Lagrangian is 
:math:`L = \frac{1}{2}m(\dot{q}_1^2 + \dot{q}_2^2 + \dot{q}_3^2) - mgq_3`.
The constraint force is perpendicular to the slope, so is described
as :math:`\bar{F} = F\hat{n}_x - F\hat{n}_y`. The appended equation is
the second time derivative of the constraint equation :math:`\ddot{q_1} - \ddot{q_2} = 0`.
Combining all, gives:

.. math::
    \begin{array}{r}
    m\ddot{q}_1= \phantom{-}F\\
    m\ddot{q}_2= -F\\
    m\ddot{q}_3 + mg = \phantom{-}0\\
    \ddot{q}_1 - \ddot{q}_2\!\! = \phantom{-}0 
    \end{array}

This can be put in matrix-form, by extracting the unknown acceleration and force magnitude:

.. math::
    \begin{bmatrix} m & 0 & 0 &-1 \\ 0 & m & 0 & 1 \\ 0 & 0 & m & 0 \\ 1 & -1 & 0 & 0\end{bmatrix}
    \begin{bmatrix} \ddot{q}_1 \\ \ddot{q}_2 \\ \ddot{q}_3 \\ F \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \\ -mg \\ 0\end{bmatrix}


It can be challenging to find the direction of the constraint force from the geometry of the system directly.
There is a trick, called the method of the `Lagrange multipliers`_, to quickly find the correct generalized
forces associated with the constraint forces. 

.. _`Lagrange multipliers`: https://en.wikipedia.org/wiki/Lagrange_multiplier 

Given a motion constraint (time derivatives of configuration constraint or a nonholonomic constraint) in the general form

.. math::

    \sum_r a_r(\bar{q}) \dot{q}_r = 0

The generalized force is found as:

.. math::

    F_r = \lambda a_r(\bar{q})

Here :math:`\lambda` is a variable encoding the magnitude of the constraint force. It is
called  the Lagrange multiplier. The same :math:`\lambda` is used for each :math:`r`, that is, 
each constraint has a single associated Lagrange multiplier.

Due to how it is constructed, the power produced by the constraint force is always zero, as expected.

.. math::

    P = \sum_r F_r\dot{q}_r = \sum \lambda a_r(\bar{q})\dot{q}_r  = \lambda \sum a_r(\bar{q})\dot{q}_r = \lambda \cdot 0

For example, consider the pointmass to be constrained to move in a bowl
:math:`q_1^2 + q_2^2 + q_3^2 -1 = 0`, :numref:`fig-lagrange-bowl`.  Taking the
time derivative gives: :math:`a_1 = 2q_1`, :math:`a_2 = 2q_2`, and :math:`a_3 =
2q_3`.  This results in generalized reaction forces :math:`F_1 = 2\lambda q_1`,
:math:`F_2 = 2\lambda q_2` and :math:`F_3 = 2\lambda q_3`.

.. _fig-lagrange-bowl:
.. figure:: figures/lagrange-bowl.svg
   :align: center

   Point mass :math:`P` constrained to the surface of a spherical bowl with
   radius :math:`1` and constraint force measure numbers :math:`F_1,F_2,F_3`.

Often, there are multiple constraints on the same system. For convenience, the handling of these constraints can be combined.
Consider the :math:`m+M` dimensional general constraint equations consisting of the time derivatives of the holonomic constraints
and/or the non-holonomic constraints:

.. math::

    \bar{f}_{hn}(\bar{q}, \bar{\dot{q}}) = \mathbf{M}_{hn}\bar{\dot{q}} = 0,

the combined constraint forces are given as:

.. math::

    \bar{F}_r = \mathbf{M}_{hn}^\text{T}\bar{\lambda},

where :math:`\bar{\lambda}` is a vector of :math:`m + M` Lagrange multipliers, one for each constraint (row in :math:`\mathbf{M}_{hn}`).


Example: turning the freely floating body discussed earlier into a rolling sphere.
----------------------------------------------------------------------------------

The non-slip condition of the rolling sphere is split into three constraints: the velocity of
the contact point (:math:`G`) is zero in the :math:`\hat{n}_x`, the :math:`\hat{n}_y` and  the :math:`\hat{n}_z`
direction. The first two constraints are non-holonomic, the last constraint is the time derivative of
a holonomic constraint. All three constraints are enforced by contact forces in their respective directions.

The contact point can be found according by :math:`\bar{r}^{G/C} = -r \hat{n}_z`. By using the :ref:`Velocity
Two Point Theorem`, the following constraints are found.

.. math::

    \begin{array}{l}
    \bar{n}_x\cdot ({}^N\bar{v}^C + {}^N\bar{\omega}^B \times (-r\hat{n}_z)) = 0 \\
    \bar{n}_y\cdot ({}^N\bar{v}^C + {}^N\bar{\omega}^B \times (-r\hat{n}_z)) = 0 \\
    \bar{n}_z\cdot ({}^N\bar{v}^C + {}^N\bar{\omega}^B \times (-r\hat{n}_z)) = 0 \\
    \end{array}

These can be used to derive the constraint force and the additional equations using the Lagrange-multiplier
method as shown below. Note that here only the first time derivative of the constraint equation is used, 
again because the second time derivatives of the generalized coordinates appear.

.. admonition:: Frames and Body Setup
   :class: dropdown

    Setting up reference frames
   
    .. jupyter-execute::

        psi,theta, phi, x, y, z = me.dynamicsymbols('psi theta phi x y z')
        N = me.ReferenceFrame('N')
        B = me.ReferenceFrame('B')
        B.orient_body_fixed(N, (psi, theta, phi), 'zxy')

        # Mass and inertia
        m, Ixx, Iyy, Izz = sm.symbols('M, I_{xx}, I_{yy}, I_{zz}')
        I_B = me.inertia(B, Ixx, Iyy, Izz)

    Finding the kinetic energy:

    .. jupyter-execute::

        omega_B = B.ang_vel_in(N)
        r_com = x*N.x + y*N.y + z*N.z
        v_com = r_com.dt(N)
        K = omega_B.dot(I_B.dot(omega_B))/2 + m*v_com.dot(v_com)/2

    Deriving equations of motion:

    .. jupyter-execute::

        t = me.dynamicsymbols._t
        q = sm.Matrix([psi, theta, phi, x, y, z])
        qd = q.diff(t)
        qdd = qd.diff(t)

        L = sm.Matrix([K])
        left_hand_side = L.jacobian(qd).diff(t) - L.jacobian(q)

        qdd_zerod = {qddr: 0 for qddr in qdd}
        Md = left_hand_side.jacobian(qdd)
        gd = left_hand_side.xreplace(qdd_zerod)

To make this free floating body a rolling wheel, three constraints are needed: the
velocity of the contact point should be zero in :math:`\hat{n}_x`, :math:`\hat{n}_y`
and :math:`\hat{n}_x` direction.

.. jupyter-execute::

    lambda1, lambda2, lambda3 = me.dynamicsymbols('lambda1, lambda2, lambda3') 
    constraint = (v_com + B.ang_vel_in(N).cross(-N.z)).to_matrix(N)
    M_hn = constraint.jacobian(qd)
    diff_constraint = constraint.diff(t)
    sm.trigsimp(constraint)

This constraint information must then be added to the original equations. To do
so, we make use of a useful fact. 

.. jupyter-execute::

    diff_constraint.jacobian(qdd) - M_hn

This equality is true for all constraints, as can easily be shown by taking the time
derivative of the constraint equation, using the chain rule.

The combined equations can now be written in a block matrix form:

.. math::
        \begin{bmatrix} \mathbf{M}_d & -\mathbf{M}_{hn}^T \\ \mathbf{
        M}_{hn} & 0\end{bmatrix}\begin{bmatrix}\ddot{\bar{q}} \\ \bar{\lambda} \end{bmatrix} = 
        \begin{bmatrix} \bar{F}_r - \bar{g}_d \\ - \frac{\partial \mathbf{M}_{hn}\dot{\bar{q}}}{\partial \bar{q}}\dot{\bar{q}} \end{bmatrix},

where :math:`\bar{g}` is the dynamic bias, and the last term on the right hand side, 
called the constraint bias, can be quickly computed as:

.. jupyter-execute::

    constraint_bias = diff_constraint.xreplace({qddr : 0 for qddr in qdd})

We call the block matrix called the extended mass matrix, and the vector on the right hand side the extended dynamic bias. 

With these `n + m + M` equations, it is possible to solve for :math:`\ddot{\bar{q}}` and :math:`\lambda`. It is therefore possible to
integrate/simulate the system directly. However, because only the second derivative of the constraint is satisfied, numerical
errors can build up, so the constraint is not satisfied. It is better to use a differential algebraic solver, as discussed 
in :ref:`Equations of Motion with Holonomic Constraints`. See `the scikit.ode documentation`_ for a worked example.

.. _`the scikit.ode documentation`: https://github.com/bmcage/odes/blob/master/ipython_examples/Planar%20Pendulum%20as%20DAE.ipynb



The method of the Lagrange multiplier can of course also be used within Kane's method. However,
it increases the number of equations, which is why the elimination approach is often
preferred there. An exception being scenarios where the constraint force itself is a useful output,
for instance to check no-slip conditions in case of limited friction.


Lagrange's vs Kane's
====================

The is book has now presented two alternatives to the Newton-Euler method: Kane's method and Lagrange's method. 
This raises the questions: when should each alternative method be used?

For constrained systems, Kane's method has the advantage that the equations of motion are given for a set of
independent generalized speeds only. In other words, Kane's method gives a minimal set of equations. This can
give rise to simplified equations, additional insight, and
numerically more efficient simulation. This also gives the benefit that Lagrange multipliers are not needed
when solving constrained systems with Kane's method.

Furthermore, the connection from Kane's method to vector mechanics, that is, Newton's laws, is clearer, which
can provide additional insight, and make it easier to encorporate non-conservative forces such as friction.

On the other hand, the Lagrange method only requires energies as input, for which only the velocities
of the bodies are needed. Therefore, it can be simpler to derive than the accelerations which are needed for Kane's
method.

Furthermore, the Lagrange method results in a set of equations with well understood structures and properties.
These structures and properties are not studied further in these materials. A starting point for further study
is `Noether's theorem`_, which extends the idea of ignorable coordinates to find conserved quantities like
momentum and energy. 

.. _`Noether's theorem`: https://en.wikipedia.org/wiki/Noether%27s_theorem_





.. (Learn more) Generalized momentum
.. =================================

.. The partial derivative of the Lagrangian with respect to generalized speed is
.. called the generalized momentum.

.. .. math::

..     p = \frac{\partial L}{\partial \dot{\bar{q}}}

.. Some ideas behind generalized momentum will be discussed with the following example,
.. which is a simplified version of the falling cat example:
.. * body A is a cylinder that can rotate wrt ground around same axis as gravity: :math:`\hat{n}_z``
.. * body B is a cylinder that can rotate wrt body A around same axis as gravity
.. * body C is a cylinder that can rotate wrt body C around a (body fixed) axis perpendicular to gravity :math:`\hat{b}_x`
.. * There are two actuators providing a torque on the joints between bodies A and B and bodies B and C respectively.

.. This example will also show how to apply motor torques at joints.

.. .. jupyter-execute::

..    t, l, r, T_b, T_c = sm.symbols('t, l, r, T_b, T_c')
..    q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')

..    N = me.ReferenceFrame('N')
..    A = me.ReferenceFrame('A')
..    B = me.ReferenceFrame('B')
..    C = me.ReferenceFrame('C')

..    A.orient_axis(N, q1, N.z)
..    B.orient_axis(A, q2, A.z)
..    C.orient_axis(B, q3, B.x) 

..    g = 1
..    rho = 1
..    m = rho*l*sm.pi*r**2
..    I_xx_or_yy = m/12*(3*r**2 + l**2)
..    I_zz= m/2*r**2
..    I_A_Ao = me.inertia(A, I_xx_or_yy , I_xx_or_yy, I_zz)
..    I_B_Bo = me.inertia(B, I_xx_or_yy , I_xx_or_yy, I_zz)
..    I_C_Co = me.inertia(C, I_xx_or_yy , I_xx_or_yy, I_zz)

..    O = me.Point('O')
..    O.set_vel(N, 0.0)
..    Ao = me.Point("A_c")
..    Ao.set_pos(O, -0.5*l*A.z)
..    Bo = me.Point("B_c")
..    Bo.set_pos(Ao, -0.5*l*A.z - 0.5*l*B.z)
..    Co = me.Point("C_c")
..    Co.set_pos(Bo, -0.5*l*B.z -0.5*l*C.z)

.. The next step is again to form the Lagrangian and find the equations of motion. As the system has no further constraints, 
.. the Lagrange multiplier method is not needed. The actuator torques are added to the right hand side of the equation, in
.. the same way as active forces are added to Kane's equations. Here the torques are represented by the variables :math:`T_b`
.. and :math:`T_c` are used to represent.

.. .. jupyter-execute::

..     T = m/2*(squarednorm(Ao.vel(N)) + squarednorm(Bo.vel(N)) + squarednorm(Co.vel(N))) + 1/2*(
..             quadraticform(I_A_Ao, A.ang_vel_in(N)) + quadraticform(I_B_Bo, B.ang_vel_in(N)) + quadraticform(I_C_Co, C.ang_vel_in(N)))
..     V = m*g*N.z.dot(Co.pos_from(O))
..     L = sm.Matrix([T - V])

..     q = sm.Matrix([q1, q2, q3])
..     qd = q.diff(t)
..     qdd = qd.diff(t)

..     p = L.jacobian(qd)
..     p.simplify()
..     left_hand_side = (p.diff(t) - L.jacobian(q)).transpose()

..     qdd_zerod = {qddr: 0 for qddr in qdd}
..     Md = left_hand_side.jacobian(qdd)
..     gd = left_hand_side.xreplace(qdd_zerod)

..     F_r = sm.Matrix([0.0, T_b, T_c])
..     qdd_sol = Md.solve(F_r - gd)


.. .. Practice problem: add a damping force or a coulomb friction force in the first joint 
.. .. (the example and this problem are inspired by a talk by A. Ruina, https://www.youtube.com/watch?v=j-wHI764dWU)


.. The generalized momenta are an invertable function of the generalized speeds. The Euler-Lagrange
.. equation can therefore be rewritten in the form:

.. .. math::

..     \dot{p_r} = \frac{\partial L}{\partial q_r} + \bar{F}_r

.. .. math::

..     \dot{q_r} = \dot{q_r}(\bar{p})  

.. which forms a `Hamiltonian System`_. Hamiltonian systems and their
.. extension Port-Hamiltonian systems are often used in physics and control theory respectively.

.. .. _`Hamiltonian System`: https://en.wikipedia.org/wiki/Hamiltonian_system

.. For the system described above, the following code derives these equations:

.. .. jupyter-execute::

..    p1, p2, p3 = me.dynamicsymbols('p1, p2, p3')
..    p_sym = sm.Matrix([p1, p2, p3])
..    qd_repl = sm.solve(p_sym - p.transpose(), qd)
..    pd = F_r - L.jacobian(q).transpose()
..    qd_solve = qd.xreplace(qd_repl)

.. There are two important realizations:

.. .. jupyter-execute::

..    pd

.. The time derivative of the first generalized momentum is zero. That means the generalized momentum
.. is conserved. This is always the case when the Lagrangian does not depend on a given generalized coordinate, and there
.. are no non-conservative active forces acting on that coordinate either. This statement is a particular case of
.. `Noether's theorem`_.

.. .. _`Noether's theorem`: https://en.wikipedia.org/wiki/Noether%27s_theorem_

.. .. jupyter-execute::

..    p.transpose().jacobian(qd) - Md

.. The Jacobian of the generalized momenta with respect to the generalized velocities is the mass matrix. This is always
.. true, because the kinetic energy can be written as :math:`\frac{1}{2}\dot{\bar{q}}^\text{T}\mathbf{M}_d\dot{\bar{q}}`. 
.. As a result

.. .. math::

..     \bar{p} = \mathbf{M}_d(q)\dot{\bar{q}},

.. which explains the name generalized momentum, as this matches the definitions of momentum and angular momentum in the case
.. of pointmasses.


.. (Learn more) Euler-Lagrange in optimization
.. ===========================================

.. The Euler-Lagrange equation also appears in a different setting: optimization. When optimizing
.. a function :math:`f` over its arguments :math:`q`, we have the well known necessary condition for an optimum:

.. .. math::

..     \frac{\partial f}{\partial q} = 0

.. It is also possible to consider optimizing not over variables, but over functions of one variable. This problem
.. is considered in the mathematical field `Calculus of Variations`_
.. To do so, there must then be a function-like thing that turns possible function into a value which we want to
.. optimize. Such a function-like thing is called a functional, and is often given as an integral. The
.. optimization problem then takes the following form:

.. .. _`Calculus of Variations`: https://en.wikipedia.org/wiki/Calculus_of_variations

.. .. math::

..     \min_{q(t)} \int_{0}^{T} L(t, q, \dot{q})\text{d}t \quad \text{subject to} \quad q(0) = 0, q(T) = q_T  

.. Examples of such optimizations are:

.. * The shortest path problem, where :math:`L = |\dot{q}|`
.. * The brachistochrone problem, that tries to find the shape of a slope, such that a ball rolling off it
..   reaches the bottom in minimal time
.. * Various optimal control problem, in which the integral over the torque squared plus the position error squared
..   should be minimized.

.. For the functional optimization problem, there is again a necessary condition:

.. .. math::

..     \frac{\text{d}}{\text{d}t}\frac{\partial L}{\partial \dot{q}} - \frac{\partial L}{\partial q}= 0,

.. which we recognize as the Euler-Lagrange equations.

.. This means that the laws of nature governing rigid body motions result in motions that minimize the integral of the
.. Lagrangian.  This is called Hamilton's principle. It turns out that 
.. `many physical laws_` take such a form of minimizing
.. the value of a function. One example is Fermat's principle, which states that light takes the path of minimum time.

.. .. _`many physical laws`: https://en.wikipedia.org/wiki/Variational_principle

.. The optimization point-of-view of the Lagrange method also gives an interpretation for the Lagrange multipliers. They
.. are the same as the Lagrange multipliers used in optimization.






