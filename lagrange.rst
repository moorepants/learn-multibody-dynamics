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

Introduction
============

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

As an alternative, the following equation also calculates the generalized inertial forces of a
system, now by starting from the `kinetic energy TODO`_ :math:`T(\dot{\bar{q}}, \bar{q})`
expressed as function of the generalized coordinates :math:`\bar{q}`, and 
their time derivatives.

.. _`kinetic energy TODO`: https://en.wikipedia.org/wiki/Work_in_process


.. math::
  :label: eq-lagrange-inertial

   -\bar{F}^*_r = \frac{\mathrm{d}}{\mathrm{d}t}\left(\frac{\partial T}{\partial u}
        \right) - \frac{\partial T}{\partial q_r}

.. warning:: Note the two minus signs in the above equation

.. note::

   In Kane's method, we are free to choose generalized speeds :math:`\bar{u}` that differ from
   the time derivatives of the generalized coordinates, that is :math:`\dot{\bar{q}}`. This
   freedom is not readily there when using the Lagrange method. Therefore we will consistently use
   :math:`\dot{\bar{q}}` in this chapter.

The generalized inertial forces computed in this manner are the same as when following
the TMT method. This can be shown by carefully matching terms in both formulations, as
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


The Lagrange Method
===================

Both the equation for computing the inertial forces from the kinetic energy, and 
the equation for computing the applied forces from a potential energy have a term
in them with the partial derivative with respect to the generalized coordinate. 
Furthermore, the potential energy does not depend on the generalized speeds. 
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


Example: 2D double pendulum with springs and sliding pointmass

Because further examples include multiple bodies, we introduce two convenience functions to
simplify the code for computing the kinetic energy:

.. jupyter-execute::

   def squarednorm(a):
       return a.dot(a)

   def quadraticform(I, v):
       return v.dot(I.dot(v))

We can then go on to define the relevant variables, constants and frames:

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

Finally, we setup the Lagrangian and derive the equations of motion:

.. jupyter-execute::

   t = sm.symbols('t')
   q = sm.Matrix([q1, q2, q3])
   qdot = q.diff(t)
   qddot = qdot.diff(t)

   T = m/2*(squarednorm(Ao.vel(N)) + squarednorm(Bo.vel(N)) + squarednorm(Q.vel(N))) + 1/2*(
       quadraticform(I_A_Ao, A.ang_vel_in(N)) + quadraticform(I_B_Bo, A.ang_vel_in(N))
   )
   V = m*g*(Ao.pos_from(O).dot(-N.x) + Bo.pos_from(O).dot(-N.x)) + kt/2*(q1**2) + kt/2*q2**2 + kl/2*q3**2

   L = sm.Matrix([T - V])
   lhs = L.jacobian(qdot).diff(t) - L.jacobian(q)
   M = lhs.transpose().jacobian(qddot)
   G = lhs.transpose() - M*qddot


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

For a particle of mass :math:`m` and position :math:`\bar{r}^{P/O} = q_1 \hat{n}_x + q_2 \hat{n}_y + q_3\hat{n}_z` on a 
slope :math:`q_1 = q_2`, we take the second time derivative of the constraint equation :math:`\ddot{q_1} - \ddot{q_2} = 0`,
and a constraint force :math:`\bar{F} = F\hat{n}_x - F\hat{n}_y`. Using the Lagrangian 
:math:`V = \frac{1}{2}m(\dot{q}_1^2 + \dot{q}_2^2 + \dot{q}_3^2) - mgq_3`, we can then derive:

.. math::
    \begin{array}{r}
    m\ddot{q}_1= \phantom{-}F\\
    m\ddot{q}_2= -F\\
    m\ddot{q}_3 + mg = \phantom{-}0\\
    \ddot{q}_1 - \ddot{q}_2\!\! = \phantom{-}0 
    \end{array}

This can be put in matrix-form, by extracting the unknown acceleration and force magnitude;

.. math::
    \begin{bmatrix} m & 0 & 0 &-1 \\ 0 & m & 0 & 1 \\ 0 & 0 & m & 0 \\ 1 & -1 & 0 & 0\end{bmatrix}
    \begin{bmatrix} \ddot{q}_1 \\ \ddot{q}_2 \\ \ddot{q}_3 \\ F \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \\ -mg \\ 0\end{bmatrix}


It can be difficult to find the direction of the constraint force from the geometric of the system directly.
There is a trick, called the method of the Lagrange multupliers, to quickly find the correct generalized
forces associated with the constraint forces. 

Given a constraint in the general form

.. math::

    \sum_r a_r(\bar{q}) \dot{q}_r = 0

We find the generalized force as:

.. math::

    F_r = \lambda a_r(\bar{q})

Here :math:`\lambda` is a variable encoding the magnitude of the constraint force. It is
called  the Lagrange multiplier. The same :math:`\lambda` is used for each :math:`r`, that is, 
each constraint has a single associated Lagrange multiplier.

Due to how it is constructed, the power produced by the constraint force is always zero, as expected.

.. math::

    P = \sum_r F_r\dot{q}_r = \sum \lambda a_r\dot{q}_r  = \lambda \sum a_r\dot{q}_r = \lambda \cdot 0

For example, consider the pointmass to be constrained to move in a bowl :math:`x^2 + y^2 + z^2 -1 = 0`.
Taking the time derivative, we find: :math:`a_1 = 2q_1` :math:`a_2 = 2q_2` and :math:`a_3 = 2q_3`.
We would find :math:`F_1 = 2\lambda q_1`, :math:`F_2 = 2\lambda q_2` and :math:`F_3 = 2\lambda q_3`.


**Example: turning the freely floating body discussed earlier into a rolling sphere.**

The non-slip condition of the rolling sphere is split into three constraints: the velocity of
the contact point (:math:`G`) is zero in both the :math:`\hat{n}_x`, :math:`\hat{n}_y` and :math:`\hat{n}_z`
direction. These constraints are enforced by contact forces in their respective directions.

The contact point can be found according by :math:`\bar{r}^{G/C} = -r \hat{n}_z`. We therefore get the
constraint:

.. math::

    \begin{array}{l}
    \bar{n}_x\cdot ({}^N\bar{v}^C + {}^N\bar{\omega}^B \times -r\hat{n}_z) = 0 \\
    \bar{n}_y\cdot ({}^N\bar{v}^C + {}^N\bar{\omega}^B \times -r\hat{n}_z) = 0 \\
    \bar{n}_z\cdot ({}^N\bar{v}^C + {}^N\bar{\omega}^B \times -r\hat{n}_z) = 0 \\
    \end{array}

These can be used to derive the constraint force and the additional equations using the Lagrange-multiplier
method as shown below. Note that here only the first time derivative of the constraint equation is used, 
again because the second time derivatives of the generalized coordinates appear.

.. container:: invisible

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
        mass_matrix = left_hand_side.jacobian(qddot)
        dynamic_bias = left_hand_side - mass_matrix*qddot

To make these free floating body a rolling wheel, three constraints are needed: the
velocity of the contact point should be zero in :math:`\hat{n}_x`, :math:`\hat{n}_y`
and :math:`\hat{n}_x` direction.

.. jupyter-execute::

    lambda1, lambda2, lambda3 = me.dynamicsymbols('lambda1, lambda2, lambda3') 
    constraint = (v_com + B.ang_vel_in(N).cross(-N.z)).to_matrix(N)
    A = constraint.jacobian(qdot)
    diff_constraint = constraint.diff(t)

This constraint information must then be added to the original equations. To do
so, we make use of a useful fact, which is true for all systems and constraints:

.. jupyter-execute::

    diff_constraint.jacobian(qddot) - A

This allows us to create our equations in a block matrix form:

.. math::
        \begin{bmatrix} M & A^T \\ A & 0\end{bmatrix}\begin{bmatrix}\ddot{\bar{q}} \\ \lambda \end{bmatrix} = 
        \begin{bmatrix} F_r - g \\ - \frac{\partial A\dot{\bar{q}}}{\partial \bar{q}}\dot{\bar{q}} \end{bmatrix},

where :math:`A` is the jacobian of the constraints, as used above,  :math:`g` is the dynamic bias, and the last term on the right hand side can be computed as;

.. jupyter-execute::

    constraint_bias = diff_constraint - diff_constraint.jacobian(qddot)*qddot





    
    

    



    
The method of the Lagrange multiplier can of course also be used within Kane's method. However,
this results in a larger system of equations, which is why the elimination approach is often
preferred there. An exception being scenarios where the constraint force itself is a useful output,
for instance to check no-slip conditions in case of limited friction.


Lagrange's vs Kane's
====================

The Lagrangian method is the  second method to deriving the equations of motion presented in this book,
after Kane's method. This raises the questions: when should each
method be used.

For constrained systems, Kane's method has the advantage that the equations of motion are given for a set of
independent generalized velocities only. This can give rise to simplified equations, additional insight, and
numerically more efficient simulation.

Furthermore, the connection from Kane's method to vector mechanics, that is, Newton's law's, is clearer, which
can provide additional insight, and make it easier to encorporate non-conservative forces such as friction.

On the other hand, the Lagrange-method results in a set of equations with well understood structures and properties.
These structures and properties are not studied further in these materials, other than the following "learn more" section.
With further study, these aspects can make it easier to generalize results across multiple systems, for example
when designing control algorithms.


(Learn more) Generalized momentum
=================================

The partial derivative of the Lagrangian with respect to generalized speed is
called the generalized momentum.

Examples showing that this matches to momentum and angular momentum in relevant 
particle cases.

If the Lagrangian does not depend on a generalized coordinates, its associated
generalized momentum is conserved.

Some ideas behind generalized momentum will be discussed with the following example,
which is a simplified version of the falling cat example:
* body A is a cylinder that can rotate wrt ground around same axis as gravity: :math:`\hat{n}_z``
* body B is a cylinder that can rotate wrt body A around same axis as gravity
* body C is a cylinder that can rotate wrt body C around a (body fixed) axis perpendicular to gravity :math:`\hat{b}_x`
* There are two actuators providing a torque on the joints between bodies A and B and bodies B and C respectively.

This example will also show how to apply motor torques at joints.

.. jupyter-execute::

   t, l, r, T_b, T_c = sm.symbols('t, l, r, T_b, T_c')
   q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')

   N = me.ReferenceFrame('N')
   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')
   C = me.ReferenceFrame('C')

   A.orient_axis(N, q1, N.z)
   B.orient_axis(A, q2, A.z)
   C.orient_axis(B, q3, B.x) 

   g = 1
   rho = 1
   m = rho*l*sm.pi*r**2
   I_xx_or_yy = m/12*(3*r**2 + l**2)
   I_zz= m/2*r**2
   I_A_Ao = me.inertia(A, I_xx_or_yy , I_xx_or_yy, I_zz)
   I_B_Bo = me.inertia(B, I_xx_or_yy , I_xx_or_yy, I_zz)
   I_C_Co = me.inertia(C, I_xx_or_yy , I_xx_or_yy, I_zz)

   O = me.Point('O')
   O.set_vel(N, 0.0)
   Ao = me.Point("A_c")
   Ao.set_pos(O, -0.5*l*A.z)
   Bo = me.Point("B_c")
   Bo.set_pos(Ao, -0.5*l*A.z - 0.5*l*B.z)
   Co = me.Point("C_c")
   Co.set_pos(Bo, -0.5*l*B.z -0.5*l*C.z)

The next step is again to form the Lagrangian and find the equations of motion. As the system has no further constraints, 
the Lagrange multiplier method is not needed. The actuator torques are added to the right hand side of the equation, in
the same way as active forces are added to Kane's equations. Here the torques are represented by the variables :math:`T_b`
and :math:`T_c` are used to represent.

.. jupyter-execute::

   T = m/2*(squarednorm(Ao.vel(N)) + squarednorm(Bo.vel(N)) + squarednorm(Co.vel(N))) + 1/2*(
           quadraticform(I_A_Ao, A.ang_vel_in(N)) + quadraticform(I_B_Bo, B.ang_vel_in(N)) + quadraticform(I_C_Co, C.ang_vel_in(N)))
   V = m*g*N.z.dot(Co.pos_from(O))
   L = sm.Matrix([T - V])

   q = sm.Matrix([q1, q2, q3])
   q_dot = q.diff(t)
   q_ddot = q_dot.diff(t)

   p = L.jacobian(q_dot)
   p.simplify()
   lhs = p.diff(t) - L.jacobian(q)
   rhs = sm.Matrix([0.0, T_b, T_c])

   M = lhs.transpose().jacobian(q_ddot)
   G = lhs.transpose() - M*q_ddot

   q_ddot_sol = M.solve(rhs - G)


.. Practice problem: add a damping force or a coulomb friction force in the first joint 
.. (the example and this problem are inspired by a talk by A. Ruina, https://www.youtube.com/watch?v=j-wHI764dWU)


The generalized momenta are an invertable function of the generalized speeds. We can therefore replace the
Lagrangian equation by:

.. math::

    \dot{p_r} = \frac{\partial L}{\partial q_r}

.. math::

    \dot{q_r} = \dot{q_r}(\bar{p})  

which are equivalent to the equations obtained using Hamilton's method. Hamiltonian systems and their
extension Port-Hamiltonian system are often used in physics and control theory respectively.

For the system described above, the following derives these equations:

.. jupyter-execute::

   p1, p2, p3 = me.dynamicsymbols('p1, p2, p3')
   p_sym = sm.Matrix([p1, p2, p3])
   J_p_wrt_qdot = p.transpose().jacobian(q_dot)
   p_dot = rhs - L.jacobian(q).transpose()
   q_dot_solve = J_p_wrt_qdot.solve(p_sym)

There are two important realizations:

.. jupyter-execute::

   p_dot

Here we see t that the time derivative of the first generalized momentum is zero. That means the generalized momentum
is conserved. This is always the case when the Lagrangian does not depend on a given generalized coordinate, and there
are no non-conservative active forces acting on that coordinate either. This statement is a particular case of the
so called Noether's theorem.

.. jupyter-execute::

   J_p_wrt_qdot - M

The jacobian of the generalized momenta with respect to the generalized coordinates is the mass matrix. This is always
true. As a result, we have:

.. math::

    p = M(q)\dot{q},

which explains the name generalized momentum, as this matches the definitions of momentum and angular momentum in the case
of pointmasses.


(Learn more) Euler-Lagrange in optimization
===========================================

The Euler-Lagrange equation also appears in a different setting: optimization. When optimizing
a function $f$ over its arguments $q$, we have the well known necessary condition for an optimum:

.. math::

    \frac{\partial f}{\partial q} = 0

It is also possible to consider optimizing not over variables, but over functions of one variable. 
To do so, there must then be a function-like thing that turns possible function into a value which we want to
optimize. Such a function-like thing is called a functional, and is often given as an integral. The
optimization problem then takes the following form:

.. math::

    \min_{q(t)} \int_{0}^{T} L(t, q, \dot{q})\text{d}t \quad \text{s.t.} q(0) = 0, q(T) = q_T  

Examples of such optimizations are:

* The shortest path problem, where :math:`L = |\dot{q}|`
* The brachistochrone problem, that tries to find the shape of a slope, such that a ball rolling off it
  reaches the bottom in minimal time
* Various optimal control problem, in which the integral over the torque squared plus the position error squared
  should be minimized.

For the functional optimization problem, there is again a necessary condition:

.. math::

    \frac{\text{d}}{\text{d}t}\frac{\partial L}{\partial \dot{q}} - \frac{\partial L}{\partial q}= 0,

which we recognize as the Euler-Lagrange equations.

This means that the laws of nature governing rigid body motions result in motions that minimize the integral of the
Lagrangian.  This is called Hamilton's principle. It turns out that many physical laws take such a form of minimizing
the value of a function. One example is Fermat's principle, which states that light takes the path of minimum time.

The optimization point-of-view of the Lagrange method also gives an interpretation for the Lagrange multipliers. They
are the same as the Lagrange multipliers used in optimization.






