======================
Holonomic  Constraints
======================

.. note::

   You can download this example as a Python script:
   :jupyter-download-script:`configuration` or Jupyter Notebook:
   :jupyter-download-notebook:`configuration`.

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

- derive and specify the configuration constraints (holonomic constraints)
  equations for a system of connected rigid bodies
- numerically solve a set of holonomic constraints for the dependent
  coordinates
- apply a point configuration constraints as a general approach to constraining
  a system
- calculate the number of generalized coordinates
- choose generalized coordinates
- calculate velocities when holonomic constraints are present

Four-Bar Linkage
================

.. todo:: This chapter would be clearer if I made a P_4* and P_4 so we can talk
   about the points that come into existence when the linkage is cut.

Consider the linkage shown below:

.. _configuration-four-bar:
.. figure:: figures/configuration-four-bar.svg
   :align: center
   :width: 600px

   a) Shows four links in a plane :math:`A`, :math:`B`, :math:`C`, and
   :math:`N` with respective lengths :math:`l_a,l_b,l_c,l_n` connected in a
   closed loop at points :math:`P_1,P_2,P_3,P_4`. b) Shows the same linkage
   that has been separated at point :math:`P_4` to make it an open chain of
   links.

This is a planar `four-bar linkage`_ with reference frames :math:`N,A,B,C`
attached to each bar. Four bar linkages are used in a wide variety of
mechanisms. One you may be familiar with is this rear suspension on a mountain
bicycle:

.. _mountain-bike-suspension:
.. figure:: https://upload.wikimedia.org/wikipedia/commons/thumb/7/7c/MtbFrameGeometry_FSR.png/320px-MtbFrameGeometry_FSR.png
   :align: center

   Four bar linkage shown in blue, red, orange, and green used in the rear
   suspension mechanism of a mountain bicycle.

   Cartemere, CC BY-SA 3.0 https://creativecommons.org/licenses/by-sa/3.0, via Wikimedia Commons

.. _four-bar linkage: https://en.wikipedia.org/wiki/Four-bar_linkage

Depending on the length of the links, different motion types are possible.
:numref:`grashof-animation` shows some of the possible motions.

.. _grashof-animation:
.. figure:: https://upload.wikimedia.org/wikipedia/commons/c/ca/Grashof_Type_I_Four-Bar_Kinematic_Inversions.gif
   :align: center
   :width: 80%

   Pasimi, CC BY-SA 4.0 https://creativecommons.org/licenses/by-sa/4.0, via Wikimedia Commons

A four bar linkage is an example of a *closed kinematic loop*. The case of
:numref:`configuration-four-bar` there are two vector paths to point
:math:`P_4` from :math:`P_1`:

.. math::
   :label: vector-loop

   \bar{r}^{P_4/P_1} & = l_n \hat{n}_x \\
   \bar{r}^{P_4/P_1} & = \bar{r}^{P_2/P_1} + \bar{r}^{P_3/P_2} + \bar{r}^{P_4/P_3} = l_a\hat{a}_x + l_b\hat{b}_x + l_c\hat{c}_x

For the loop to close, the two vector paths must equate. We can resolve this by
disconnecting the loop at some location, :math:`P_4` in our case, and forming
the *open loop* vector equations to points that should coincide. Keep in mind
that we assume that the lengths are constant and the angles change with time.

Setup the variables, reference frames, and points:

.. jupyter-execute::

   q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')
   la, lb, lc, ln = sm.symbols('l_a, l_b, l_c, l_n')

   N = me.ReferenceFrame('N')
   A = me.ReferenceFrame('A')
   B = me.ReferenceFrame('B')
   C = me.ReferenceFrame('C')

   A.orient_axis(N, q1, N.z)
   B.orient_axis(A, q2, A.z)
   C.orient_axis(B, q3, B.z)

   P1 = me.Point('P1')
   P2 = me.Point('P2')
   P3 = me.Point('P3')
   P4 = me.Point('P4')

SymPy Mechanics will warn you if you try to establish a closed loop among a set
of points and you should not do that because functions that use points have no
way to know which vector path you desire to use. Instead you will establish
positions among points on one open leg of the chain:

.. jupyter-execute::

   P2.set_pos(P1, la*A.x)
   P3.set_pos(P2, lb*B.x)
   P4.set_pos(P3, lc*C.x)

   P4.pos_from(P1)

Now, declare a vector for the other path to :math:`P_4`:

.. jupyter-execute::

   r_P1_P4 = ln*N.x

With both vector paths written, we can form the left hand side of the following
equation:

.. math::
   :label: constraint-expression

   \bar{r}^{P_4/P_1} - \left( \bar{r}^{P_2/P_1} + \bar{r}^{P_3/P_2} + \bar{r}^{P_4/P_3} \right) = 0

Using :external:py:meth:`~sympy.physics.vector.point.Point.pos_from` for the
open loop leg made of points and the additional vector:

.. jupyter-execute::

   loop = P4.pos_from(P1) - r_P1_P4
   loop

This "loop" vector expression must equate to zero for our linkage to always be
a closed loop. We have a planar mechanism, so we can extract two scalar
equations associated with a pair of unit vectors in the plane of the mechanism.
We can pick any two non-parallel unit vectors to express the components in, with
the intuitive choice being :math:`\hat{n}_x` and :math:`\hat{y}`.

.. jupyter-execute::

   fhx = sm.trigsimp(loop.dot(N.x))
   fhx

.. jupyter-execute::

   fhy = sm.trigsimp(loop.dot(N.y))
   fhy

For the loop to close, these two expressions must equal zero for all values
:math:`q_1,q_2,q_3`. These are two nonlinear equations in three time varying
variables called coordinates. The solution can be found if we solve for two of
the time varying variables. For example, :math:`q_2` and :math:`q_3` can be
solved for in terms of :math:`q_1`. We would then say that :math:`q_2` and
:math:`q_3` depend on :math:`q_1`. These two equations are called holonomic
constraints, or configuration constraints, because they constrain the kinematic
configuration to be a loop. Holonomic constraints take the form of a real
valued vector function:

.. math::
   :label: configuration-constraint

   \bar{f}_h(q_1, \ldots, q_N, t) = 0 \textrm{ where } \bar{f}_h \in \mathbb{R}^M

:math:`N` is number of coordinates that you have used to describe the system
and :math:`M` is the number of scalar constraint equations.

.. warning::

   Holonomic constraints are defined strictly as equations that are function of
   the :math:`N` time varying coordinates. It is true that these equations are
   only valid for a limited set of ranges for the constants in the equations,
   e.g. the lengths of the bars, but the range and combination constraints on
   the constants are not what we are considering here. Secondly, Eq.
   :math:numref:`configuration-constraint` does not represent inequality
   constraints. A coordinate may be constrained to a specific range, e.g.
   :math:`-\pi<q_1<\pi`, but these are not holonomic constraints in the sense
   definied here. Inequality constraints are generally dealt with using
   collision models to capture the real dynamics of forcefully limiting motion.

The four-bar linkage constraints are functions of configuration variables: time
varying angles and distances. In our case the constraint equations are:

.. math::
   :label: four-bar-constraints

   \bar{f}_h(q_1, q_2, q_3) = 0 \textrm{ where } \bar{f}_h \in \mathbb{R}^2

and :math:`N=3` and :math:`M=2`.

In SymPy, we'll typically form this column vector as so:

.. jupyter-execute::

   fh = sm.Matrix([fhx, fhy])
   fh

.. admonition:: Exercise

   `Watt's Linkage`_ is a four-bar linkage that can generate almost straight
   line motion of the center point of the middle coupler link. Write the
   holonomic constraints for the Watt's Linkage. The coupler link has a length
   of :math:`2a`, the left and right links have length :math:`b`. Make the
   vertical distance between the fixed points of the left and right lengths
   :math:`2a` and the horizontal distance :math:`(2-1/20)b`. Use the same
   reference frame and angle definitions as the four-bar linkage above.

   .. figure:: https://upload.wikimedia.org/wikipedia/commons/9/9e/Watts_Linkage.gif
      :width: 60%
      :align: center

      Arglin Kampling, CC BY-SA 4.0 https://creativecommons.org/licenses/by-sa/4.0, via Wikimedia Commons

   .. _Watt's Linkage: https://en.wikipedia.org/wiki/Watt%27s_linkage

.. admonition:: Solution
   :class: dropdown

   .. jupyter-execute::

      q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')
      a, b = sm.symbols('a, b')

      N = me.ReferenceFrame('N')
      A = me.ReferenceFrame('A')
      B = me.ReferenceFrame('B')
      C = me.ReferenceFrame('C')

      A.orient_axis(N, q1, N.z)
      B.orient_axis(A, q2, A.z)
      C.orient_axis(B, q3, B.z)

      P1 = me.Point('P1')
      P2 = me.Point('P2')
      P3 = me.Point('P3')
      P4 = me.Point('P4')

      P2.set_pos(P1, b*A.x)
      P3.set_pos(P2, 2*a*B.x)
      P4.set_pos(P3, b*C.x)

      P4.pos_from(P1)

      r_P1_P4 = (2 - sm.S(1)/20)*b*N.x - 2*a*N.y

      loop = P4.pos_from(P1) - r_P1_P4

      fh_watts = sm.trigsimp(sm.Matrix([loop.dot(N.x), loop.dot(N.y)]))
      fh_watts

Solving Holonomic Constraints
=============================

Only the simplest of holonomic constraint equations may be solved symbolically
due to their nonlinear nature, so you will in general need to solve them
numerically. In :ref:`Equations of Motion with Holonomic Constraints` we will
show how to solve them for simulation purposes, but for now SymPy's
:external:py:func:`~sympy.solvers.solvers.nsolve` can be used to numerically solve
the equations. If we choose :math:`q_2` and :math:`q_3` to be the dependent
coordinates, we need to select numerical values for all other variables. Note
that not all link length combinations result in a valid linkage geometry.
Starting with the replacements,

.. jupyter-execute::

   import math  # provides pi as a float

   repl = {
       la: 1.0,
       lb: 4.0,
       lc: 3.0,
       ln: 5.0,
       q1: 30.0/180.0*math.pi,  # 30 degrees in radians
   }
   repl

we can then formulate the constraint equations such that only :math:`q_2` and
:math:`q_3` are variables:

.. jupyter-execute::

   fh.xreplace(repl)

Generally, there may be multiple numerical solutions for the unknowns and the
underlying algorithms require a guess to return a specific result. If we make an
educated guess for the unknowns, then we can find the specific solution with
``nsolve()``:

.. jupyter-execute::

   q2_guess = -75.0/180.0*math.pi  # -75 degrees in radians
   q3_guess = 100.0/180.0*math.pi  # 100 degrees in radians

   sol = sm.nsolve(fh.xreplace(repl), (q2, q3), (q2_guess, q3_guess))
   sol/math.pi*180.0  # to degrees

.. admonition:: Exercise

   Find the angles of the remaining links in `Watt's Linkage`_ if the middle
   linkage is rotated clockwise 5 degrees, :math:`a=1`, and :math:`b=4`.

   .. _Watt's Linkage: https://en.wikipedia.org/wiki/Watt%27s_linkage

.. admonition:: Solution
   :class: dropdown

   The angle relative to vertical of the middle link is
   :math:`3\pi/2-(q_1+q_2)`, which we can use to solve for :math:`q_2`.

   .. jupyter-execute::

      repl = {
          a: 1.0,
          b: 4.0,
          q2: 3.0*math.pi/2.0 - 5.0/180.0*math.pi - q1,
      }
      repl

   .. jupyter-execute::

      fh_watts.xreplace(repl)

   .. jupyter-execute::

      q1_guess = 10.0/180.0*math.pi
      q3_guess = 100.0/180.0*math.pi

      sol = sm.nsolve(fh_watts.xreplace(repl), (q1, q3), (q1_guess, q3_guess))
      sol/math.pi*180.0  # to degrees

..
   .. jupyter-execute::

      # code to plot the linkage
      coordinates = P1.pos_from(P1).to_matrix(N)
      for point in [P2, P3, P4]:
          coordinates = coordinates.row_join(point.pos_from(P1).to_matrix(N))
      eval_point_coords = sm.lambdify((q1, q2, q3, a, b), coordinates)
      eval_point_coords(1.0, 2.0, 3.0, 4.0, 5.0)
      x, y, _ = eval_point_coords(
          float(sol[0, 0]),
          float(repl[q2].xreplace({q1: sol[0, 0]})),
          float(sol[1, 0]),
          repl[a], repl[b])
      import matplotlib.pyplot as plt
      plt.plot(x, y)
      plt.grid()
      plt.axis('equal')

General Holonomic Constraints
=============================

If you consider a set of :math:`v` points, :math:`P_1,P_2,\ldots,P_v` that can
move unconstrained in Euclidean 3D space, then one would need :math:`3v`
constraint equations to fix the points (fully constrain the motion) in that
Euclidean space. For the four points in the four-bar linkage, we would then
need :math:`3(4)=12` constraints to lock all the points fully in place. The
figure below will be used to illustrate the general idea of constraining the
configuration of the four bar linkage.

.. _configuration-constraints:
.. figure:: figures/configuration-constraints.svg
   :align: center
   :width: 400px

   a) Four points in 3D space, b) four points constrained to 2D space, c)
   points are fixed to adjacent points by a fixed length, d) the first point is
   fixed at :math:`O` in two dimensions, e) the fourth point is fixed in the
   :math:`y` coordinate relative to :math:`O`.

Starting with a), there are the four points in 3D Euclidean space that are free
to move. Moving to b), each of the four points can be then constrained to be in
a plane with:

.. math::
   :label: planar-constraints

   \bar{r}^{P_1/O}\cdot\hat{n}_z = 0 \\
   \bar{r}^{P_2/O}\cdot\hat{n}_z = 0 \\
   \bar{r}^{P_3/O}\cdot\hat{n}_z = 0 \\
   \bar{r}^{P_4/O}\cdot\hat{n}_z = 0

where :math:`O` is a point fixed in :math:`N`. This applies four constraints
leaving 8 coordinates for the planar location of the points. Now at c) we
constrain the points with:

.. math::
   :label: length-constraints

   |\bar{r}^{P_2/P_1}| = l_a \\
   |\bar{r}^{P_3/P_2}| = l_b \\
   |\bar{r}^{P_4/P_3}| = l_c \\
   |\bar{r}^{P_4/P_1}| = l_n

These four constraint equations keep the points within the specified distances
from each other leaving 4 coordinates free. In d) point :math:`P_1` is fixed
relative to :math:`O` with 2 scalar constraints:

.. math::
   :label: p1-constraint

   \bar{r}^{P_1/O}\cdot\hat{n}_x = 0 \\
   \bar{r}^{P_1/O}\cdot\hat{n}_y = 0

Finally in e), :math:`P_4` is constrained with the single scalar:

.. math::
   :label: p4-constraint

   \bar{r}^{P_4/P_1} \cdot \hat{n}_y = 0

Notice that we did not need :math:`\bar{r}^{P_4/P_1} \cdot \hat{n}_x = 0`,
because :math:numref:`length-constraints` ensures the :math:`x` coordinate of
:math:`P_4` is in the correct location.

These 11 constraints leave a single free coordinate to describe the orientation
of :math:`A`, :math:`B`, and :math:`C` in :math:`N`. When we originally
sketched :numref:`configuration-four-bar` most of these constraints were
implied, i.e. we drew a planar mechanism with points :math:`P_1` and
:math:`P_4` fixed in :math:`N`, but formally there are 12 coordinates needed to
locate the four points and 11 constraints that constrain them to have the
configuration of a four-bar linkage.

A general holonomic constraint for a set of :math:`v` points with Cartesian
coordinates is then ([Kane1985]_ pg. 35):

.. math::
   :label: holonomic-cartesian

   f_h(x_1, y_1, z_1, \ldots, x_v, y_v, z_v, t) = 0

We include :math:`t` as it is possible that the constraint is an explicit
function of time (instead of only implicit, as seen above in the four-bar
linkage example).

Generalized Coordinates
=======================

If a set of :math:`v` points are constrained with :math:`M` holonomic
constraints then only :math:`n` of the Cartesian coordinates are independent of
each other. The number of independent coordinates is then defined as
([Kane1985]_ pg. 37):

.. math::
   :label: num-gen-coord

   n := 3v - M

These :math:`n` independent Cartesian coordinates can also be expressed as
:math:`n` functions of time :math:`q_1(t),q_2(t),\ldots,q_n(t)` in such a way
that the constraint equations are always satisfied. These functions
:math:`q_1(t),q_2(t),\ldots,q_n(t)` are called *generalized coordinates* and it
is possible to find :math:`n` independent coordinates that minimize the number
of explicit constraint equations needed to describe the system's configuration
at all times :math:`t`. These generalized coordinates are typically determined
by inspection of the system and there is a bit of an art to choosing the best
set. But you can always fall back to the formal process of constraining each
relevant point. If you describe your system with :math:`N\leq3v` coordinates
then:

.. math::

   n := N - M

Take this simple pendulum with points :math:`O` and :math:`P` as an example:

.. figure:: figures/configuration-pendulum.svg
   :align: center
   :width: 400px

If the pendulum length :math:`l` is constant and the orientation between
:math:`A` and :math:`N` can change, then the location of :math:`P` relative to
:math:`O` can be described with the Cartesian coordinates :math:`x` and
:math:`y`. It should be clear that :math:`x` and :math:`y` depend on each other
for this system. The constraint relationship between those two coordinates is:

.. math::
   :label: pendulum-length-constraint

   x^2 + y^2 = l^2

This implies that only one coordinate is independent, i.e. :math:`n=1`. More
formally, the two points give :math:`3v=3(2)=6` and there are 2 constraints for
the planar motion of each point, 2 constraints fixing :math:`O` in :math:`N`,
and 1 constraint fixing the distance from :math:`O` to :math:`P`, making
:math:`M=5` and thus confirming our intuition :math:`n=6-5=1`.

But there may be functions of time that relieve us from having to consider Eq.
:math:numref:`pendulum-length-constraint`. For example, these two coordinates
can also be written as as functions of the angle :math:`q`:

.. math::
   :label: xy-func-of-q

   x = l\cos q \\
   y = l\sin q

and if we describe the configuration with only :math:`q`, the constraint is
implicitly satisfied. :math:`q` is then a generalized coordinate because it
satisfies :math:`n=1` and we do not have to explicitly define a constraint
equation.

Now, let's return to the four-bar linkage example in
:numref:`configuration-four-bar` and think about what the generalized
coordinates of this system are. We know, at least intuitively, that :math:`n=1`
for the four bar linkage. The four-bar linkage in
:numref:`configuration-four-bar` is described in a way that assumes a number of
constraints are fulfilled, such as Eqs. :math:numref:`planar-constraints` and
:math:numref:`p1-constraint`, so we do not have to formally consider them.

.. admonition:: Exercise

   Are :math:`q_1,q_2,q_3` generalized coordinates of the four-bar linkage? If
   not, why?

.. admonition:: Solution
   :class: dropdown

   Any one of the :math:`q_1,q_2,q_3` can be a generalized coordinate, but only
   one. The other two are depdendent due to the two constraints. We started
   with three coordinates :math:`q_1,q_2,q_3` describing the open chain
   :math:`P_1` to :math:`P_2` to :math:`P_3` to :math:`P_4`. Then we have two
   scalar constraint equations, leaving :math:`n=1`. Thus we can choose
   :math:`q_1`, :math:`q_2`, **or** :math:`q_3` to be the indepdendent
   generalized coordinate.

If we take the general approach, starting with four unconstrained points, we
need 11 constraints to describe the system, but if we select generalized
coordinates to describe the system we only need 2 constraint equations (Eq.
:math:numref:`four-bar-constraints`)! This simplifies the mathematical problem
description and, as we will later see, is essential for obtaining the simplest
forms of the equations of motion of a multibody system.

Calculating Additional Kinematic Quantities
===========================================

You will often need to calculate velocities and accelerations of points and
reference frames of systems with holonomic constraints. Due to the
differentiation chain rule, velocities will be linear in the time derivatives
of the coordinates and accelerations will be linear in the double time
derivatives of the coordinates. Our holonomic constraints dictate that there is
no relative motion between points or reference frames, implying that the
relevant positions, velocities, and accelerations will all equate to zero.

Start by setting up the points for the four-bar linkage again:

.. jupyter-execute::

   P1 = me.Point('P1')
   P2 = me.Point('P2')
   P3 = me.Point('P3')
   P4 = me.Point('P4')
   P2.set_pos(P1, la*A.x)
   P3.set_pos(P2, lb*B.x)
   P4.set_pos(P3, lc*C.x)

In the four-bar linkage, :math:`{}^N\bar{v}^{P_4}` must be zero. We can
calculate the unconstrained velocity like so:

.. jupyter-execute::

   P1.set_vel(N, 0)
   P4.vel(N)

The scalar velocity constraints can be formed in a similar fashion as the
configuration constraints:

.. math::

   {}^N\bar{v}^{P_4}\cdot\hat{n}_x = 0 \\
   {}^N\bar{v}^{P_4}\cdot\hat{n}_y = 0

.. jupyter-execute::

   sm.trigsimp(P4.vel(N).dot(N.x))

.. jupyter-execute::

   sm.trigsimp(P4.vel(N).dot(N.y))

Notice that this is identical to taking the time derivative of the constraint
vector function :math:`\bar{f}_h`:

.. jupyter-execute::

   t = me.dynamicsymbols._t
   fhd = fh.diff(t)
   fhd

We can see that the expressions are linear in :math:`\dot{q}_1,\dot{q}_2` and
:math:`\dot{q}_3`. If we select :math:`\dot{q}_2` and :math:`\dot{q}_3` to be
dependent, we can solve the linear system :math:`\mathbf{A}\bar{x}=\bar{b}` for
those variables using the technique shown in :ref:`Solving Linear Systems`.
First we define a column vector holding the dependent variables:

.. jupyter-execute::

   x = sm.Matrix([q2.diff(t), q3.diff(t)])
   x

then extract the linear terms:

.. jupyter-execute::

   A = fhd.jacobian(x)
   A

find the terms not linear in the dependent variables:

.. jupyter-execute::

   b = -fhd.xreplace({q2.diff(t): 0, q3.diff(t): 0})
   b

and finally solve for the dependent variables:

.. jupyter-execute::

   x_sol = sm.simplify(A.LUsolve(b))
   x_sol

Now we can write any velocity strictly in terms of the independent speed
:math:`\dot{q}_1` and all of the other coordinates.
:external:py:meth:`~sympy.physics.vector.vector.Vector.free_dynamicsymbols`
shows us what coordinates and their time derivatives present an any vector:

.. jupyter-execute::

   P4.vel(N).free_dynamicsymbols(N)

Using the results in ``x_sol`` above we can write the velocity in terms of only
the independent :math:`\dot{q}_1`:

.. math::

   {}^N\bar{v}^A =
   v_x(\dot{q}_1, q_1, q_2, q_3)\hat{n}_x +
   v_y(\dot{q}_1, q_1, q_2, q_3)\hat{n}_y +
   v_z(\dot{q}_1, q_1, q_2, q_3)\hat{n}_z

Making the substitutions gives the desired result:

.. jupyter-execute::

   qd_dep_repl = {
     q2.diff(t): x_sol[0, 0],
     q3.diff(t): x_sol[1, 0],
   }
   qd_dep_repl

.. jupyter-execute::

   P4.vel(N).xreplace(qd_dep_repl)

.. jupyter-execute::

   P4.vel(N).xreplace(qd_dep_repl).free_dynamicsymbols(N)

The holonomic constraints will have to be solved numerically as described in
:ref:`Solving Holonomic Constraints`, but once done only the independent
:math:`\dot{q}_1` is needed.

.. todo:: Add exercise to numerically calculate the velocity of the center
   point of the Watt's Linkage.
