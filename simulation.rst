==========
Simulation
==========

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`simulation` or Jupyter Notebook:
   :jupyter-download:notebook:`simulation`.

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   me.init_vprinting(use_latex='mathjax')

Numerical Integration
=====================

As mentioned at the end of the prior chapter, we will need to numerically
integrate the equations of motion. If they are in explicit form, this integral
describes how we can arrive at trajectories in time for the state variables by
integrating with respect to time from an initial time :math:`t_0` to a final
time :math:`t_f`.

.. math::
   :label: eq-eom-integral-2

   \bar{x}(t) = \int^{t_f}_{t_0} \bar{f}_m(\bar{x}, t) dt

Recall that :math:`\bar{f}_m` is:

.. math::
   :label: eq-state-form-explicit-2

   \bar{f}_m(\bar{x}, t)
   =
   -\mathbf{M}_m^{-1}
   \bar{g}_m

It is possible to form :math:`-\mathbf{M}_m^{-1}` symbolically and it may be
suitable or preferrable for a given problem, but there are some possible
drawbacks. For example, if the degrees of freedom are quite large, the
resulting symbolic equations become exponetially more complex. Thus, we will
now move from symbolics to numerics.

The NumPy_ library is the defacto base library for numeric computing with
Python. NumPy allows us to do `array programming`_ with Python by providing
floating point array data types and vectorized operators to enable repeat
operations across arrays of values. In Sec.
:ref:`sec-evaluating-symbolic-expressions` we introduced the SymPy function
:external:py:func:`~sympy.utilities.lambdify.lambdify`. ``lambdify()`` will be
our way to bridge the symbolic world of SymPy with the numeric world of NumPy.

.. _NumPy: https://numpy.org
.. _array programming: https://en.wikipedia.org/wiki/Array_programming

We will import NumPy like so:

.. jupyter-execute::

   import numpy as np

.. warning::

   Beware that mixing SymPy and NumPy datatypes will rarely, if at all, provide
   you with functioning code. Be careful because sometimes it may look like the
   two libraries mix. For example, you can do this:

   .. jupyter-execute::

      a, b, c, d = sm.symbols('a, b, c, d')

      mat = np.array([[a, b], [c, d]])
      mat

   But that will almost certainly cause you problems as you move forward. The
   process should always be:

   .. jupyter-execute::

      sym_mat = sm.Matrix([[a, b], [c, d]])
      eval_sym_mat = sm.lambdify((a, b, c, d), sym_mat)
      num_mat = eval_sym_mat(1.0, 2.0, 3.0, 4.0)
      num_mat

   Also, be careful because NumPy and SymPy have many functions that are named
   the same and you don't want to mix them up:

   .. jupyter-execute::

      np.cos(5) + sm.cos(5)

Numerical Evaluation
====================

Returning to the example of the two rods and the sliding mass from the previous
chapter, we regenerate the symbolic equations of motion and stop when we have
:math:`\bar{q}`, :math:`\bar{u}`, :math:`\mathbf{M}_k`, :math:`\bar{g}_k`,
:math:`\mathbf{M}_d`, and :math:`\bar{g}_d`. The following dropdown has the
SymPy code to generate these symbolic vectors and matrices.

.. admonition:: Symbolic Setup Code
   :class: dropdown

   .. jupyter-execute::

      m, g, kt, kl, l = sm.symbols('m, g, k_t, k_l, l')
      q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')
      u1, u2, u3 = me.dynamicsymbols('u1, u2, u3')

      N = me.ReferenceFrame('N')
      A = me.ReferenceFrame('A')
      B = me.ReferenceFrame('B')

      A.orient_axis(N, q1, N.z)
      B.orient_axis(A, q2, A.x)

      A.set_ang_vel(N, u1*N.z)
      B.set_ang_vel(A, u2*A.x)

      O = me.Point('O')
      Ao = me.Point('A_O')
      Bo = me.Point('B_O')
      Q = me.Point('Q')

      Ao.set_pos(O, l/2*A.x)
      Bo.set_pos(O, l*A.x)
      Q.set_pos(Bo, q3*B.y)

      O.set_vel(N, 0)
      Ao.v2pt_theory(O, N, A)
      Bo.v2pt_theory(O, N, A)
      Q.set_vel(B, u3*B.y)
      Q.v1pt_theory(Bo, N, B)

      t = me.dynamicsymbols._t

      qdot_repl = {q1.diff(t): u1,
                   q2.diff(t): u2,
                   q3.diff(t): u3}

      Q.set_acc(N, Q.acc(N).xreplace(qdot_repl))

      R_Ao = m*g*N.x
      R_Bo = m*g*N.x + kl*q3*B.y
      R_Q = m/4*g*N.x - kl*q3*B.y
      T_A = -kt*q1*N.z + kt*q2*A.x
      T_B = -kt*q2*A.x

      I = m*l**2/12
      I_A_Ao = I*me.outer(A.y, A.y) + I*me.outer(A.z, A.z)
      I_B_Bo = I*me.outer(B.x, B.x) + I*me.outer(B.z, B.z)

      points = [Ao, Bo, Q]
      forces = [R_Ao, R_Bo, R_Q]
      masses = [m, m, m/4]

      frames = [A, B]
      torques = [T_A, T_B]
      inertias = [I_A_Ao, I_B_Bo]

      Fr = []
      Frs = []

      for ur in [u1, u2, u3]:

         Fri = 0
         Frsi = 0

         for Pi, Ri, mi in zip(points, forces, masses):
            vr = Pi.vel(N).diff(ur, N)
            Fri += vr.dot(Ri)
            Rs = -mi*Pi.acc(N)
            Frsi += vr.dot(Rs)

         for Bi, Ti, Ii in zip(frames, torques, inertias):
            wr = Bi.ang_vel_in(N).diff(ur, N)
            Fri += wr.dot(Ti)
            Ts = -(Bi.ang_acc_in(N).dot(Ii) +
                   me.cross(Bi.ang_vel_in(N), Ii).dot(Bi.ang_vel_in(N)))
            Frsi += wr.dot(Ts)

         Fr.append(Fri)
         Frs.append(Frsi)

      Fr = sm.Matrix(Fr)
      Frs = sm.Matrix(Frs)

      q = sm.Matrix([q1, q2, q3])
      u = sm.Matrix([u1, u2, u3])

      qd = q.diff(t)
      ud = u.diff(t)

      ud_zerod = {udr: 0 for udr in ud}

      Mk = -sm.eye(3)
      gk = u

      Md = Frs.jacobian(ud)
      gd = Frs.xreplace(ud_zerod) + Fr

.. jupyter-execute::

   q, u, qd, ud

.. jupyter-execute::

   Mk, gk

.. jupyter-execute::

   Md, gd

Additionally, we will define a column vector :math:`\bar{p}` that contains all
of the constant parameters in the equations of motion. We should know these
from our problem definition but they can also be found with:

.. jupyter-execute::

   Mk.free_symbols | gk.free_symbols | Md.free_symbols | gd.free_symbols

The ``|`` operator does the union of Python sets, which is the date type that
``free_symbols`` returns. :math:`t` is not a constant parameter, but the rest
are. We can then define the symbolic :math:`p` as:

.. jupyter-execute::

   p = sm.Matrix([g, kl, kt, l, m])
   p

Now we will create a function to evaluate :math:`\mathbf{M}_k`,
:math:`\bar{g}_k`, :math:`\mathbf{M}_d`, and :math:`\bar{g}_d`. given
:math:`\bar{q}`, :math:`\bar{u}` and :math:`\bar{p}`.

.. jupyter-execute::

   eval_eom = sm.lambdify((q, u, p), [Mk, gk, Md, gd])

To test out the function ``eval_eom()`` we need some NumPy 1D arrays for
:math:`\bar{q}`, :math:`\bar{u}` and :math:`\bar{p}`.

.. warning:: Make sure to use consistent units when you introduce numbers! I
   recommend always using
   :math:`\textrm{force}=\textrm{mass}\cdot\textrm{acceleration}\rightarrow
   N=kg \cdot m s^{-2}` and :math:`\textrm{torque}=\textrm{inertia} \times
   \textrm{angular acceleration}\rightarrow N \cdot m = kg \cdot m^2 \cdot rad
   \cdot s^{-2}`.

The :external:py:func:`~numpy.deg2rad` and :external:py:func:`~numpy.rad2deg`
are helpful for angle conversions. All SymPy and NumPy trigonometric functions
operate on radians, so you'll have to convert if you prefer thinking in
degrees. My recommendation is to only use degrees when displaying the outputs,
so keep any calls to these two functions at the input and output of your whole
computation pipeline.

.. jupyter-execute::

   q_vals = np.array([
       np.deg2rad(25.0),  # q1, rad
       np.deg2rad(5.0),  # q2, rad
       0.1,  # q3, m
   ])
   q_vals, type(q_vals), q_vals.shape

.. jupyter-execute::

   u_vals = np.array([
       0.1,  # u1, rad/s
       2.2,  # u2, rad/s
       0.3,  # u3, m/s
   ])
   u_vals, type(u_vals), u_vals.shape

.. jupyter-execute::

   p_vals = np.array([
       9.81,  # g, m/s**2
       2.0,  # kl, N/m
       0.01,  # kt, Nm/rad
       0.6,  # l, m
       1.0,  # m, kg
   ])
   p_vals, type(p_vals), p_vals.shape

.. jupyter-execute::

   Mk_vals, gk_vals, Md_vals, gd_vals = eval_eom(q_vals, u_vals, p_vals)
   Mk_vals, gk_vals, Md_vals, gd_vals

Now that :external:py:func:`~numpy.linalg.solve` can be used to solve the
system of linear equations (:math:`\mathbf{A}\bar{x}=\bar{b}` type systems).

.. note:: Note the use of :external:py:func:`~numpy.squeeze`. This forces
   ``gk_vals`` and ``gd_vals`` to be a 1D array with shape(3,) instead of a 2D
   array of shape(3, 1). This then causes ``qd_vals`` and ``ud_vals`` to be 1D
   arrays.

   .. jupyter-execute::

      np.linalg.solve(-Mk_vals, gk_vals)

.. jupyter-execute::

   qd_vals = np.linalg.solve(-Mk_vals, np.squeeze(gk_vals))
   qd_vals

.. jupyter-execute::

   ud_vals = np.linalg.solve(-Md_vals, np.squeeze(gd_vals))
   ud_vals

Simulate
========

To simulate the system forward in time, we solve the `initial value problem`_
of the ordinary differential equations.

.. _initial value problem: https://en.wikipedia.org/wiki/Initial_value_problem

A simple way to do so, is to use `Euler's Method`_.

.. _Euler's Method: https://en.wikipedia.org/wiki/Euler_method

.. math::

   \bar{x}_{i + 1} = \bar{x}_i + \Delta t \bar{f}_m(t_i, \bar{x}_i, \bar{p})

.. jupyter-execute::

   def euler_integrate(rhs_func, tspan, initial_cond, p_vals, delt=0.01):
       """Returns state trajectory and corresponding values of time resulting
       from integrating the ordinary differential equations with Euler's
       Method.

       delt = 0.01  # seconds/sample

       """
       num_samples = int((tspan[1] - tspan[0])/delt)
       ts = np.linspace(tspan[0], tspan[1], num=num_samples + 1)

       x = np.empty((len(ts), len(initial_cond)))

       # Set the initial conditions to the first element.
       x[0, :] = initial_cond

       # Use a for loop to sequentially calculate each new x.
       for i, ti in enumerate(ts[:-1]):
           x[i + 1, :] = x[i, :] + delt*rhs_func(ti, x[i, :], p_vals)

       return ts, x

Now we need a Python function that represents :math:`\bar{f}_m(t_i, \bar{x}_i,
\bar{p})`. This function evaluates the right hand side of the explicity
ordinary differential equations and calculated the time derivatives of the
state.

.. jupyter-execute::

   def eval_rhs(t, x, p):
       """Return the right hand side of the explicit ordinary differential
       equations which evaluates the time derivative of the state ``x`` at time
       ``t``.

       Parameters
       ==========
       t : float
          Time in seconds.
       x : array_like, shape(6,)
          State at time t: [q1, q2, q3, u1, u2, u3]
       p : array_like, shape(5,)
          Constant parameters: [g, kl, kt, l, m]

       Returns
       =======
       xd : ndarray, shape(6,)
           Derivative of the state with respect to time.

       """

       # unpack the q and u vectors from x
       q = x[:3]
       u = x[3:]

       # evaluate the equations of motion matrices with the values of q, u, p
       Mk, gk, Md, gd = eval_eom(q, u, p)

       # solve for q' and u'
       qd = np.linalg.solve(-Mk, gk.squeeze())
       ud = np.linalg.solve(-Md, gd.squeeze())

       # pack q' and u' into a new state time derivative vector x'
       xd = np.empty_like(x)
       xd[:3] = qd
       xd[3:] = ud

       return xd

With the function evaluated and numerical values already defined above we can
check to see if it works. First combine :math:`\bar{q}` and :math:`\bar{u}`
into a single column vector ``x_0`` and pick an arbitrary of time.

.. jupyter-execute::

   x0 = np.empty(6)
   x0[:3] = q_vals
   x0[3:] = u_vals

   t0 = 0.1

Now execute the function:

.. jupyter-execute::

   eval_rhs(t0, x0, p_vals)

It seems to work, giving a result for the time derivative of the state vector.
Now we can try out the the ``euler_integrate`` function to integration from
``t0`` to ``tf``:

.. jupyter-execute::

   tf = 2.0

   ts, xs = euler_integrate(eval_rhs, (t0, tf), x0, p_vals)

Our ``euler_integrate()`` function returns the state trajectory and the
corresponding time. They look like:

.. jupyter-execute::

   ts

.. jupyter-execute::

   type(ts), ts.shape

.. jupyter-execute::

   xs

.. jupyter-execute::

   type(xs), xs.shape

Plotting Simulation Trajectories
================================

Matplotlib_ is the most widely used library for making plots. Browse `their
example gallery`_ to get an idea of the library's capabilities. We will import
matplotlib like so:

.. jupyter-execute::

   import matplotlib.pyplot as plt

.. _Matplotlib: https://matplotlib.org
.. _their example gallery: https://matplotlib.org/stable/gallery/index.html

The :external:py:func:`~matplotlib.pyplot.plot` function offers the simplest
way to plot a chart of :math:`x` values versus :math:`y` values. I designed the
output of ``euler_integrate()`` to work well with this plotting function. To
make a basic plot use:

.. jupyter-execute::

   plt.plot(ts, xs);

.. note:: The closing semicolon at the end of the statement supressesses the
   display of the returned objects from the function. See the difference here:

   .. jupyter-execute::

      plt.plot(ts, xs)

This plot shows that the state trajectory changes with respect to time, but
without any more information it is hard to interpret. The following function
uses :external:py:func:`~matplotlib.pyplot.subplots` to make a figure with
panels for the different state variables. I use
:external:py:func:`~sympy.physics.vector.printing.mlatex` to include the
symbolic symbol names in the legends.

.. jupyter-execute::

   def plot_results(ts, xs):
       """Returns the array of axes of a 4 panel plot of the state trajectory
       versus time.

       Parameters
       ==========
       ts : array_like, shape(n,)
          Values of time.
       xs : array_like, shape(n, 6)
          Values of the state trajectories corresponding to ``ts``.

       Returns
       =======
       axes : ndarray, shape(4,)
          Matplotlib axes for each panel.

       """

       fig, axes = plt.subplots(4, 1, sharex=True)

       fig.set_size_inches((10.0, 6.0))

       axes[0].plot(ts, np.rad2deg(xs[:, :2]))
       axes[1].plot(ts, xs[:, 2])
       axes[2].plot(ts, np.rad2deg(xs[:, 3:5]))
       axes[3].plot(ts, xs[:, 5])

       axes[0].legend([me.mlatex(q[0], mode='inline'),
                       me.mlatex(q[1], mode='inline')])
       axes[1].legend([me.mlatex(q[2], mode='inline')])
       axes[2].legend([me.mlatex(u[0], mode='inline'),
                       me.mlatex(u[1], mode='inline')])
       axes[3].legend([me.mlatex(u[2], mode='inline')])

       axes[0].set_ylabel('Angle [deg]')
       axes[1].set_ylabel('Distance [m]')
       axes[2].set_ylabel('Angular Rate [deg/s]')
       axes[3].set_ylabel('Speed [m/s]')

       axes[3].set_xlabel('Time [s]')

       fig.tight_layout()

       return axes

Our function now gives an interpretable view of the results:

.. jupyter-execute::

   plot_results(ts, xs);

.. todo:: Describe the results.

Integrating with SciPy
======================

Our ``euler_integrate()`` function seems to do the trick, but it all numerical
integrators suffer from numerical errors. Careful attention to `truncation
error`_ is needed for to keep the error trajectories within some acceptable
tolerance for your purposes. Euler's Method has poor error properties and there
is a large number of other numerical integration methods that provide better
results, at the cost of more complexity in their calculations.

.. _truncation error: https://en.wikipedia.org/wiki/Truncation_error_(numerical_integration)

SciPy is built on top of NumPy and provides a large assortment of battletested
numerical methods, including numerical methods for integration. We are solving
the initial problem of oridinary differential equations and SciPy includes the
function :external:py:func:`~scipy.integrate.solve_ivp` as an alternative to
our ``euler_integrate()`` function. ``solve_ivp()`` provides access to a
several different integration methods that are sutiable for different problems.
The default method is a `Runga-Kutta method` that works well for many types of
problems.

.. _Runga-Kutta method: https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods

We will only be using this function from SciPy so we can import it directly
with:

.. jupyter-execute::

   from scipy.integrate import solve_ivp

We can use ``solve_ivp()`` in much the same way as our ``euler_integrate()``
function. The difference is that ``solve_ivp()`` takes a function that
evaluates the right hand side of the ordinary differential equations that is of
the form ``f(t, x)``. Our parameter vector ``p`` must be passed to the
``args=`` optional keyword argument. If we only have one extra argument, as we
do ``f(t, x, p)``, then we must make a single element tuple ``(p_vals,)``.
Other than that, the inputs are the same as ``euler_integrate()``.
``solve_ivp()`` returns a solution object that contains quite a bit of
information (other than the trajectories). See the documentation for
``solve_ivp()`` for all the details.

Here is how we use the integrator with our previously defined system:

.. jupyter-execute::

   result = solve_ivp(eval_rhs, (t0, tf), x0, args=(p_vals,))

The time values are in the ``result.t`` attribute:

.. jupyter-execute::

   result.t

and the state trajectory is in the ``result.y`` attribute:

.. jupyter-execute::

   result.y

Note the shape of the trajectory array:

.. jupyter-execute::

   np.shape(result.y)

It is the transpose of our ``xs`` above. Knowing that we can use our
``plot_results()`` function to view the results. I use
:external:py:func:`~numpy.transpose` to transpose the array before passing it
into the plot function.

.. jupyter-execute::

   plot_results(result.t, np.transpose(result.y));

The default result is very coarse in time. This is because the underlying
integration algorthim adaptively selects the necessary time steps to stay
within the desired maximum truncation error. If you want to specify which time
values you'd like the result presented at you can do so by interpolating the
results by providing the time values with the keyword argumetn ``t_eval=``.

.. jupyter-execute::

   result = solve_ivp(eval_rhs, (t0, tf), x0, args=(p_vals,), t_eval=ts)

.. jupyter-execute::

   plot_results(result.t, np.transpose(result.y));

Now let's compare the results from ``euler_inegrate()`` with ``solve_ivp()``,
the later of which uses a Runga-Kutta method that has lower truncation error.
We'll plot only :math:`q_1`.

.. jupyter-execute::

   fig, ax = plt.subplots()
   fig.set_size_inches((10.0, 6.0))

   ax.plot(ts, np.rad2deg(xs[:, 0]), 'k',
           result.t, np.rad2deg(result.y.T[:, 0]), 'b');
   ax.legend(['euler_integrate', 'solve_ivp'])
   ax.set_xlabel('Time [s]')
   ax.set_ylabel('Angle [deg]')

You can clearly see that the Euler Method deviates from the more accurate
Runga-Kutta method. You'll need to learn more about truncation error and the
various integration methods to ensure you are getting the results you desire,
but that is all I'll go over for the purposes of this chapter.

Now set ``xs`` equal to the ``solve_ivp()`` result for use in the next
sectionn:

.. jupyter-execute::

   xs = result.y.T

Animation with Matplotlib
=========================

.. todo:: Sample time for 30 fps

Matplotlib provides tools to make animations by iterating over data and
updating the plot. I'll create a very simple set of plots that give 4 views of
points on the two bodies moving in space.

First create a function that calculates the :math:`xyz` coordinates relative to
point :math:`O`.

.. jupyter-execute::

   M = me.ReferenceFrame('M')
   M.orient_axis(N, sm.pi, N.y)

   By1 = me.Point('By1')
   By2 = me.Point('By2')
   By1.set_pos(Bo, l/2*B.y)
   By2.set_pos(Bo, -l/2*B.y)

   coordinates = O.pos_from(O).to_matrix(M)
   for point in [Bo, Q, By1, By2]:
      coordinates = coordinates.row_join(point.pos_from(O).to_matrix(M))

   eval_point_coords = sm.lambdify((q, p), coordinates)
   eval_point_coords(q_vals, p_vals)

Now create the desired figure with the initial conditions shown:

.. jupyter-execute::

   fig = plt.figure()
   fig.set_size_inches((8.0, 8.0))

   axes = []
   axes.append(fig.add_subplot(2, 2, 1))
   axes.append(fig.add_subplot(2, 2, 2, projection='3d'))
   axes.append(fig.add_subplot(2, 2, 3))
   axes.append(fig.add_subplot(2, 2, 4))

   x, y, z = eval_point_coords(q_vals, p_vals)

   line_prop = {
       'color': 'black',
       'marker': 'o',
       'markerfacecolor': 'blue',
       'markersize': 10,
   }

   # top
   top_lines, = axes[0].plot(y, z, **line_prop)
   axes[0].set_xlim((-0.5, 0.5))
   axes[0].set_ylim((-0.5, 0.5))
   axes[0].set_title('Top View')
   axes[0].set_aspect('equal')

   # 3d
   lines_3d, = axes[1].plot(y, z, x, **line_prop)
   axes[1].set_xlim((-0.5, 0.5))
   axes[1].set_ylim((-0.5, 0.5))
   axes[1].set_zlim((-0.8, 0.2))

   # front
   front_lines, = axes[2].plot(y, x, **line_prop)
   axes[2].set_xlim((-0.5, 0.5))
   axes[2].set_ylim((-0.8, 0.2))
   axes[2].set_title('Front View')
   axes[2].set_aspect('equal')

   # right
   right_lines, = axes[3].plot(z, x, **line_prop)
   axes[3].set_xlim((-0.5, 0.5))
   axes[3].set_ylim((-0.8, 0.2))
   axes[3].set_title('Right View')
   axes[3].set_aspect('equal')

   fig.tight_layout()

Create the animation update function.

.. jupyter-execute::

   import matplotlib.animation as animation

   def animate(i):
      x, y, z = eval_point_coords(xs[i, :3], p_vals)
      top_lines.set_data(y, z)
      lines_3d.set_data_3d(y, z, x)
      front_lines.set_data(y, x)
      right_lines.set_data(z, x)

   ani = animation.FuncAnimation(fig, animate, len(ts))

Display the resuls:

.. jupyter-execute::

   from IPython.display import HTML
   HTML(ani.to_jshtml(fps=30))
