============================
Simulation and Visualization
============================

.. note::

   You can download this example as a Python script:
   :jupyter-download-script:`simulation` or Jupyter Notebook:
   :jupyter-download-notebook:`simulation`.

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

- evaluate equations of motion numerically
- numerically integrate the ordinary differential equations of a multibody
  system
- plot the system's state trajectories versus time
- compare integration methods to observe integration error
- create a simple animation of the motion of the multibody system

Numerical Integration
=====================

As mentioned at the end of the prior chapter, we will need to numerically
integrate the equations of motion. If they are in explicit form, this integral
describes how we can arrive at trajectories in time for the state variables by
integrating with respect to time from an initial time :math:`t_0` to a final
time :math:`t_f`. Recall that the time derivative of the state :math:`\bar{x}`
is:

.. math::
   :label: eq-state-form-explicit-2

   \dot{\bar{x}}(t) = \bar{f}_m(\bar{x}, t) = -\mathbf{M}_m^{-1}\bar{g}_m

We can then find :math:`\bar{x}` by integration with respect to time:

.. math::
   :label: eq-eom-integral-2

   \bar{x}(t) = \int^{t_f}_{t_0} \bar{f}_m(\bar{x}, t) dt

It is possible to form :math:`-\mathbf{M}_m^{-1}\bar{g}_m` symbolically and it
may be suitable or preferable for a given problem, but there are some possible
drawbacks. For example, if the degrees of freedom are quite large, the
resulting symbolic equations become exponentially more complex. Thus, it is
generally best to move from symbolics to numerics before formulating the
explicit ordinary differential equations.

Numerical Evaluation
====================

The NumPy_ library is the de facto base library for numeric computing with
Python. NumPy allows us to do `array programming`_ with Python by providing
floating point array data types and vectorized operators to enable repeat
operations across arrays of values. In Sec.
:ref:`Evaluating Symbolic Expressions` we introduced the SymPy function
:external:py:func:`~sympy.utilities.lambdify.lambdify`. ``lambdify()`` will be
our way to bridge the symbolic world of SymPy with the numeric world of NumPy.

.. _NumPy: https://numpy.org
.. _array programming: https://en.wikipedia.org/wiki/Array_programming

We will import NumPy like so, by convention:

.. jupyter-execute::

   import numpy as np

.. warning::

   Beware that mixing SymPy and NumPy data types will rarely, if at all,
   provide you with functioning code. Be careful because sometimes it may look
   like the two libraries mix. For example, you can do this:

   .. jupyter-execute::

      a, b, c, d = sm.symbols('a, b, c, d')

      mat = np.array([[a, b], [c, d]])
      mat

   which gives a NumPy array containing SymPy symbols. But this will almost
   certainly cause you problems as you move forward. The process you should
   always follow for the purposes of this text is:

   .. jupyter-execute::

      sym_mat = sm.Matrix([[a, b], [c, d]])
      eval_sym_mat = sm.lambdify((a, b, c, d), sym_mat)
      num_mat = eval_sym_mat(1.0, 2.0, 3.0, 4.0)
      num_mat

   Also, be careful because NumPy and SymPy have many functions that are named
   the same and you likley don't want to mix them up:

   .. jupyter-execute::

      np.cos(5) + sm.cos(5)

   We import NumPy as ``np`` and SymPy as ``sm`` to ensure functions with the
   same names can coexist.

Returning to the example of the two rods and the sliding mass from the previous
chapter, we regenerate the symbolic equations of motion and stop when we have
:math:`\bar{q}`, :math:`\bar{u}`, :math:`\mathbf{M}_k`, :math:`\bar{g}_k`,
:math:`\mathbf{M}_d`, and :math:`\bar{g}_d`. The following drop down has the
SymPy code to generate these symbolic vectors and matrices take from the prior
chapter.

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

      Fr_bar = []
      Frs_bar = []

      for ur in [u1, u2, u3]:

          Fr = 0
          Frs = 0

          for Pi, Ri, mi in zip(points, forces, masses):
              vr = Pi.vel(N).diff(ur, N)
              Fr += vr.dot(Ri)
              Rs = -mi*Pi.acc(N)
              Frs += vr.dot(Rs)

          for Bi, Ti, Ii in zip(frames, torques, inertias):
              wr = Bi.ang_vel_in(N).diff(ur, N)
              Fr += wr.dot(Ti)
              Ts = -(Bi.ang_acc_in(N).dot(Ii) +
                     me.cross(Bi.ang_vel_in(N), Ii).dot(Bi.ang_vel_in(N)))
              Frs += wr.dot(Ts)

          Fr_bar.append(Fr)
          Frs_bar.append(Frs)

      Fr = sm.Matrix(Fr_bar)
      Frs = sm.Matrix(Frs_bar)

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
from our problem definition but they can also be found using
``free_symbols()``:

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
   :math:`\textrm{force}=\textrm{mass}\times\textrm{acceleration}\rightarrow
   N=kg \ m \cdot s^{-2}` and :math:`\textrm{torque}=\textrm{inertia} \times
   \textrm{angular acceleration}\rightarrow N \ m = kg \ m^2 \cdot rad
   \ s^{-2}`.

The :external:py:func:`~numpy.deg2rad` and :external:py:func:`~numpy.rad2deg`
are helpful for angle conversions. All SymPy and NumPy trigonometric functions
operate on radians, so you'll have to convert if you prefer thinking in
degrees. My recommendation is to only use degrees when displaying the outputs,
so keep any calls to these two functions at the input and output of your whole
computation pipeline.

Here I introduce ``q_vals``, ``u_vals``, and ``p_vals``, each a 1D NumPy array.
Make sure to use a different variable name than your symbols so you can
distinguish the symbolic and numeric matrices and arrays.

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

Now we can call ``eval_eom`` with the numeric inputs to get the numerical
values of all of the equation of motion matrices and vectors:

.. jupyter-execute::

   Mk_vals, gk_vals, Md_vals, gd_vals = eval_eom(q_vals, u_vals, p_vals)
   Mk_vals, gk_vals, Md_vals, gd_vals

Now we can solve for the state derivatives, :math:`\dot{\bar{q}}` and
:math:`\dot{\bar{u}}`, numerically using NumPy's
:external:py:func:`~numpy.linalg.solve` function (not the same as SymPy's
``solve()``!) for linear systems of equations
(:math:`\mathbf{A}\bar{x}=\bar{b}` type systems).

We first numerically solve the kinematical differential equations for
:math:`\dot{\bar{q}}`:

.. jupyter-execute::

   qd_vals = np.linalg.solve(-Mk_vals, np.squeeze(gk_vals))
   qd_vals

In this case, :math:`\dot{\bar{q}}=\bar{u}` but for nontrivial generalized
speed definitions that will not be so. This next linear system solve gives the
accelerations :math:`\dot{\bar{u}}`:

.. jupyter-execute::

   ud_vals = np.linalg.solve(-Md_vals, np.squeeze(gd_vals))
   ud_vals

.. note:: Note the use of :external:py:func:`~numpy.squeeze`. This forces
   ``gk_vals`` and ``gd_vals`` to be a 1D array with shape(3,) instead of a 2D
   array of shape(3, 1). This then causes ``qd_vals`` and ``ud_vals`` to be 1D
   arrays instead of 2D.

   .. jupyter-execute::

      np.linalg.solve(-Mk_vals, gk_vals)

Simulation
==========

To simulate the system forward in time, we solve the `initial value problem`_
of the ordinary differential equations by numerically integrating
:math:`\bar{f}_m(t, \bar{x}, \bar{p})`. A simple way to do so, is to use
`Euler's Method`_:

.. math::
   :label: eq-eulers-method

   \bar{x}_{i + 1} = \bar{x}_i + \Delta t \bar{f}_m(t_i, \bar{x}_i, \bar{p})

Starting with :math:`t_i=t_0` and some initial values of the states
:math:`\bar{x}_i=\bar{x}_0`, the state at :math:`\Delta t` in the future is
computed. We repeat this until :math:`t_i=t_f` to find the trajectories of
:math:`\bar{x}` with respect to time.

.. _initial value problem: https://en.wikipedia.org/wiki/Initial_value_problem
.. _Euler's Method: https://en.wikipedia.org/wiki/Euler_method

The following function implements Euler's Method:

.. jupyter-execute::

   def euler_integrate(rhs_func, tspan, x0_vals, p_vals, delt=0.03):
       """Returns state trajectory and corresponding values of time resulting
       from integrating the ordinary differential equations with Euler's
       Method.

       Parameters
       ==========
       rhs_func : function
          Python function that evaluates the derivative of the state and takes
          this form ``dxdt = f(t, x, p)``.
       tspan : 2-tuple of floats
          The initial time and final time values: (t0, tf).
       x0_vals : array_like, shape(2*n,)
          Values of the state x at t0.
       p_vals : array_like, shape(o,)
          Values of constant parameters.
       delt : float
          Integration time step in seconds/step.

       Returns
       =======
       ts : ndarray(m, )
          Monotonically increasing values of time.
       xs : ndarray(m, 2*n)
          State values at each time in ts.

       """
       # generate monotonically increasing values of time.
       duration = tspan[1] - tspan[0]
       num_samples = round(duration/delt) + 1
       ts = np.arange(tspan[0], tspan[0] + delt*num_samples, delt)

       # create an empty array to hold the state values.
       x = np.empty((len(ts), len(x0_vals)))

       # set the initial conditions to the first element.
       x[0, :] = x0_vals

       # use a for loop to sequentially calculate each new x.
       for i, ti in enumerate(ts[:-1]):
           x[i + 1, :] = x[i, :] + delt*rhs_func(ti, x[i, :], p_vals)

       return ts, x

I used :external:py:func:`~numpy.linspace` to generate equally spaced values
between :math:`t_0` and :math:`t_f`. Now we need a Python function that
represents :math:`\bar{f}_m(t_i, \bar{x}_i, \bar{p})`. This function evaluates
the right hand side of the explicitly ordinary differential equations which
calculates the time derivatives of the state.

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
           Derivative of the state with respect to time at time ``t``.

       """

       # unpack the q and u vectors from x
       q = x[:3]
       u = x[3:]

       # evaluate the equations of motion matrices with the values of q, u, p
       Mk, gk, Md, gd = eval_eom(q, u, p)

       # solve for q' and u'
       qd = np.linalg.solve(-Mk, np.squeeze(gk))
       ud = np.linalg.solve(-Md, np.squeeze(gd))

       # pack dq/dt and du/dt into a new state time derivative vector dx/dt
       xd = np.empty_like(x)
       xd[:3] = qd
       xd[3:] = ud

       return xd

With the function evaluated and numerical values already defined above we can
check to see if it works. First combine :math:`\bar{q}` and :math:`\bar{u}`
into a single column vector of the initial conditions ``x0`` and pick an
arbitrary value for time.

.. jupyter-execute::

   x0 = np.empty(6)
   x0[:3] = q_vals
   x0[3:] = u_vals

   t0 = 0.1

Now execute the function:

.. jupyter-execute::

   eval_rhs(t0, x0, p_vals)

It seems to work, giving a result for the time derivative of the state vector,
matching the results we had above. Now we can try out the ``euler_integrate()``
function to integration from ``t0`` to ``tf``:

.. jupyter-execute::

   tf = 2.0

.. jupyter-execute::

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

Matplotlib_ is the most widely used Python library for making plots. Browse
`their example gallery`_ to get an idea of the library's capabilities. We will
use matplotlib to visualize the state trajectories and animate our system. The
convention for importing the main functionality of matplotlib is:

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

.. note:: The closing semicolon at the end of the statement suppresses the
   display of the returned objects in Jupyter. See the difference here:

   .. jupyter-execute::

      plt.plot(ts, xs)

This plot shows that the state trajectory changes with respect to time, but
without some more information it is hard to interpret. The following function
uses :external:py:func:`~matplotlib.pyplot.subplots` to make a figure with
panels for the different state variables. I use
:external:py:func:`~sympy.physics.vector.printing.vlatex` to include the
symbolic symbol names in the legends. The other matplotlib functions and
methods I use are:

- :external:py:meth:`~matplotlib.figure.Figure.set_size_inches`
- :external:py:meth:`~matplotlib.axes.Axes.plot`
- :external:py:meth:`~matplotlib.axes.Axes.legend`
- :external:py:meth:`~matplotlib.axes.Axes.set_ylabel`
- :external:py:meth:`~matplotlib.axes.Axes.set_xlabel`
- :external:py:meth:`~matplotlib.figure.Figure.tight_layout`

I also make use of array slicing notation to select which rows and columns I
want from each array. See the NumPy documentation `Indexing on ndarrays`_ for
information on how this works.

.. _Indexing on ndarrays: https://numpy.org/doc/stable/user/basics.indexing.html

.. jupyter-execute::

   def plot_results(ts, xs):
       """Returns the array of axes of a 4 panel plot of the state trajectory
       versus time.

       Parameters
       ==========
       ts : array_like, shape(m,)
          Values of time.
       xs : array_like, shape(m, 6)
          Values of the state trajectories corresponding to ``ts`` in order
          [q1, q2, q3, u1, u2, u3].

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

       axes[0].legend([me.vlatex(q[0], mode='inline'),
                       me.vlatex(q[1], mode='inline')])
       axes[1].legend([me.vlatex(q[2], mode='inline')])
       axes[2].legend([me.vlatex(u[0], mode='inline'),
                       me.vlatex(u[1], mode='inline')])
       axes[3].legend([me.vlatex(u[2], mode='inline')])

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

We now see that :math:`q_1` oscillates between :math:`\pm 40 \textrm{deg}` with
a single period. :math:`q_2` grows to around :math:`\pm 100 \textrm{deg}`, and
:math:`q_3` has half an oscillation between -0.2 and 0.2 meters. For the
initial conditions and constants we choose, this seems physically feasible.

Integration with SciPy
======================

Our ``euler_integrate()`` function seems to do the trick, but all numerical
integrators suffer from numerical errors. Careful attention to `truncation
error`_ is needed to keep the error in the resulting trajectories within some
acceptable tolerance for your problem's needs. Euler's Method has poor
truncation error unless very small time steps are chosen. But more time steps
results in longer computation time. There are a large number of other numerical
integration methods that provide better results with fewer time steps, but at
the cost of more complexity in the integration algorithm.

.. _truncation error: https://en.wikipedia.org/wiki/Truncation_error_(numerical_integration)

SciPy_ is built on top of NumPy and provides a large assortment of battle
tested numerical methods for NumPy arrays, including numerical methods for
integration. We are solving the initial value problem of ordinary differential
equations and SciPy includes the function
:external:py:func:`~scipy.integrate.solve_ivp` for this purpose.
``solve_ivp()`` provides access to a several different integration methods that
are suitable for different problems. The default method used is a `Runge-Kutta
method`_ that works well for non-stiff problems.

.. _SciPy: https://www.scipy.org
.. _Runge-Kutta method: https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods

We will only be using ``solve_ivp()`` from SciPy so we can import it directly
with:

.. jupyter-execute::

   from scipy.integrate import solve_ivp

We can use ``solve_ivp()`` in much the same way as our ``euler_integrate()``
function (in fact I designed ``euler_integrate()`` to mimic ``solve_ivp()``).
The difference is that ``solve_ivp()`` takes a function that evaluates the
right hand side of the ordinary differential equations that is of the form
``f(t, x)`` (no ``p``!). Our parameter vector ``p`` must be passed to the
``args=`` optional keyword argument in ``solve_ivp()`` to get things to work.
If we only have one extra argument, as we do ``f(t, x, p)``, then we must make
a 1-tuple ``(p_vals,)``.  Other than that, the inputs are the same as
``euler_integrate()``.  ``solve_ivp()`` returns a solution object that contains
quite a bit of information (other than the trajectories). See the documentation
for :external:py:func:`~scipy.integrate.solve_ivp` for all the details and more
examples.

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

It is the transpose of our ``xs`` computed above. Knowing that we can use our
``plot_results()`` function to view the results. I use
:external:py:func:`~numpy.transpose` to transpose the array before passing it
into the plot function.

.. jupyter-execute::

   plot_results(result.t, np.transpose(result.y));

The default result is very coarse in time (only 10 steps!). This is because the
underlying integration algorithm adaptively selects the necessary time steps to
stay within the desired maximum truncation error. The Runge-Kutta method gives
good accuracy with fewer integration steps in this case.

If you want to specify which time values you'd like the result presented at you
can do so by interpolating the results by providing the time values with the
keyword argument ``t_eval=``.

.. jupyter-execute::

   result = solve_ivp(eval_rhs, (t0, tf), x0, args=(p_vals,), t_eval=ts)

.. jupyter-execute::

   plot_results(result.t, np.transpose(result.y));

Lastly, let's compare the results from ``euler_inegrate()`` with
``solve_ivp()``, the later of which uses a Runge-Kutta method that has lower
truncation error.  We'll plot only :math:`q_1` for this comparison.

.. jupyter-execute::

   fig, ax = plt.subplots()
   fig.set_size_inches((10.0, 6.0))

   ax.plot(ts, np.rad2deg(xs[:, 0]), 'k',
           result.t, np.rad2deg(np.transpose(result.y)[:, 0]), 'b');
   ax.legend(['euler_integrate', 'solve_ivp'])
   ax.set_xlabel('Time [s]')
   ax.set_ylabel('Angle [deg]');

You can clearly see that the Euler Method deviates from the more accurate
Runge-Kutta method. You'll need to learn more about truncation error and the
various integration methods to ensure you are getting the results you desire.
For now, be aware that truncation error and `floating point arithmetic
error`_ can give you inaccurate results.

.. _floating point arithmetic error: https://en.wikipedia.org/wiki/Floating-point_arithmetic

Now set ``xs`` equal to the ``solve_ivp()`` result for use in the next section:

.. jupyter-execute::

   xs = np.transpose(result.y)

Animation with Matplotlib
=========================

Matplotlib also provides tools to make animations by iterating over data and
updating the plot. I'll create a very simple set of plots that give 4 views of
interesting points in our system.

Matplotlib's plot axes default to displaying the abscissa (:math:`x`)
horizontal and positive towards the right and the ordinate (:math:`y`) vertical
and positive upwards. The coordinate system in
:numref:`fig-eom-double-rod-pendulum` has :math:`\hat{n}_x` positive downwards
and :math:`\hat{n}_y` positive to the right. We can create a viewing reference
frame :math:`M` that matches matplotlib's axes like so:

.. jupyter-execute::

   M = me.ReferenceFrame('M')
   M.orient_axis(N, sm.pi/2, N.z)

Now :math:`\hat{m}_x` is positive to the right, :math:`\hat{m}_y` is positive
upwards, and :math:`\hat{m}_z` points out of the screen.

I'll also introduce a couple of points on each end of the rod :math:`B`, just
for visualization purposes:

.. jupyter-execute::

   Bl = me.Point('B_l')
   Br = me.Point('B_r')
   Bl.set_pos(Bo, -l/2*B.y)
   Br.set_pos(Bo, l/2*B.y)

Now, we can project the four points :math:`B_o,Q,B_l,B_r` onto the unit vectors
of :math:`M` using ``lambdify()`` to get the Cartesian coordinates of each
point relative to point :math:`O`. I use
:external:py:meth:`~sympy.matrices.common.MatrixCommon.row_join` to stack the
matrices together to build a single matrix with all points' coordinates.

.. jupyter-execute::

   coordinates = O.pos_from(O).to_matrix(M)
   for point in [Bo, Q, Bl, Br]:
       coordinates = coordinates.row_join(point.pos_from(O).to_matrix(M))

   eval_point_coords = sm.lambdify((q, p), coordinates)
   eval_point_coords(q_vals, p_vals)

The first row are the :math:`x` coordinates of each point, the second row has
the :math:`y` coordinates, and the last the :math:`z` coordinates.

Now create the desired 4 panel figure with three 2D views of the system and one
with a 3D view using the initial conditions and constant parameters shown. I
make use of :external:py:meth:`~matplotlib.figure.Figure.add_subplot` to
control if the panel is 2D or 3D.
:external:py:meth:`~matplotlib.axes.Axes.set_aspect` ensures that the abscissa
and ordinate dimensions display in a 1:1 ratio.

.. jupyter-execute::

   # initial configuration of the points
   x, y, z = eval_point_coords(q_vals, p_vals)

   # create a figure
   fig = plt.figure()
   fig.set_size_inches((10.0, 10.0))

   # setup the subplots
   ax_top = fig.add_subplot(2, 2, 1)
   ax_3d = fig.add_subplot(2, 2, 2, projection='3d')
   ax_front = fig.add_subplot(2, 2, 3)
   ax_right = fig.add_subplot(2, 2, 4)

   # common line and marker properties for each panel
   line_prop = {
       'color': 'black',
       'marker': 'o',
       'markerfacecolor': 'blue',
       'markersize': 10,
   }

   # top view
   lines_top, = ax_top.plot(x, z, **line_prop)
   ax_top.set_xlim((-0.5, 0.5))
   ax_top.set_ylim((0.5, -0.5))
   ax_top.set_title('Top View')
   ax_top.set_xlabel('x')
   ax_top.set_ylabel('z')
   ax_top.set_aspect('equal')

   # 3d view
   lines_3d, = ax_3d.plot(x, z, y, **line_prop)
   ax_3d.set_xlim((-0.5, 0.5))
   ax_3d.set_ylim((0.5, -0.5))
   ax_3d.set_zlim((-0.8, 0.2))
   ax_3d.set_xlabel('x')
   ax_3d.set_ylabel('z')
   ax_3d.set_zlabel('y')

   # front view
   lines_front, = ax_front.plot(x, y, **line_prop)
   ax_front.set_xlim((-0.5, 0.5))
   ax_front.set_ylim((-0.8, 0.2))
   ax_front.set_title('Front View')
   ax_front.set_xlabel('x')
   ax_front.set_ylabel('y')
   ax_front.set_aspect('equal')

   # right view
   lines_right, = ax_right.plot(z, y, **line_prop)
   ax_right.set_xlim((0.5, -0.5))
   ax_right.set_ylim((-0.8, 0.2))
   ax_right.set_title('Right View')
   ax_right.set_xlabel('z')
   ax_right.set_ylabel('y')
   ax_right.set_aspect('equal')

   fig.tight_layout()

Now we will use :external:py:class:`~matplotlib.animation.FuncAnimation` to
generate an animation. See the `animation examples`_ for more information on
creating animations with matplotib.

.. _animation examples: https://matplotlib.org/3.5.1/gallery/index.html#animation

First import ``FuncAnimation()``:

.. jupyter-execute::

   from matplotlib.animation import FuncAnimation

Now create a function that takes an frame index ``i``, calculates the
configuration of the points for the i\ :sup:`th` state in ``xs``, and updates
the data for the lines we have already plotted with
:external:py:meth:`~matplotlib.lines.Line2D.set_data` and
:external:py:meth:`~mpl_toolkits.mplot3d.art3d.Line3D.set_data_3d`.

.. jupyter-execute::

   def animate(i):
       x, y, z = eval_point_coords(xs[i, :3], p_vals)
       lines_top.set_data(x, z)
       lines_3d.set_data_3d(x, z, y)
       lines_front.set_data(x, y)
       lines_right.set_data(z, y)

Now provide the figure, the animation update function, and the number of frames
to ``FuncAnimation``:

.. jupyter-execute::

   ani = FuncAnimation(fig, animate, len(ts))

``FuncAnimation`` can create an interactive animation, movie files, and other
types of outputs. Here I take advantage of IPython's HTML display function and
the :external:py:meth:`~matplotlib.animation.Animation.to_jshtml` method to
create a web browser friendly visualization of the animation.

.. jupyter-execute::

   from IPython.display import HTML

   HTML(ani.to_jshtml(fps=30))

If we've setup our animation correctly and our equations of motion are correct,
we should see physically believable motion of our system. In this case, it
looks like we've successfully simulated and visualized our first multibody
system!
