#!/usr/bin/env python
# coding: utf-8

# In[1]:


from IPython.display import HTML
from matplotlib.animation import FuncAnimation
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np
import sympy as sm
import sympy.physics.mechanics as me

me.init_vprinting(use_latex='mathjax')


# In[2]:


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


# In[3]:


q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')
u1, u2, u3 = me.dynamicsymbols('u1, u2, u3')
la, lb, lc, ln = sm.symbols('l_a, l_b, l_c, l_n')
m, g = sm.symbols('m, g')
t = me.dynamicsymbols._t

p = sm.Matrix([la, lb, lc, ln, m, g])

q = sm.Matrix([q1])
qr = sm.Matrix([q2, q3])
qN = q.col_join(qr)

u = sm.Matrix([u1])
ur = sm.Matrix([u2, u3])
uN = u.col_join(ur)

qdN = qN.diff(t)
ud = u.diff(t)

p, q, qr, qN, u, ur, uN, qdN, ud


# In[4]:


ur_zero = {ui: 0 for ui in ur}
uN_zero = {ui: 0 for ui in uN}
qdN_zero = {qdi: 0 for qdi in qdN}
ud_zero = {udi: 0 for udi in ud}


# In[5]:


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

P2.set_pos(P1, la*A.x)
P3.set_pos(P2, lb*B.x)
P4.set_pos(P3, lc*C.x)


# In[6]:


loop = P4.pos_from(P1) - ln*N.x

fh = sm.Matrix([loop.dot(N.x), loop.dot(N.y)])
fh = sm.trigsimp(fh)
fh


# In[7]:


me.find_dynamicsymbols(fh)


# In[8]:


fk = sm.Matrix([
    q1.diff(t) - u1,
    q2.diff(t) - u2,
    q3.diff(t) - u3,
])
Mk = fk.jacobian(qdN)
gk = fk.xreplace(qdN_zero)
qdN_sol = -Mk.LUsolve(gk)
qd_repl = dict(zip(qdN, qdN_sol))
qd_repl


# In[9]:


fhd = fh.diff(t).xreplace(qd_repl)
fhd = sm.trigsimp(fhd)
fhd


# In[10]:


me.find_dynamicsymbols(fhd)


# In[11]:


Mhd = fhd.jacobian(ur)
ghd = fhd.xreplace(ur_zero)
ur_sol = sm.trigsimp(-Mhd.LUsolve(ghd))
ur_repl = dict(zip(ur, ur_sol))
ur_repl[u2]


# In[12]:


ur_repl[u3]


# In[13]:


A.set_ang_vel(N, u1*N.z)
B.set_ang_vel(A, ur_repl[u2]*A.z)
C.set_ang_vel(B, ur_repl[u3]*B.z)


# In[14]:


P1.set_vel(N, 0)
P2.v2pt_theory(P1, N, A)
P3.v2pt_theory(P2, N, B)
P4.v2pt_theory(P3, N, C)

(me.find_dynamicsymbols(P2.vel(N), reference_frame=N) |
 me.find_dynamicsymbols(P3.vel(N), reference_frame=N) |
 me.find_dynamicsymbols(P4.vel(N), reference_frame=N))


# In[15]:


gk = gk.xreplace(ur_repl)


# In[16]:


R_P2 = -m*g*N.y
R_P3 = -m*g*N.y


# In[17]:


Fr = sm.Matrix([
    P2.vel(N).diff(u1, N).dot(R_P2) + P3.vel(N).diff(u1, N).dot(R_P3)
])
Fr


# In[18]:


me.find_dynamicsymbols(Fr)


# In[19]:


me.find_dynamicsymbols(P2.acc(N), reference_frame=N)


# In[20]:


me.find_dynamicsymbols(P3.acc(N), reference_frame=N)


# In[21]:


Rs_P2 = -m*P2.acc(N)
Rs_P3 = -m*P3.acc(N).xreplace(qd_repl).xreplace(ur_repl)


# In[22]:


Frs = sm.Matrix([
    P2.vel(N).diff(u1, N).dot(Rs_P2) + P3.vel(N).diff(u1, N).dot(Rs_P3)
])
me.find_dynamicsymbols(Frs)


# In[23]:


Md = Frs.jacobian(ud)
gd = Frs.xreplace(ud_zero) + Fr


# In[24]:


me.find_dynamicsymbols(Mk), me.find_dynamicsymbols(gk)


# In[25]:


me.find_dynamicsymbols(Md), me.find_dynamicsymbols(gd)


# In[26]:


eval_k = sm.lambdify((qN, u, p), (Mk, gk))
eval_d = sm.lambdify((qN, u, p), (Md, gd))


def eval_rhs(t, x, p):
    """Return the derivative of the state at time t.

    Parameters
    ==========
    t : float
    x : array_like, shape(4,)
       x = [q1, q2, q3, u1]
    p : array_like, shape(6,)
       p = [la, lb, lc, ln, m, g]

    Returns
    =======
    xd : ndarray, shape(4,)
       xd = [q1d, q2d, q3d, u1d]

    """

    qN = x[:3]  # shape(3,)
    u = x[3:]   # shape(1,)

    Mk, gk = eval_k(qN, u, p)
    qNd = -np.linalg.solve(Mk, np.squeeze(gk))

    # Md, gd, and ud are each shape(1,1)
    Md, gd = eval_d(qN, u, p)
    ud = -np.linalg.solve(Md, gd)[0]

    return np.hstack((qNd, ud))


# In[27]:


p_vals = np.array([
    0.8,  # la [m]
    2.0,  # lb [m]
    1.0,  # lc [m]
    2.0,  # ln [m]
    1.0,  # m [kg]
    9.81,  # g [m/s^2]
])


# In[28]:


from scipy.optimize import fsolve


# In[29]:


eval_fh = sm.lambdify((qr, q1, p), fh)


# In[30]:


q1_val = np.deg2rad(10.0)
qr_guess = np.deg2rad([10.0, -150.0])


# In[31]:


q2_val, q3_val = fsolve(
    lambda qr, q1, p: np.squeeze(eval_fh(qr, q1, p)),  # squeeze to a 1d array
    qr_guess,  # initial guess for q2 and q3
    args=(q1_val, p_vals)) # known values in fh


# In[32]:


qN_vals = np.array([q1_val, q2_val, q3_val])
np.rad2deg(qN_vals)


# In[33]:


eval_fh(qN_vals[1:], qN_vals[0], p_vals)


# In[34]:


u1_val = 0.0
x0 = np.hstack((qN_vals, u1_val))
x0


# In[35]:


t0, tf, fps = 0.0, 30.0, 20


# In[36]:


eval_rhs(t0, x0, p_vals)


# In[37]:


def eval_constraints(xs, p):
    """Returns the value of the left hand side of the holonomic constraints
    at each time instance.

    Parameters
    ==========
    xs : ndarray, shape(n, 4)
        States at each of n time steps.
    p : ndarray, shape(6,)
        Constant parameters.

    Returns
    =======
    con : ndarray, shape(n, 2)
        fh evaluated at each xi in xs.

    """
    con = []
    for xi in xs:  # xs is shape(n, 4)
       con.append(eval_fh(xi[1:3], xi[0], p).squeeze())
    return np.array(con)


# In[38]:


def simulate(eval_rhs, t0, tf, fps, q1_0, u1_0, q2_0g, q3_0g, p):
    """Returns the simulation results.

    Parameters
    ==========
    eval_rhs : function
       Function that returns the derivatives of the states in the form:
       ``eval_rhs(t, x, p)``.
    t0 : float
       Initial time in seconds.
    tf : float
       Final time in seconds.
    fps : integer
       Number of "frames" per second to output.
    q1_0 : float
       Initial q1 angle in radians.
    u1_0 : float
       Initial u1 rate in radians/s.
    q2_0g : float
       Guess for the initial q2 angle in radians.
    q3_0g : float
       Guess for the initial q3 angle in radians.
    p : array_like, shape(6,)
       Constant parameters p = [la, lb, lc, ln, m, g].

    Returns
    =======
    ts : ndarray, shape(n,)
       Time values.
    xs : ndarray, shape(n, 4)
       State values at each time.
    con : ndarray, shape(n, 2)
       Constraint violations at each time in meters.

    """

    # generate the time steps
    ts = np.linspace(t0, tf, num=int(fps*(tf - t0)))

    # solve for the dependent coordinates
    q2_val, q3_val = fsolve(
        lambda qr, q1, p: np.squeeze(eval_fh(qr, q1, p)),
        [q2_0g, q3_0g],
        args=(q1_0, p))

    # establish the initial conditions
    x0 = np.array([q1_val, q2_val, q3_val, u1_0])

    # integrate the equations of motion
    sol = solve_ivp(eval_rhs, (ts[0], ts[-1]), x0, args=(p,), t_eval=ts,
                    rtol=1e-3, atol=1e-6)
    xs = np.transpose(sol.y)
    ts = sol.t

    # evaluate the constraints
    con = eval_constraints(xs, p)

    return ts, xs, con


# In[39]:


def plot_results(ts, xs, con):
    """Returns the array of axes of a 4 panel plot of the state trajectory
    versus time.

    Parameters
    ==========
    ts : array_like, shape(n,)
       Values of time.
    xs : array_like, shape(n, 4)
       Values of the state trajectories corresponding to ``ts`` in order
       [q1, q2, q3, u1].
    con : array_like, shape(n, 2)
       x and y constraint residuals of P4 at each time in ``ts``.

    Returns
    =======
    axes : ndarray, shape(3,)
       Matplotlib axes for each panel.

    """
    fig, axes = plt.subplots(3, 1, sharex=True)

    fig.set_size_inches((10.0, 6.0))

    axes[0].plot(ts, np.rad2deg(xs[:, :3]))  # q1(t), q2(t), q3(t)
    axes[1].plot(ts, np.rad2deg(xs[:, 3]))  # u1(t)
    axes[2].plot(ts, np.squeeze(con))  # fh(t)

    axes[0].legend(['$q_1$', '$q_2$', '$q_3$'])
    axes[1].legend(['$u_1$'])
    axes[2].legend([r'$\cdot\hat{n}_x$', r'$\cdot\hat{n}_y$'])

    axes[0].set_ylabel('Angle [deg]')
    axes[1].set_ylabel('Angular Rate [deg/s]')
    axes[2].set_ylabel('Distance [m]')
    axes[2].set_xlabel('Time [s]')

    fig.tight_layout()

    return axes


# In[40]:


ts, xs, con = simulate(
    eval_rhs,
    t0=t0,
    tf=tf,
    fps=fps,
    q1_0=np.deg2rad(10.0),
    u1_0=0.0,
    q2_0g=np.deg2rad(10.0),
    q3_0g=np.deg2rad(-150.0),
    p=p_vals,
)
plot_results(ts, xs, con);


# In[41]:


coordinates = P2.pos_from(P1).to_matrix(N)
for point in [P3, P4, P1, P2]:
   coordinates = coordinates.row_join(point.pos_from(P1).to_matrix(N))
eval_point_coords = sm.lambdify((qN, p), coordinates)


# In[42]:


def setup_animation_plot(ts, xs, p):
    """Returns objects needed for the animation.

    Parameters
    ==========
    ts : array_like, shape(n,)
       Values of time.
    xs : array_like, shape(n, 4)
       Values of the state trajectories corresponding to ``ts`` in order
       [q1, q2, q3, u1].
    p : array_like, shape(6,)

    """

    x, y, z = eval_point_coords(xs[0, :3], p)

    fig, ax = plt.subplots()
    fig.set_size_inches((10.0, 10.0))
    ax.set_aspect('equal')
    ax.grid()

    lines, = ax.plot(x, y, color='black',
                     marker='o', markerfacecolor='blue', markersize=10)

    title_text = ax.set_title('Time = {:1.1f} s'.format(ts[0]))
    ax.set_xlim((-1.0, 3.0))
    ax.set_ylim((-1.0, 1.0))
    ax.set_xlabel('$x$ [m]')
    ax.set_ylabel('$y$ [m]')

    return fig, ax, title_text, lines

setup_animation_plot(ts, xs, p_vals);


# In[43]:


def animate_linkage(ts, xs, p):
    """Returns an animation object.

    Parameters
    ==========
    ts : array_like, shape(n,)
    xs : array_like, shape(n, 4)
       x = [q1, q2, q3, u1]
    p : array_like, shape(6,)
       p = [la, lb, lc, ln, m, g]

    """
    # setup the initial figure and axes
    fig, ax, title_text, lines = setup_animation_plot(ts, xs, p)

    # precalculate all of the point coordinates
    coords = []
    for xi in xs:
        coords.append(eval_point_coords(xi[:3], p))
    coords = np.array(coords)

    # define the animation update function
    def update(i):
        title_text.set_text('Time = {:1.1f} s'.format(ts[i]))
        lines.set_data(coords[i, 0, :], coords[i, 1, :])

    # close figure to prevent premature display
    plt.close()

    # create and return the animation
    return FuncAnimation(fig, update, len(ts))


# In[44]:


HTML(animate_linkage(ts, xs, p_vals).to_jshtml(fps=fps))


# In[45]:


def eval_rhs_fsolve(t, x, p):
    """Return the derivative of the state at time t.

    Parameters
    ==========
    t : float
    x : array_like, shape(4,)
       x = [q1, q2, q3, u1]
    p : array_like, shape(6,)
       p = [la, lb, lc, ln, m, g]

    Returns
    =======
    xd : ndarray, shape(4,)
       xd = [q1d, q2d, q3d, u1d]

    Notes
    =====

    Includes a holonomic constraint correction.

    """
    qN = x[:3]
    u = x[3:]

    # correct the dependent coordinates
    qN[1:] = fsolve(lambda qr, q1, p: np.squeeze(eval_fh(qr, q1, p)),
                    qN[1:],  # guess with current solution for q2 and q3
                    args=(qN[0], p_vals))

    Mk, gk = eval_k(qN, u, p)
    qNd = -np.linalg.solve(Mk, np.squeeze(gk))

    Md, gd = eval_d(qN, u, p)
    ud = -np.linalg.solve(Md, gd)[0]

    return np.hstack((qNd, ud))


# In[46]:


ts_fsolve, xs_fsolve, con_fsolve = simulate(
    eval_rhs_fsolve,
    t0=t0,
    tf=tf,
    fps=fps,
    q1_0=np.deg2rad(10.0),
    u1_0=0.0,
    q2_0g=np.deg2rad(20.0),
    q3_0g=np.deg2rad(-150.0),
    p=p_vals,
)

plot_results(ts_fsolve, xs_fsolve, con_fsolve);


# In[47]:


HTML(animate_linkage(ts_fsolve, xs_fsolve, p_vals).to_jshtml(fps=fps))


# In[48]:


from scikits.odes import dae


# In[49]:


def eval_eom(t, x, xd, residual, p):
    """Returns the residual vector of the equations of motion.

    Parameters
    ==========
    t : float
       Time at evaluation.
    x : ndarray, shape(4,)
       State vector at time t: x = [q1, q2, q3, u1].
    xd : ndarray, shape(4,)
       Time derivative of the state vector at time t: xd = [q1d, q2d, q3d, u1d].
    residual : ndarray, shape(4,)
       Vector to store the residuals in: residuals = [fk, fd, fh1, fh2].
    p : ndarray, shape(6,)
       Constant parameters: p = [la, lb, lc, ln, m, g]

    """

    q1, q2, q3, u1 = x
    q1d, _, _, u1d = xd  # ignore the q2d and q3d values

    Md, gd = eval_d([q1, q2, q3], [u1], p)

    residual[0] = -q1d + u1  # fk, float
    residual[1] = Md[0]*u1d + gd[0]  # fd, float
    residual[2:] = eval_fh([q2, q3], [q1], p).squeeze()  # fh, shape(2,)


# In[50]:


Md_vals, gd_vals = eval_d(x0[:3], x0[3:], p_vals)

xd0 = np.array([
   0.0,  # q1d [rad/s]
   0.0,  # q2d [rad/s]
   0.0,  # q3d [rad/s]
   -np.linalg.solve(Md_vals, gd_vals)[0][0],  # u1d [rad/s^2]
])
xd0


# In[51]:


residual = np.empty(4)
residual


# In[52]:


eval_eom(t0, x0, xd0, residual, p_vals)
residual


# In[53]:


solver = dae('ida',
             eval_eom,
             rtol=1e-3,
             atol=1e-6,
             algebraic_vars_idx=[2, 3],
             user_data=p_vals,
             old_api=False)


# In[54]:


solution = solver.solve(ts, x0, xd0)

ts_dae = solution.values.t
xs_dae = solution.values.y
con_dae = eval_constraints(xs_dae, p_vals)


# In[55]:


plot_results(ts_dae, xs_dae, con_dae);


# In[56]:


HTML(animate_linkage(ts_dae, xs_dae, p_vals).to_jshtml(fps=fps))


# In[57]:


fig, ax = plt.subplots()
fig.set_size_inches((10.0, 6.0))

ax.plot(
    ts_dae, np.rad2deg(xs_dae[:, -1]), 'black',
    ts, np.rad2deg(xs[:, -1]), 'C0',
    ts_fsolve, np.rad2deg(xs_fsolve[:, -1]), 'C1',
)
ax.set_xlabel('Time [s]')
ax.set_ylabel('$u_1$ [deg/s]')
ax.legend(['IDA', 'solve_ivp', 'solve_ivp + fsolve']);


# In[58]:


solver = dae('ida',
             eval_eom,
             rtol=1e-10,
             atol=1e-10,
             algebraic_vars_idx=[2, 3],
             user_data=p_vals,
             old_api=False)

solution = solver.solve(ts, x0, xd0)

ts_dae = solution.values.t
xs_dae = solution.values.y
con_dae = eval_constraints(xs_dae, p_vals)

plot_results(ts_dae, xs_dae, con_dae);

