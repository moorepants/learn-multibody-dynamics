#!/usr/bin/env python
# coding: utf-8

# In[1]:


from scipy.integrate import solve_ivp
import numpy as np
import sympy as sm
import sympy.physics.mechanics as me


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
p = sm.Matrix([g, kl, kt, l, m])

qd = q.diff(t)
ud = u.diff(t)

ud_zerod = {udr: 0 for udr in ud}

Mk = -sm.eye(3)
gk = u

Md = Frs.jacobian(ud)
gd = Frs.xreplace(ud_zerod) + Fr

eval_eom = sm.lambdify((q, u, p), [Mk, gk, Md, gd])

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

q_vals = np.array([
    np.deg2rad(25.0),  # q1, rad
    np.deg2rad(5.0),  # q2, rad
    0.1,  # q3, m
])

u_vals = np.array([
    0.1,  # u1, rad/s
    2.2,  # u2, rad/s
    0.3,  # u3, m/s
])

p_vals = np.array([
    9.81,  # g, m/s**2
    3.0,  # kl, N/m
    0.01,  # kt, Nm/rad
    0.6,  # l, m
    1.0,  # m, kg
])

x0 = np.empty(6)
x0[:3] = q_vals
x0[3:] = u_vals

fps = 20
t0, tf = 0.0, 10.0
ts = np.linspace(t0, tf, num=int(fps*(tf - t0)))
result = solve_ivp(eval_rhs, (t0, tf), x0, args=(p_vals,), t_eval=ts)
xs = result.y.T


# In[4]:


ts.shape, xs.shape


# In[5]:


import pythreejs as p3js


# In[6]:


cyl_geom = p3js.CylinderGeometry(radiusTop=2.0, radiusBottom=10.0, height=50.0)
cyl_geom


# In[7]:


red_material = p3js.MeshStandardMaterial(color='red')

cyl_mesh = p3js.Mesh(geometry=cyl_geom, material=red_material)

cyl_mesh


# In[8]:


cyl_geom = p3js.CylinderGeometry(radiusTop=0.1, radiusBottom=0.5, height=2.0)
cyl_material = p3js.MeshStandardMaterial(color='orange', wireframe=True)
cyl_mesh = p3js.Mesh(geometry=cyl_geom, material=cyl_material)
axes = p3js.AxesHelper()
cyl_mesh.add(axes)
cyl_mesh.position = (3.0, 3.0, 3.0)


# In[9]:


view_width = 600
view_height = 400

camera = p3js.PerspectiveCamera(position=[10.0, 10.0, 10.0],
                                aspect=view_width/view_height)
dir_light = p3js.DirectionalLight(position=[0.0, 10.0, 10.0])
ambient_light = p3js.AmbientLight()

axes = p3js.AxesHelper()
scene = p3js.Scene(children=[cyl_mesh, axes, camera, dir_light, ambient_light])
controller = p3js.OrbitControls(controlling=camera)
renderer = p3js.Renderer(camera=camera,
                         scene=scene,
                         controls=[controller],
                         width=view_width,
                         height=view_height)


# In[10]:


renderer


# In[11]:


cyl_mesh.matrix


# In[12]:


len(cyl_mesh.matrix)


# In[13]:


np.array(cyl_mesh.matrix).reshape(4, 4)


# In[14]:


np.array(cyl_mesh.matrix).reshape(4, 4).flatten()


# In[15]:


Ac = me.ReferenceFrame('Ac')
Ac.orient_axis(A, sm.pi/2, A.z)


# In[16]:


TA = sm.eye(4)
TA[:3, :3] = Ac.dcm(N)
TA[3, :3] = sm.transpose(Ao.pos_from(O).to_matrix(N))
TA


# In[17]:


TB = sm.eye(4)
TB[:3, :3] = B.dcm(N)
TB[3, :3] = sm.transpose(Bo.pos_from(O).to_matrix(N))
TB


# In[18]:


TQ = sm.eye(4)
TQ[:3, :3] = B.dcm(N)
TQ[3, :3] = sm.transpose(Q.pos_from(O).to_matrix(N))
TQ


# In[19]:


TA = TA.reshape(16, 1)
TB = TB.reshape(16, 1)
TQ = TQ.reshape(16, 1)


# In[20]:


TA


# In[21]:


eval_transform = sm.lambdify((q, p), (TA, TB, TQ))
eval_transform(q_vals, p_vals)


# In[22]:


TAs = []
TBs = []
TQs = []

for xi in xs:
    TAi, TBi, TQi = eval_transform(xi[:3], p_vals)
    TAs.append(TAi.squeeze().tolist())
    TBs.append(TBi.squeeze().tolist())
    TQs.append(TQi.squeeze().tolist())


# In[23]:


TAs[:2]


# In[24]:


rod_radius = p_vals[3]/20  # l/20
sphere_radius = p_vals[3]/16  # l/16

geom_A = p3js.CylinderGeometry(
    radiusTop=rod_radius,
    radiusBottom=rod_radius,
    height=p_vals[3],  # l
)

geom_B = p3js.CylinderGeometry(
    radiusTop=rod_radius,
    radiusBottom=rod_radius,
    height=p_vals[3],  # l
)

geom_Q = p3js.SphereGeometry(radius=sphere_radius)


# In[25]:


arrow_length = 0.2

mesh_A = p3js.Mesh(
    geometry=geom_A,
    material=p3js.MeshStandardMaterial(color='red'),
    name='mesh_A',
)
mesh_A.matrixAutoUpdate = False
mesh_A.add(p3js.AxesHelper(arrow_length))
mesh_A.matrix = TAs[0]

mesh_B = p3js.Mesh(
    geometry=geom_B,
    material=p3js.MeshStandardMaterial(color='blue'),
    name='mesh_B',
)
mesh_B.matrixAutoUpdate = False
mesh_B.add(p3js.AxesHelper(arrow_length))
mesh_B.matrix = TBs[0]

mesh_Q = p3js.Mesh(
    geometry=geom_Q,
    material=p3js.MeshStandardMaterial(color='green'),
    name='mesh_Q',
)
mesh_Q.matrixAutoUpdate = False
mesh_Q.add(p3js.AxesHelper(arrow_length))
mesh_Q.matrix = TQs[0]


# In[26]:


view_width = 600
view_height = 400

camera = p3js.PerspectiveCamera(position=[1.5, 0.6, 1],
                                up=[-1.0, 0.0, 0.0],
                                aspect=view_width/view_height)

key_light = p3js.DirectionalLight(position=[0, 10, 10])
ambient_light = p3js.AmbientLight()

axes = p3js.AxesHelper()

children = [mesh_A, mesh_B, mesh_Q, axes, camera, key_light, ambient_light]

scene = p3js.Scene(children=children)

controller = p3js.OrbitControls(controlling=camera)
renderer = p3js.Renderer(camera=camera, scene=scene, controls=[controller],
                         width=view_width, height=view_height)


# In[27]:


track_A = p3js.VectorKeyframeTrack(
    name="scene/mesh_A.matrix",
    times=ts,
    values=TAs
)

track_B = p3js.VectorKeyframeTrack(
    name="scene/mesh_B.matrix",
    times=ts,
    values=TBs
)

track_Q = p3js.VectorKeyframeTrack(
    name="scene/mesh_Q.matrix",
    times=ts,
    values=TQs
)


# In[28]:


tracks = [track_B, track_A, track_Q]
duration = ts[-1] - ts[0]
clip = p3js.AnimationClip(tracks=tracks, duration=duration)
action = p3js.AnimationAction(p3js.AnimationMixer(scene), clip, scene)


# In[29]:


renderer


# In[30]:


action

