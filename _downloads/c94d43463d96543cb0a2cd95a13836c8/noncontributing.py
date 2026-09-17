#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


m1, m2, l1, l2, g = sm.symbols('m1, m2, l1, l2, g')
q1, q2, u1, u2, T1, T2 = me.dynamicsymbols('q1, q2, u1, u2, T1, T2')
t = me.dynamicsymbols._t

p = sm.Matrix([m1, m2, l1, l2, g])
q = sm.Matrix([q1, q2])
u = sm.Matrix([u1, u2])
r = sm.Matrix([T1, T2])

ud = u.diff(t)

p, q, u, r, ud


# In[4]:


N = me.ReferenceFrame('N')
A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')

A.orient_axis(N, q1, N.z)
B.orient_axis(N, q2, N.z)

A.set_ang_vel(N, u1*N.z)
B.set_ang_vel(N, u2*N.z)


# In[5]:


O = me.Point('O')
P1 = O.locatenew('P1', -l1*A.y)
P2 = P1.locatenew('P2', -l2*B.y)

O.set_vel(N, 0)
P1.v2pt_theory(O, N, A)


# In[6]:


P2.v2pt_theory(P1, N, B)


# In[7]:


P1.a2pt_theory(O, N, A)


# In[8]:


P2.a2pt_theory(P1, N, B)


# In[9]:


F_P1 = T1*A.y - T2*B.y - m1*g*N.y
F_P1.express(N)


# In[10]:


F_P2 = T2*B.y - m2*g*N.y
F_P2.express(N)


# In[11]:


zero_P1 = F_P1 - m1*P1.acc(N)
zero_P2 = F_P2 - m2*P2.acc(N)


# In[12]:


fd = sm.Matrix([
    zero_P1.dot(N.x),
    zero_P1.dot(N.y),
    zero_P2.dot(N.x),
    zero_P2.dot(N.y),
])
fd


# In[13]:


(me.find_dynamicsymbols(fd[0]), me.find_dynamicsymbols(fd[1]),
 me.find_dynamicsymbols(fd[2]), me.find_dynamicsymbols(fd[3]))


# In[14]:


ud, r


# In[15]:


udr = ud.col_join(r)
udr_zero = {v: 0 for v in udr}

Md = fd.jacobian(udr)
gd = fd.xreplace(udr_zero)

Md, udr, gd


# In[16]:


u3, u4 = me.dynamicsymbols('u3, u4')

N_v_P1a = P1.vel(N) - u3*A.y
N_v_P1a


# In[17]:


N_v_P2a = N_v_P1a + me.cross(B.ang_vel_in(N), P2.pos_from(P1)) - u4*B.y
N_v_P2a


# In[18]:


R_P1 = -m1*g*N.y
R_P2 = -m2*g*N.y


# In[19]:


F1 = P1.vel(N).diff(u1, N).dot(R_P1) + P2.vel(N).diff(u1, N).dot(R_P2)
F1


# In[20]:


F2 = P1.vel(N).diff(u2, N).dot(R_P1) + P2.vel(N).diff(u2, N).dot(R_P2)
F2


# In[21]:


R_P1_aux = R_P1 + T1*A.y - T2*B.y
R_P2_aux = R_P2 + T2*B.y


# In[22]:


F3 = N_v_P1a.diff(u3, N).dot(R_P1_aux) + N_v_P2a.diff(u3, N).dot(R_P2_aux)
F3


# In[23]:


F4 = N_v_P1a.diff(u4, N).dot(R_P1_aux) + N_v_P2a.diff(u4, N).dot(R_P2_aux)
F4


# In[24]:


Fr = sm.Matrix([F1, F2, F3, F4])
Fr


# In[25]:


Rs_P1 = -m1*P1.acc(N)
Rs_P2 = -m2*P2.acc(N)


# In[26]:


F1s = P1.vel(N).diff(u1, N).dot(Rs_P1) + P2.vel(N).diff(u1, N).dot(Rs_P2)
F1s


# In[27]:


F2s = P1.vel(N).diff(u2, N).dot(Rs_P1) + P2.vel(N).diff(u2, N).dot(Rs_P2)
F2s


# In[28]:


F3s = N_v_P1a.diff(u3, N).dot(Rs_P1) + N_v_P2a.diff(u3, N).dot(Rs_P2)
F3s


# In[29]:


F4s = N_v_P1a.diff(u4, N).dot(Rs_P1) + N_v_P2a.diff(u4, N).dot(Rs_P2)
F4s


# In[30]:


Frs = sm.Matrix([F1s, F2s, F3s, F4s])
Frs = sm.trigsimp(Frs)
Frs


# In[31]:


fa = Frs + Fr
me.find_dynamicsymbols(fa)


# In[32]:


Ma = fa.jacobian(udr)
ga = fa.xreplace(udr_zero)

Ma, udr, ga


# In[33]:


udr_sol = -Ma.LUsolve(ga)


# In[34]:


T1_sol = sm.trigsimp(udr_sol[2])
T1_sol


# In[35]:


T2_sol = sm.trigsimp(udr_sol[3])
T2_sol


# In[36]:


q0 = np.array([
    np.deg2rad(15.0),  # q1 [rad]
    np.deg2rad(25.0),  # q2 [rad]
])

u0 = np.array([
    np.deg2rad(123.0),  # u1 [rad/s]
    np.deg2rad(-41.0),  # u2 [rad/s]
])

p_vals = np.array([
    1.2,  # m1 [kg]
    5.6,  # m2 [kg]
    1.34,  # l1 [m]
    6.7,  # l2 [m]
    9.81,  # g [m/2^2]
])


# In[37]:


eval_d = sm.lambdify((q, u, p), (Md, gd))
eval_a = sm.lambdify((q, u, p), (Ma, ga))

Md_vals, gd_vals = eval_d(q0, u0, p_vals)
Ma_vals, ga_vals = eval_a(q0, u0, p_vals)


# In[38]:


-np.linalg.solve(Md_vals, np.squeeze(gd_vals))


# In[39]:


-np.linalg.solve(Ma_vals, np.squeeze(ga_vals))


# In[40]:


eval_forces = sm.lambdify((q, u, p), (T1_sol, T2_sol))
eval_forces(q0, u0, p_vals)


# In[41]:


def eval_rhs_newton(t, x, p):

    q = x[:2]
    u = x[2:]

    Md, gd = eval_d(q, u, p)
    udr = -np.linalg.solve(Md, np.squeeze(gd))

    qd = u
    ud = sol[:2]
    r = sol[2:]

    return np.hstack((qd, ud))

