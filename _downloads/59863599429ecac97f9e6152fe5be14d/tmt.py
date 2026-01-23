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


m, g, kt, kl, l = sm.symbols('m, g, k_t, k_l, l')
q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')
u1, u2, u3 = me.dynamicsymbols('u1, u2, u3')
t = me.dynamicsymbols._t

q = sm.Matrix([q1, q2, q3])
u = sm.Matrix([u1, u2, u3])
p = sm.Matrix([g, kl, kt, l, m])
q, u, p


# In[4]:


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

Ao.vel(N), A.ang_vel_in(N), Bo.vel(N), B.ang_vel_in(N), Q.vel(N)


# In[5]:


R_Ao = m*g*N.x
R_Bo = m*g*N.x + kl*q3*B.y
R_Q = m/4*g*N.x - kl*q3*B.y
T_A = -kt*q1*N.z + kt*q2*A.x
T_B = -kt*q2*A.x


# In[6]:


I = m*l**2/12
I_A_Ao = I*me.outer(A.y, A.y) + I*me.outer(A.z, A.z)
I_B_Bo = I*me.outer(B.x, B.x) + I*me.outer(B.z, B.z)


# In[7]:


v = sm.Matrix([
    Ao.vel(N).dot(N.x),
    Ao.vel(N).dot(N.y),
    Ao.vel(N).dot(N.z),
    A.ang_vel_in(N).dot(N.x),
    A.ang_vel_in(N).dot(N.y),
    A.ang_vel_in(N).dot(N.z),
    Bo.vel(N).dot(N.x),
    Bo.vel(N).dot(N.y),
    Bo.vel(N).dot(N.z),
    B.ang_vel_in(N).dot(N.x),
    B.ang_vel_in(N).dot(N.y),
    B.ang_vel_in(N).dot(N.z),
    Q.vel(N).dot(N.x),
    Q.vel(N).dot(N.y),
    Q.vel(N).dot(N.z),
])
v


# In[8]:


MA = sm.diag(m, m, m).col_join(sm.zeros(3)).row_join(sm.zeros(3).col_join(I_A_Ao.to_matrix(N)))
MA


# In[9]:


MB = sm.diag(m, m, m).col_join(sm.zeros(3)).row_join(sm.zeros(3).col_join(I_B_Bo.to_matrix(N)))
sm.trigsimp(MB)


# In[10]:


MQ = sm.diag(m/4, m/4, m/4)
MQ


# In[11]:


M = sm.diag(MA, MB, MQ)


# In[12]:


F = sm.Matrix([
    R_Ao.dot(N.x),
    R_Ao.dot(N.y),
    R_Ao.dot(N.z),
    T_A.dot(N.x),
    T_A.dot(N.y),
    T_A.dot(N.z),
    R_Bo.dot(N.x),
    R_Bo.dot(N.y),
    R_Bo.dot(N.z),
    T_B.dot(N.x),
    T_B.dot(N.y),
    T_B.dot(N.z),
    R_Q.dot(N.x),
    R_Q.dot(N.y),
    R_Q.dot(N.z),
])
F


# In[13]:


T = v.jacobian(u)
T


# In[14]:


qd_repl = dict(zip(q.diff(t), u))
ud_repl = {udi: 0 for udi in u.diff(t)}
gbar = (M*v).diff(t).xreplace(qd_repl).xreplace(ud_repl)
sm.trigsimp(gbar)


# In[15]:


Md = sm.trigsimp(-T.transpose()*M*T)
Md


# In[16]:


gd = sm.trigsimp(T.transpose()*(F - gbar))
gd


# In[17]:


u_vals = np.array([
    0.1,  # u1, rad/s
    2.2,  # u2, rad/s
    0.3,  # u3, m/s
])

q_vals = np.array([
    np.deg2rad(25.0),  # q1, rad
    np.deg2rad(5.0),  # q2, rad
    0.1,  # q3, m
])

p_vals = np.array([
    9.81,  # g, m/s**2
    2.0,  # kl, N/m
    0.01,  # kt, Nm/rad
    0.6,  # l, m
    1.0,  # m, kg
])


# In[18]:


eval_d = sm.lambdify((u, q, p), (Md, gd))

Md_vals, gd_vals = eval_d(u_vals, q_vals, p_vals)
Md_vals, gd_vals


# In[19]:


eval_d(u_vals, q_vals, p_vals)
ud_vals = -np.linalg.solve(Md_vals, np.squeeze(gd_vals))
ud_vals

