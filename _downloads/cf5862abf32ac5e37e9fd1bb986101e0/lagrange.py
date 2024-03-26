#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# In[4]:


N_w_B = B.ang_vel_in(N)
r_O_P = x*N.x + y*N.y + z*N.z
N_v_C = r_O_P.dt(N)
K = N_w_B.dot(I_B.dot(N_w_B))/2 + m*N_v_C.dot(N_v_C)/2
K


# In[5]:


psid = psi.diff(t)
F_psi_s = K.diff(psid).diff(t) - K.diff(psi)


# In[6]:


K_as_matrix = sm.Matrix([K])
Fs_transposed = K_as_matrix.jacobian(qd).diff(t) - K_as_matrix.jacobian(q)
Fs = Fs_transposed.transpose()
Fs


# In[7]:


qdd_zerod = {qddr: 0 for qddr in qdd}
Md = Fs.jacobian(qdd)
gd = Fs.xreplace(qdd_zerod)
Md.simplify()
gd.simplify()
Md, gd


# In[8]:


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


# In[9]:


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


# In[10]:


left_hand_side = L.jacobian(qd).diff(t) - L.jacobian(q)
qdd_zerod = {qddr: 0 for qddr in qdd}
Md = left_hand_side.jacobian(qdd)
gd = left_hand_side.xreplace(qdd_zerod)
me.find_dynamicsymbols(Md), me.find_dynamicsymbols(gd)


# In[11]:


p = L.jacobian(qd).transpose()
sm.trigsimp(p)


# In[12]:


psi,theta, phi, x, y, z = me.dynamicsymbols('psi theta phi x y z')
N = me.ReferenceFrame('N')
B = me.ReferenceFrame('B')
B.orient_body_fixed(N, (psi, theta, phi), 'zxy')

# Mass and inertia
m, Ixx, Iyy, Izz = sm.symbols('M, I_{xx}, I_{yy}, I_{zz}')
I_B = me.inertia(B, Ixx, Iyy, Izz)


# In[13]:


omega_B = B.ang_vel_in(N)
r_com = x*N.x + y*N.y + z*N.z
v_com = r_com.dt(N)
K = omega_B.dot(I_B.dot(omega_B))/2 + m*v_com.dot(v_com)/2


# In[14]:


t = me.dynamicsymbols._t
q = sm.Matrix([psi, theta, phi, x, y, z])
qd = q.diff(t)
qdd = qd.diff(t)

L = sm.Matrix([K])
left_hand_side = L.jacobian(qd).diff(t) - L.jacobian(q)

qdd_zerod = {qddr: 0 for qddr in qdd}
Md = left_hand_side.jacobian(qdd)
gd = left_hand_side.xreplace(qdd_zerod)


# In[15]:


lambda1, lambda2, lambda3 = me.dynamicsymbols('lambda1, lambda2, lambda3')
constraint = (v_com + B.ang_vel_in(N).cross(-N.z)).to_matrix(N)
M_hn = constraint.jacobian(qd)
diff_constraint = constraint.diff(t)
sm.trigsimp(constraint)


# In[16]:


diff_constraint.jacobian(qdd) - M_hn


# In[17]:


constraint_bias = diff_constraint.xreplace({qddr : 0 for qddr in qdd})

