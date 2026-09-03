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


m, Ixx, Iyy, Izz = sm.symbols('m, I_{xx}, I_{yy}, I_{zz}')
Ixy, Iyz, Ixz = sm.symbols('I_{xy}, I_{yz}, I_{xz}')
Fx, Fy, Fz, Mx, My, Mz = me.dynamicsymbols('F_x, F_y, F_z, M_x, M_y, M_z')
u1, u2, u3, u4, u5, u6 = me.dynamicsymbols('u1, u2, u3, u4, u5, u6')

A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')

Bo = me.Point('Bo')


# In[4]:


A_w_B = u4*B.x + u5*B.y + u6*B.z
B.set_ang_vel(A, A_w_B)

A_v_Bo = u1*B.x + u2*B.y + u3*B.z
Bo.set_vel(A, A_v_Bo)


# In[5]:


v_Bo_1 = A_v_Bo.diff(u1, A, var_in_dcm=False)
v_Bo_2 = A_v_Bo.diff(u2, A, var_in_dcm=False)
v_Bo_3 = A_v_Bo.diff(u3, A, var_in_dcm=False)
v_Bo_4 = A_v_Bo.diff(u4, A, var_in_dcm=False)
v_Bo_5 = A_v_Bo.diff(u5, A, var_in_dcm=False)
v_Bo_6 = A_v_Bo.diff(u6, A, var_in_dcm=False)

v_Bo_1, v_Bo_2, v_Bo_3, v_Bo_4, v_Bo_5, v_Bo_6


# In[6]:


w_B_1 = A_w_B.diff(u1, A, var_in_dcm=False)
w_B_2 = A_w_B.diff(u2, A, var_in_dcm=False)
w_B_3 = A_w_B.diff(u3, A, var_in_dcm=False)
w_B_4 = A_w_B.diff(u4, A, var_in_dcm=False)
w_B_5 = A_w_B.diff(u5, A, var_in_dcm=False)
w_B_6 = A_w_B.diff(u6, A, var_in_dcm=False)

w_B_1, w_B_2, w_B_3, w_B_4, w_B_5, w_B_6


# In[7]:


par_vels = me.partial_velocity([A_v_Bo, A_w_B], [u1, u2, u3, u4, u5, u6], A)

par_vels


# In[8]:


T = Mx*B.x + My*B.y + Mz*B.z
R = Fx*B.x + Fy*B.y + Fz*B.z

F1 = v_Bo_1.dot(R) + w_B_1.dot(T)
F2 = v_Bo_2.dot(R) + w_B_2.dot(T)
F3 = v_Bo_3.dot(R) + w_B_3.dot(T)
F4 = v_Bo_4.dot(R) + w_B_4.dot(T)
F5 = v_Bo_5.dot(R) + w_B_5.dot(T)
F6 = v_Bo_6.dot(R) + w_B_6.dot(T)

Fr = sm.Matrix([F1, F2, F3, F4, F4, F6])
Fr


# In[9]:


I = me.inertia(B, Ixx, Iyy, Izz, Ixy, Iyz, Ixz)

Rs = -m*Bo.acc(A)
Ts = -(B.ang_acc_in(A).dot(I) + me.cross(A_w_B, I).dot(A_w_B))

F1s = v_Bo_1.dot(Rs) + w_B_1.dot(Ts)
F2s = v_Bo_2.dot(Rs) + w_B_2.dot(Ts)
F3s = v_Bo_3.dot(Rs) + w_B_3.dot(Ts)
F4s = v_Bo_4.dot(Rs) + w_B_4.dot(Ts)
F5s = v_Bo_5.dot(Rs) + w_B_5.dot(Ts)
F6s = v_Bo_6.dot(Rs) + w_B_6.dot(Ts)

Frs = sm.Matrix([F1s, F2s, F3s, F4s, F5s, F6s])
Frs


# In[10]:


Fr + Frs


# In[11]:


u = sm.Matrix([u1, u2, u3, u4, u5, u6])
t = me.dynamicsymbols._t
ud = u.diff(t)


# In[12]:


M = -Frs.jacobian(ud)
M


# In[13]:


C = -Frs.xreplace({udi: 0 for udi in ud})
C


# In[14]:


F = Fr
F


# In[15]:


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

Ao.set_pos(O, l/2*A.x)
Bo.set_pos(O, l*A.x)

O.set_vel(N, 0)
Ao.v2pt_theory(O, N, A)
Bo.v2pt_theory(O, N, A)

Ao.vel(N), Bo.vel(N), A.ang_vel_in(N), B.ang_vel_in(N)


# In[16]:


Q = me.Point('Q')
Q.set_pos(Bo, q3*B.y)
Q.set_vel(B, u3*B.y)
Q.v1pt_theory(Bo, N, B)

Q.vel(N)


# In[17]:


Ao.acc(N), Bo.acc(N), A.ang_acc_in(N), B.ang_acc_in(N)


# In[18]:


Q.acc(N)


# In[19]:


t = me.dynamicsymbols._t

qdot_repl = {q1.diff(t): u1,
             q2.diff(t): u2,
             q3.diff(t): u3}

Q.set_acc(N, Q.acc(N).xreplace(qdot_repl))
Q.acc(N)


# In[20]:


expr = m*q1.diff(t, 2) + kt*q1.diff(t) + kl*q1
expr


# In[21]:


expr1 = expr.xreplace({q1.diff(t, 2): u1.diff(t)}).xreplace({q1.diff(t): u1}).xreplace({q1: q2/q1})
expr1


# In[22]:


expr2 = expr.xreplace({q1: q2/q1}).xreplace({q1.diff(t): u1}).xreplace({q1.diff(t, 2): u1.diff(t)})
expr2


# In[23]:


print(expr1)


# In[24]:


print(expr2)


# In[25]:


expr.xreplace({q1: q2/q1, q1.diff(t): u1, q1.diff(t, 2): u1.diff(t)})


# In[26]:


expr.xreplace({q1.diff(t, 2): u1.diff(t), q1.diff(t): u1, q1: q2/q1})


# In[27]:


R_Ao = m*g*N.x
R_Bo = m*g*N.x + kl*q3*B.y
R_Q = m/4*g*N.x - kl*q3*B.y
T_A = -kt*q1*N.z + kt*q2*A.x
T_B = -kt*q2*A.x


# In[28]:


I = m*l**2/12
I_A_Ao = I*me.outer(A.y, A.y) + I*me.outer(A.z, A.z)
I_B_Bo = I*me.outer(B.x, B.x) + I*me.outer(B.z, B.z)


# In[29]:


v_Ao_1 = Ao.vel(N).diff(u1, N)
v_Bo_1 = Bo.vel(N).diff(u1, N)
v_Q_1 = Q.vel(N).diff(u1, N)

v_Ao_2 = Ao.vel(N).diff(u2, N)
v_Bo_2 = Bo.vel(N).diff(u2, N)
v_Q_2 = Q.vel(N).diff(u2, N)

v_Ao_3 = Ao.vel(N).diff(u3, N)
v_Bo_3 = Bo.vel(N).diff(u3, N)
v_Q_3 = Q.vel(N).diff(u3, N)

w_A_1 = A.ang_vel_in(N).diff(u1, N)
w_B_1 = B.ang_vel_in(N).diff(u1, N)

w_A_2 = A.ang_vel_in(N).diff(u2, N)
w_B_2 = B.ang_vel_in(N).diff(u2, N)

w_A_3 = A.ang_vel_in(N).diff(u3, N)
w_B_3 = B.ang_vel_in(N).diff(u3, N)


# In[30]:


F1 = v_Ao_1.dot(R_Ao) + v_Bo_1.dot(R_Bo) + v_Q_1.dot(R_Q) + w_A_1.dot(T_A) + w_B_1.dot(T_B)
F2 = v_Ao_2.dot(R_Ao) + v_Bo_2.dot(R_Bo) + v_Q_2.dot(R_Q) + w_A_2.dot(T_A) + w_B_2.dot(T_B)
F3 = v_Ao_3.dot(R_Ao) + v_Bo_3.dot(R_Bo) + v_Q_3.dot(R_Q) + w_A_3.dot(T_A) + w_B_3.dot(T_B)


# In[31]:


Fr = sm.Matrix([F1, F2, F3])
Fr


# In[32]:


TAs = -(A.ang_acc_in(N).dot(I_A_Ao) + me.cross(A.ang_vel_in(N), I_A_Ao).dot(A.ang_vel_in(N)))
TBs = -(B.ang_acc_in(N).dot(I_B_Bo) + me.cross(B.ang_vel_in(N), I_B_Bo).dot(B.ang_vel_in(N)))

F1s = v_Ao_1.dot(-m*Ao.acc(N)) + v_Bo_1.dot(-m*Bo.acc(N)) + v_Q_1.dot(-m/4*Q.acc(N))
F1s += w_A_1.dot(TAs) + w_B_1.dot(TBs)

F2s = v_Ao_2.dot(-m*Ao.acc(N)) + v_Bo_2.dot(-m*Bo.acc(N)) + v_Q_2.dot(-m/4*Q.acc(N))
F2s += w_A_2.dot(TAs) + w_B_2.dot(TBs)

F3s = v_Ao_3.dot(-m*Ao.acc(N)) + v_Bo_3.dot(-m*Bo.acc(N)) + v_Q_3.dot(-m/4*Q.acc(N))
F3s += w_A_3.dot(TAs) + w_B_3.dot(TBs)


# In[33]:


Frs = sm.Matrix([F1s, F2s, F3s])
Frs


# In[34]:


me.find_dynamicsymbols(Fr)


# In[35]:


me.find_dynamicsymbols(Frs)


# In[36]:


u = sm.Matrix([u1, u2, u3])
ud = u.diff(t)
ud


# In[37]:


Md = Frs.jacobian(ud)
Md


# In[38]:


ud_zerod = {udr: 0 for udr in ud}
ud_zerod


# In[39]:


gd = Frs.xreplace(ud_zerod) + Fr
gd


# In[40]:


me.find_dynamicsymbols(Md)


# In[41]:


me.find_dynamicsymbols(gd)

