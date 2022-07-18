#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sympy as sm
import sympy.physics.mechanics as me
me.init_vprinting(use_latex='mathjax')


# In[2]:


m, Ixx, Iyy, Izz = sm.symbols('m, I_{xx}, I_{yy}, I_{zz}')
Ixy, Iyz, Ixz = sm.symbols('I_{xy}, I_{yz}, I_{xz}')
Fx, Fy, Fz, Mx, My, Mz = me.dynamicsymbols('F_x, F_y, F_z, M_x, M_y, M_z')
u1, u2, u3, u4, u5, u6 = me.dynamicsymbols('u1, u2, u3, u4, u5, u6')

A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')

Bo = me.Point('Bo')


# In[3]:


A_w_B = u4*B.x + u5*B.y + u6*B.z
B.set_ang_vel(A, A_w_B)

A_v_Bo = u1*B.x + u2*B.y + u3*B.z
Bo.set_vel(A, A_v_Bo)


# In[4]:


v_Bo_1 = A_v_Bo.diff(u1, A, var_in_dcm=False)
v_Bo_2 = A_v_Bo.diff(u2, A, var_in_dcm=False)
v_Bo_3 = A_v_Bo.diff(u3, A, var_in_dcm=False)
v_Bo_4 = A_v_Bo.diff(u4, A, var_in_dcm=False)
v_Bo_5 = A_v_Bo.diff(u5, A, var_in_dcm=False)
v_Bo_6 = A_v_Bo.diff(u6, A, var_in_dcm=False)

v_Bo_1, v_Bo_2, v_Bo_3, v_Bo_4, v_Bo_5, v_Bo_6


# In[5]:


w_B_1 = A_w_B.diff(u1, A, var_in_dcm=False)
w_B_2 = A_w_B.diff(u2, A, var_in_dcm=False)
w_B_3 = A_w_B.diff(u3, A, var_in_dcm=False)
w_B_4 = A_w_B.diff(u4, A, var_in_dcm=False)
w_B_5 = A_w_B.diff(u5, A, var_in_dcm=False)
w_B_6 = A_w_B.diff(u6, A, var_in_dcm=False)

w_B_1, w_B_2, w_B_3, w_B_4, w_B_5, w_B_6


# In[6]:


par_vels = me.partial_velocity([A_v_Bo, A_w_B], [u1, u2, u3, u4, u5, u6], A)

par_vels


# In[7]:


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


# In[8]:


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


# In[9]:


Fr + Frs


# In[10]:


u = sm.Matrix([u1, u2, u3, u4, u5, u6])
t = me.dynamicsymbols._t
ud = u.diff(t)


# In[11]:


M = -Frs.jacobian(ud)
M


# In[12]:


C = -Frs.xreplace({udi: 0 for udi in ud})
C


# In[13]:


F = Fr
F


# In[14]:


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


# In[15]:


Q = me.Point('Q')
Q.set_pos(Bo, q3*B.y)
Q.set_vel(B, u3*B.y)
Q.v1pt_theory(Bo, N, B)

Q.vel(N)


# In[16]:


Ao.acc(N), Bo.acc(N), A.ang_acc_in(N), B.ang_acc_in(N)


# In[17]:


Q.acc(N)


# In[18]:


t = me.dynamicsymbols._t

qdot_repl = {q1.diff(t): u1,
             q2.diff(t): u2,
             q3.diff(t): u3}

Q.set_acc(N, Q.acc(N).xreplace(qdot_repl))
Q.acc(N)


# In[19]:


expr = m*q1.diff(t, 2) + kt*q1.diff(t) + kl*q1
expr


# In[20]:


expr1 = expr.xreplace({q1.diff(t, 2): u1.diff(t)}).xreplace({q1.diff(t): u1}).xreplace({q1: q2/q1})
expr1


# In[21]:


expr2 = expr.xreplace({q1: q2/q1}).xreplace({q1.diff(t): u1}).xreplace({q1.diff(t, 2): u1.diff(t)})
expr2


# In[22]:


print(expr1)


# In[23]:


print(expr2)


# In[24]:


expr.xreplace({q1: q2/q1, q1.diff(t): u1, q1.diff(t, 2): u1.diff(t)})


# In[25]:


expr.xreplace({q1.diff(t, 2): u1.diff(t), q1.diff(t): u1, q1: q2/q1})


# In[26]:


R_Ao = m*g*N.x
R_Bo = m*g*N.x + kl*q3*B.y
R_Q = m/4*g*N.x - kl*q3*B.y
T_A = -kt*q1*N.z + kt*q2*A.x
T_B = -kt*q2*A.x


# In[27]:


I = m*l**2/12
I_A_Ao = I*me.outer(A.y, A.y) + I*me.outer(A.z, A.z)
I_B_Bo = I*me.outer(B.x, B.x) + I*me.outer(B.z, B.z)


# In[28]:


points = [Ao, Bo, Q]
forces = [R_Ao, R_Bo, R_Q]
masses = [m, m, m/4]

frames = [A, B]
torques = [T_A, T_B]
inertias = [I_A_Ao, I_B_Bo]

Fr_bar = []
Frs_bar = []

# loop over the three generalized speeds
for ur in [u1, u2, u3]:

    # initialize the rth GAF and GIF
    Fr = 0
    Frs = 0

    # for the rth generalized speed, loop though each point to find it's
    # contribution to the generalized forces
    for Pi, Ri, mi in zip(points, forces, masses):
        vr = Pi.vel(N).diff(ur, N)  # rth partial velocity
        Fr += vr.dot(Ri)  # sum in Pi's contribution to GAF
        Rs = -mi*Pi.acc(N)  # rth inertia force
        Frs += vr.dot(Rs)  # sum in Pi's contribution to GIF

    # for the rth generalized speed, loop though each reference frame to find
    # it's contribution to the generalized forces
    for Bi, Ti, Ii in zip(frames, torques, inertias):
        wr = Bi.ang_vel_in(N).diff(ur, N)  # rth partial velocity
        Fr += wr.dot(Ti)  # sum in Bi's contribution to the GIF
        Ts = -(Bi.ang_acc_in(N).dot(Ii) +  # rth inertia torque
               me.cross(Bi.ang_vel_in(N), Ii).dot(Bi.ang_vel_in(N)))
        Frs += wr.dot(Ts)  # sum in Bi's contribution to the GAF

    Fr_bar.append(Fr)
    Frs_bar.append(Frs)


# In[29]:


Fr = sm.Matrix(Fr_bar)
Fr


# In[30]:


Frs = sm.Matrix(Frs_bar)
Frs


# In[31]:


me.find_dynamicsymbols(Fr)


# In[32]:


me.find_dynamicsymbols(Frs)


# In[33]:


u = sm.Matrix([u1, u2, u3])
ud = u.diff(t)
ud


# In[34]:


Md = Frs.jacobian(ud)
Md


# In[35]:


ud_zerod = {udr: 0 for udr in ud}
ud_zerod


# In[36]:


gd = Frs.xreplace(ud_zerod) + Fr
gd


# In[37]:


me.find_dynamicsymbols(Md)


# In[38]:


me.find_dynamicsymbols(gd)

