#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sympy as sm
import sympy.physics.mechanics as me
me.init_vprinting(use_latex='mathjax')


# In[2]:


L = sm.symbols('L')
q1, q2, u1, u2 = me.dynamicsymbols('q1, q2, u1, u2')

N = me.ReferenceFrame('N')
R = me.ReferenceFrame('R')

R.orient_axis(N, q2, N.z)


# In[3]:


N_v_A = u1*N.x
N_v_A


# In[4]:


N_w_R = u2*N.z
N_w_R


# In[5]:


r_A_B = -L*R.x
N_v_B = N_v_A + me.cross(N_w_R, r_A_B)

N_v_B.express(N)


# In[6]:


v_A_1 = N_v_A.diff(u1, N)
v_A_2 = N_v_A.diff(u2, N)

v_A_1, v_A_2


# In[7]:


v_B_1 = N_v_B.diff(u1, N)
v_B_2 = N_v_B.diff(u2, N)

v_B_1, v_B_2


# In[8]:


w_R_1 = N_w_R.diff(u1, N)
w_R_2 = N_w_R.diff(u2, N)

w_R_1, w_R_2


# In[9]:


me.partial_velocity((N_v_A, N_v_B, N_w_R), (u1, u2), N)


# In[10]:


l = sm.symbols('l')
q1, q2, u1, u2 = me.dynamicsymbols('q1, q2, u1, u2')

N = me.ReferenceFrame('N')
A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')

A.orient_axis(N, q1, N.z)
B.orient_axis(N, q2, N.z)

O = me.Point('O')
P1 = me.Point('P1')
P2 = me.Point('P2')

O.set_vel(N, 0)

P1.set_pos(O, -l*A.y)
P2.set_pos(P1, -l*B.y)

P1.v2pt_theory(O, N, A)
P2.v2pt_theory(P1, N, B)

P1.vel(N), P2.vel(N)


# In[11]:


repl = {q1.diff(): u1, q2.diff(): u2}

N_v_P1 = P1.vel(N).xreplace(repl)
N_v_P2 = P2.vel(N).xreplace(repl)

N_v_P1, N_v_P2


# In[12]:


v_P1_1 = N_v_P1.diff(u1, N)
v_P1_2 = N_v_P1.diff(u2, N)
v_P2_1 = N_v_P2.diff(u1, N)
v_P2_2 = N_v_P2.diff(u2, N)
v_P1_1, v_P1_2, v_P2_1, v_P2_2


# In[13]:


T1, T2 = me.dynamicsymbols('T1, T2')
m1, m2, g = sm.symbols('m1, m2, g')

R1 = -m1*g*N.y + T1*A.y - T2*B.y
R1


# In[14]:


R2 = -m2*g*N.y + T2*B.y
R2


# In[15]:


F1 = me.dot(v_P1_1, R1) + me.dot(v_P2_1, R2)
F1


# In[16]:


F2 = me.dot(v_P1_2, R1) + me.dot(v_P2_2, R2)
F2


# In[17]:


m, g, k, l = sm.symbols('m, g, k, l')
q1, q2, u1, u2 = me.dynamicsymbols('q1, q2, u1, u2')

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


# In[18]:


v_Ao_1 = Ao.vel(N).diff(u1, N)
v_Ao_2 = Ao.vel(N).diff(u2, N)
v_Bo_1 = Bo.vel(N).diff(u1, N)
v_Bo_2 = Bo.vel(N).diff(u2, N)

v_Ao_1, v_Ao_2, v_Bo_1, v_Bo_2


# In[19]:


w_A_1 = A.ang_vel_in(N).diff(u1, N)
w_A_2 = A.ang_vel_in(N).diff(u2, N)
w_B_1 = B.ang_vel_in(N).diff(u1, N)
w_B_2 = B.ang_vel_in(N).diff(u2, N)

w_A_1, w_A_2, w_B_1, w_B_2


# In[20]:


R_Ao = m*g*N.x
R_Bo = m*g*N.x

R_Ao, R_Bo


# In[21]:


T_A = -k*q1*N.z + k*q2*A.x
T_B = -k*q2*A.x

T_A, T_B


# In[22]:


F1_A = v_Ao_1.dot(R_Ao) + w_A_1.dot(T_A)
F1_B = v_Bo_1.dot(R_Bo) + w_B_1.dot(T_B)
F2_A = v_Ao_2.dot(R_Ao) + w_A_2.dot(T_A)
F2_B = v_Bo_2.dot(R_Bo) + w_B_2.dot(T_B)

F1_A, F1_B, F2_A, F2_B


# In[23]:


F1 = F1_A + F1_B
F2 = F2_A + F2_B

Fr = sm.Matrix([F1, F2])
Fr


# In[24]:


m, g, k, l = sm.symbols('m, g, k, l')
q1, q2, u1, u2 = me.dynamicsymbols('q1, q2, u1, u2')

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

v_Ao_1 = Ao.vel(N).diff(u1, N)
v_Ao_2 = Ao.vel(N).diff(u2, N)
v_Bo_1 = Bo.vel(N).diff(u1, N)
v_Bo_2 = Bo.vel(N).diff(u2, N)

w_A_1 = A.ang_vel_in(N).diff(u1, N)
w_A_2 = A.ang_vel_in(N).diff(u2, N)
w_B_1 = B.ang_vel_in(N).diff(u1, N)
w_B_2 = B.ang_vel_in(N).diff(u2, N)


# In[25]:


Ao.acc(N), Bo.acc(N)


# In[26]:


A.ang_acc_in(N), B.ang_acc_in(N)


# In[27]:


I = m*l**2/12
I


# In[28]:


I_A_Ao = I*me.outer(A.y, A.y) + I*me.outer(A.z, A.z)
I_B_Bo = I*me.outer(B.x, B.x) + I*me.outer(B.z, B.z)
I_A_Ao, I_B_Bo


# In[29]:


Rs_Ao = -m*Ao.acc(N)
Rs_Bo = -m*Bo.acc(N)

Rs_Ao, Rs_Bo


# In[30]:


Ts_A = -(A.ang_acc_in(N).dot(I_A_Ao) +
         me.cross(A.ang_vel_in(N), I_A_Ao).dot(A.ang_vel_in(N)))
Ts_A


# In[31]:


Ts_B = -(B.ang_acc_in(N).dot(I_B_Bo) +
         me.cross(B.ang_vel_in(N), I_B_Bo).dot(B.ang_vel_in(N)))
Ts_B


# In[32]:


F1s_A = v_Ao_1.dot(Rs_Ao) + w_A_1.dot(Ts_A)
F1s_B = v_Bo_1.dot(Rs_Bo) + w_B_1.dot(Ts_B)
F2s_A = v_Ao_2.dot(Rs_Ao) + w_A_2.dot(Ts_A)
F2s_B = v_Bo_2.dot(Rs_Bo) + w_B_2.dot(Ts_B)


# In[33]:


F1s = F1s_A + F1s_B
F2s = F2s_A + F2s_B

Frs = sm.Matrix([F1s, F2s])
Frs

