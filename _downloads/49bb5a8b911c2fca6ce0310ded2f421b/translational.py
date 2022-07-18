#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sympy as sm
import sympy.physics.mechanics as me
me.init_vprinting(use_latex='mathjax')


# In[2]:


alpha, beta = me.dynamicsymbols('alpha, beta')

N = me.ReferenceFrame('N')
A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')

A.orient_axis(N, alpha, N.z)
B.orient_axis(A, beta, A.x)


# In[3]:


h, d, w, c, l = sm.symbols('h, d, w, c, l')

r_O_P = h*N.z
r_P_S = -d*A.x
r_S_Q = -w*B.x - (c + l/2)*B.z

r_O_P, r_P_S, r_S_Q


# In[4]:


(r_O_P + r_P_S).dt(A)


# In[5]:


A.ang_vel_in(N)


# In[6]:


me.cross(A.ang_vel_in(N), r_O_P + r_P_S)


# In[7]:


N_v_S = (r_O_P + r_P_S).dt(A) + me.cross(A.ang_vel_in(N), r_O_P + r_P_S)
N_v_S


# In[8]:


(r_O_P + r_P_S + r_S_Q).dt(B)


# In[9]:


me.cross(B.ang_vel_in(N), r_O_P + r_P_S + r_S_Q)


# In[10]:


N_v_Q = (r_O_P + r_P_S + r_S_Q).dt(B) + me.cross(B.ang_vel_in(N), r_O_P + r_P_S + r_S_Q)
N_v_Q


# In[11]:


O = me.Point('O')
P = me.Point('P')
S = me.Point('S')
Q = me.Point('Q')

P.set_pos(O, h*N.z)
S.set_pos(P, -d*A.x)
Q.set_pos(S, -w*B.x - (c + l/2)*B.z)


# In[12]:


Q.pos_from(O)


# In[13]:


O.set_vel(N, 0)


# In[14]:


Q.vel(N)


# In[15]:


N_v_P = 0*N.z


# In[16]:


N_v_S = N_v_P +  me.cross(A.ang_vel_in(N), S.pos_from(P))
N_v_S


# In[17]:


P.set_vel(N, 0)
S.v2pt_theory(P, N, A)


# In[18]:


S.vel(N)


# In[19]:


N_v_Q = N_v_S +  me.cross(B.ang_vel_in(N), Q.pos_from(S))
N_v_Q


# In[20]:


Q.v2pt_theory(S, N, B)


# In[21]:


Bc = me.Point('B_c')
Bc.set_pos(S, -c*B.z - w/2*A.x)
Bc.v2pt_theory(S, N, B)


# In[22]:


s = me.dynamicsymbols('s')
t = me.dynamicsymbols._t

R = me.Point('R')
R.set_pos(Q, l*B.z + s*B.x)


# In[23]:


B_v_R = s.diff(t)*B.x
B_v_R


# In[24]:


r_S_R = R.pos_from(S)
r_S_R


# In[25]:


N_v_T = N_v_S + me.cross(B.ang_vel_in(N), r_S_R)
N_v_T


# In[26]:


N_v_R = B_v_R + N_v_T
N_v_R


# In[27]:


S.set_vel(B, 0)
R.v1pt_theory(S, N, B)


# In[28]:


S.acc(N)


# In[29]:


me.cross(A.ang_acc_in(N), S.pos_from(P))


# In[30]:


me.cross(A.ang_vel_in(N), me.cross(A.ang_vel_in(N), S.pos_from(P)))


# In[31]:


S.a2pt_theory(P, N, A)


# In[32]:


Q.a2pt_theory(S, N, B)


# In[33]:


B_a_R = R.acc(B)
B_a_R


# In[34]:


N_a_T = R.a2pt_theory(S, N, B)
N_a_T


# In[35]:


2*me.cross(B.ang_vel_in(N), R.vel(B))


# In[36]:


R.a1pt_theory(S, N, B)

