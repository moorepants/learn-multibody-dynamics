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


alpha, beta = me.dynamicsymbols('alpha, beta')

N = me.ReferenceFrame('N')
A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')

A.orient_axis(N, alpha, N.z)
B.orient_axis(A, beta, A.x)


# In[4]:


h, d, w, c, l = sm.symbols('h, d, w, c, l')

r_O_P = h*N.z
r_P_S = -d*A.x
r_S_Q = -w*B.x - (c + l/2)*B.z

r_O_P, r_P_S, r_S_Q


# In[5]:


(r_O_P + r_P_S).dt(A)


# In[6]:


A.ang_vel_in(N)


# In[7]:


me.cross(A.ang_vel_in(N), r_O_P + r_P_S)


# In[8]:


N_v_S = (r_O_P + r_P_S).dt(A) + me.cross(A.ang_vel_in(N), r_O_P + r_P_S)
N_v_S


# In[9]:


(r_P_S + r_S_Q).dt(B)


# In[10]:


me.cross(B.ang_vel_in(N), r_P_S + r_S_Q)


# In[11]:


N_v_Q = (r_P_S + r_S_Q).dt(B) + me.cross(B.ang_vel_in(N), r_P_S + r_S_Q)
N_v_Q


# In[12]:


O = me.Point('O')
P = me.Point('P')
S = me.Point('S')
Q = me.Point('Q')

P.set_pos(O, h*N.z)
S.set_pos(P, -d*A.x)
Q.set_pos(S, -w*B.x - (c + l/2)*B.z)


# In[13]:


Q.pos_from(O)


# In[14]:


O.set_vel(N, 0)


# In[15]:


Q.vel(N)


# In[16]:


Bc = me.Point('B_c')
Bc.set_pos(S, -c*B.z - w/2*A.x)
me.cross(B.ang_vel_in(A), Bc.pos_from(S))


# In[17]:


N_v_P = 0*N.z


# In[18]:


N_v_S = N_v_P +  me.cross(A.ang_vel_in(N), S.pos_from(P))
N_v_S


# In[19]:


P.set_vel(N, 0)
S.v2pt_theory(P, N, A)


# In[20]:


S.vel(N)


# In[21]:


N_v_Q = N_v_S +  me.cross(B.ang_vel_in(N), Q.pos_from(S))
N_v_Q


# In[22]:


Q.v2pt_theory(S, N, B)


# In[23]:


Bc = me.Point('B_c')
Bc.set_pos(S, -c*B.z - w/2*A.x)
Bc.v2pt_theory(S, N, B)


# In[24]:


s = me.dynamicsymbols('s')
t = me.dynamicsymbols._t

R = me.Point('R')
R.set_pos(Q, l*B.z + s*B.x)


# In[25]:


B_v_R = s.diff(t)*B.x
B_v_R


# In[26]:


r_S_R = R.pos_from(S)
r_S_R


# In[27]:


N_v_T = N_v_S + me.cross(B.ang_vel_in(N), r_S_R)
N_v_T


# In[28]:


N_v_R = B_v_R + N_v_T
N_v_R


# In[29]:


S.set_vel(B, 0)
R.v1pt_theory(S, N, B)


# In[30]:


S.acc(N)


# In[31]:


me.cross(A.ang_acc_in(N), S.pos_from(P))


# In[32]:


me.cross(A.ang_vel_in(N), me.cross(A.ang_vel_in(N), S.pos_from(P)))


# In[33]:


S.a2pt_theory(P, N, A)


# In[34]:


Q.a2pt_theory(S, N, B)


# In[35]:


B_a_R = R.acc(B)
B_a_R


# In[36]:


N_a_T = R.a2pt_theory(S, N, B)
N_a_T


# In[37]:


2*me.cross(B.ang_vel_in(N), R.vel(B))


# In[38]:


R.a1pt_theory(S, N, B)

