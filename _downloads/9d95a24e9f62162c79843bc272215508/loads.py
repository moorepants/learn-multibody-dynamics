#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sympy as sm
import sympy.physics.mechanics as me
me.init_vprinting(use_latex='mathjax')


# In[2]:


N = me.ReferenceFrame('N')

F1 = 2*N.x + 3*N.y
F2 = -4*N.x + 5*N.y

r_O_P1 = 2*N.x
r_O_P2 = 3*N.x


# In[3]:


r_O_P = -5*N.x

M_S_P = me.cross(r_O_P1 - r_O_P, F1) + me.cross(r_O_P2 - r_O_P, F2)
M_S_P


# In[4]:


r_O_Q = 5*N.y
M_S_Q = me.cross(r_O_P1 - r_O_Q, F1) + me.cross(r_O_P2 - r_O_Q, F2)

M_S_P = M_S_Q + me.cross(r_O_Q - r_O_P, F1 + F2)
M_S_P


# In[5]:


l, w = sm.symbols('l, w')
Ffl, Ffr, Frl, Frr = me.dynamicsymbols('F_{fl}, F_{fr}, F_{rl}, F_{rr}')
alphafl, alphafr = me.dynamicsymbols(r'\alpha_{fl}, \alpha_{fr}')
alpharl, alpharr = me.dynamicsymbols(r'\alpha_{rl}, \alpha_{rr}')
delta = me.dynamicsymbols('delta')


# In[6]:


B = me.ReferenceFrame('B')
W = me.ReferenceFrame('W')
FR = me.ReferenceFrame('F_R')
FL = me.ReferenceFrame('F_L')
RR = me.ReferenceFrame('R_R')
RL = me.ReferenceFrame('R_L')

W.orient_axis(B, delta, B.z)
FR.orient_axis(W, alphafr, W.z)
FL.orient_axis(W, alphafl, W.z)
RR.orient_axis(B, alpharr, B.z)
RL.orient_axis(B, alpharl, B.z)


# In[7]:


R = Ffl*FL.x + Ffr*FR.x + Frl*RL.x + Frr*RR.x
R.express(B).simplify()


# In[8]:


T = (me.cross(l/2*B.x - w/2*B.y, Ffl*FL.x) +
     me.cross(l/2*B.x + w/2*B.y, Ffr*FR.x) +
     me.cross(-l/2*B.x - w/2*B.y, Frl*RL.x) +
     me.cross(-l/2*B.x + w/2*B.y, Frr*RR.x))
T = T.express(B).simplify()
T


# In[9]:


Bo = me.Point('Bo')
force = (Bo, R)
force


# In[10]:


torque = (B, T)
torque


# In[11]:


T, q = me.dynamicsymbols('T, q')

N = me.ReferenceFrame('N')
B = me.ReferenceFrame('B')

Tm = T*N.z


# In[12]:


(B, Tm), (N, -Tm)


# In[13]:


m, g = sm.symbols('m, g')
Fg = -m*g*N.y
Fg


# In[14]:


q0, k = sm.symbols('q0, k')
q1, q2 = me.dynamicsymbols('q1, q2')

displacement = q2 - q1 - q0
displacement


# In[15]:


Fs = -k*displacement*N.x
Fs


# In[16]:


c = sm.symbols('c')
t = me.dynamicsymbols._t

Fc = -c*displacement.diff(t)*N.x
Fc


# In[17]:


mu, m, g = sm.symbols('mu, m, g')

Fn = m*g

displacement = q2 - q1

Ff = sm.Piecewise((mu*Fn, displacement.diff(t) < 0),
                  (-mu*Fn, displacement.diff(t) > 0),
                  (0, True))*N.x
Ff


# In[18]:


Ff = -mu*Fn*sm.sign(displacement.diff(t))*N.x
Ff


# In[19]:


A, Cd, rho = sm.symbols('A, C_d, rho')
ux, uy, uz = me.dynamicsymbols('u_x, u_y, u_z', real=True)

N_v_P = ux*N.x + uy*N.y + uz*N.z

Fd = -N_v_P.normalize()*Cd*A*rho/2*N_v_P.dot(N_v_P)
Fd


# In[20]:


Fd.xreplace({uy: 0, uz:0})


# In[21]:


x, y, z = me.dynamicsymbols('x, y, z', real=True)

r_O_P = x*N.x + y*N.y + z*N.z

zh = r_O_P.dot(N.z)

zp = (sm.Abs(zh) - zh)/2
zp


# In[22]:


k, c = sm.symbols('k, c')

Fz = (k*zp**3 + c*sm.Piecewise((zh.diff(), zh < 0), (0, True)))*N.z
Fz


# In[23]:


mu = sm.symbols('mu')

vx = r_O_P.dot(N.x).diff(t)
vy = r_O_P.dot(N.y).diff(t)

Fx = -sm.Abs(vx)/vx*mu*Fz.dot(N.z)*N.x
Fx


# In[24]:


Fy = -sm.Abs(vy)/vy*mu*Fz.dot(N.z)*N.y
Fy


# In[25]:


vz = me.dynamicsymbols('v_z', negative=True)

repl = {z.diff(): vz, z: 0}

Fx.xreplace(repl), Fy.xreplace(repl), Fz.xreplace(repl)


# In[26]:


repl = {z.diff(): vz, z: 2}

Fx.xreplace(repl), Fy.xreplace(repl), Fz.xreplace(repl)


# In[27]:


repl = {z.diff(): vz, z: -2}

Fx.xreplace(repl), Fy.xreplace(repl), Fz.xreplace(repl)


# In[28]:


Fc = Fx + Fy + Fz
Fc

