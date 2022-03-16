#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sympy as sm
import sympy.physics.mechanics as me
me.init_vprinting(use_latex='mathjax')


# In[2]:


cxx, cyy, czz = me.dynamicsymbols('c_{xx}, c_{yy}, c_{zz}')
cxy, cxz, cyx = me.dynamicsymbols('c_{xy}, c_{xz}, c_{yx}')
cyz, czx, czy = me.dynamicsymbols('c_{yz}, c_{zx}, c_{zy}')

B_C_A = sm.Matrix([[cxx, cxy, cxz],
                   [cyx, cyy, cyz],
                   [czx, czy, czz]])


# In[3]:


A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
B.orient_explicit(A, B_C_A.transpose())
B.dcm(A)


# In[4]:


B.x.express(A)


# In[5]:


B.y.express(A)


# In[6]:


B.z.express(A)


# In[7]:


B.y.express(A).dt(A)


# In[8]:


mnx = me.dot(B.y.express(A).dt(A), B.z)
mnx


# In[9]:


mny = me.dot(B.z.express(A).dt(A), B.x)
mny


# In[10]:


mnz = me.dot(B.x.express(A).dt(A), B.y)
mnz


# In[11]:


A_w_B = mnx*B.x + mny*B.y + mnz*B.z
A_w_B


# In[12]:


theta = me.dynamicsymbols('theta')

B_C_A = sm.Matrix([[sm.cos(theta), sm.sin(theta), 0],
                   [-sm.sin(theta), sm.cos(theta), 0],
                   [0, 0, 1]])
B_C_A


# In[13]:


A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
A.orient_explicit(B, B_C_A)

mnx = me.dot(B.y.express(A).dt(A), B.z)
mny = me.dot(B.z.express(A).dt(A), B.x)
mnz = me.dot(B.x.express(A).dt(A), B.y)

A_w_B = mnx*B.x + mny*B.y + mnz*B.z
A_w_B.simplify()


# In[14]:


A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
B.orient_axis(A, theta, A.z)
B.ang_vel_in(A)


# In[15]:


theta = me.dynamicsymbols('theta')

A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
B.orient_axis(A, theta, A.x + A.y)
B.ang_vel_in(A)


# In[16]:


B.ang_vel_in(A).magnitude()


# In[17]:


psi, theta, phi = me.dynamicsymbols('psi, theta, varphi')

A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
B.orient_body_fixed(A, (psi, theta, phi), 'ZXZ')

mnx = me.dot(B.y.express(A).dt(A), B.z)
mny = me.dot(B.z.express(A).dt(A), B.x)
mnz = me.dot(B.x.express(A).dt(A), B.y)

A_w_B = mnx*B.x + mny*B.y + mnz*B.z
A_w_B.simplify()


# In[18]:


B.ang_vel_in(A).simplify()


# In[19]:


A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
B.orient_body_fixed(A, (psi, theta, 0), 'ZXZ')

u1, u2, u3 = me.dynamicsymbols('u1, u2, u3')

u = u1*B.x + u2*B.y + u3*B.z
u


# In[20]:


u.express(A)


# In[21]:


u.express(A).dt(A)


# In[22]:


u.dt(B)


# In[23]:


A_w_B = B.ang_vel_in(A)
A_w_B


# In[24]:


u.dt(B) + me.cross(A_w_B, u)


# In[25]:


u.express(A).dt(A).express(B).simplify()


# In[26]:


psi, theta, phi = me.dynamicsymbols('psi, theta, varphi')

A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
C = me.ReferenceFrame('C')
D = me.ReferenceFrame('D')

B.orient_axis(A, psi, A.y)
C.orient_axis(B, theta, B.x)
D.orient_axis(C, phi, C.y)


# In[27]:


A_w_B = B.ang_vel_in(A)
A_w_B


# In[28]:


B_w_C = C.ang_vel_in(B)
B_w_C


# In[29]:


C_w_D = D.ang_vel_in(C)
C_w_D


# In[30]:


A_w_D = A_w_B + B_w_C + C_w_D
A_w_D


# In[31]:


A2 = me.ReferenceFrame('A')
D2 = me.ReferenceFrame('D')
D2.orient_body_fixed(A2, (psi, theta, phi), 'YXY')
D2.ang_vel_in(A2).simplify()


# In[32]:


A_w_D.express(D)


# In[33]:


theta = me.dynamicsymbols('theta')

A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
B.orient_axis(A, theta, A.z)
B.ang_acc_in(A)


# In[34]:


B.ang_vel_in(A).dt(A)


# In[35]:


B.ang_vel_in(A).dt(B)


# In[36]:


psi, theta, phi = me.dynamicsymbols('psi, theta, varphi')

A = me.ReferenceFrame('A')
D = me.ReferenceFrame('D')
D.orient_body_fixed(A, (psi, theta, phi), 'YXY')

D.ang_acc_in(A).simplify()


# In[37]:


D.ang_vel_in(A).dt(A).simplify()


# In[38]:


D.ang_vel_in(A).dt(D).simplify()


# In[39]:


psi, theta, phi = me.dynamicsymbols('psi, theta, varphi')

A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
C = me.ReferenceFrame('C')
D = me.ReferenceFrame('D')

B.orient_axis(A, psi, A.y)
C.orient_axis(B, theta, B.x)
D.orient_axis(C, phi, C.y)


# In[40]:


A_alp_B = B.ang_acc_in(A)
A_alp_B


# In[41]:


B_alp_C = C.ang_acc_in(B)
B_alp_C


# In[42]:


C_alp_D = D.ang_acc_in(C)
C_alp_D


# In[43]:


A_alp_D = A_alp_B + B_alp_C + C_alp_D
A_alp_D.express(D).simplify()


# In[44]:


D.ang_acc_in(A).express(D).simplify()

