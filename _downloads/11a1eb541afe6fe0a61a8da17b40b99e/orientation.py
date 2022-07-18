#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sympy as sm
sm.init_printing(use_latex='mathjax')


# In[2]:


theta = sm.symbols('theta')

A_C_N = sm.Matrix([[sm.cos(theta), sm.sin(theta), 0],
                   [-sm.sin(theta), sm.cos(theta), 0],
                   [0, 0, 1]])
A_C_N


# In[3]:


sm.trigsimp(A_C_N.inv())


# In[4]:


A_C_N.transpose()


# In[5]:


A_C_N


# In[6]:


alpha = sm.symbols('alpha')

B_C_A = sm.Matrix([[sm.cos(alpha), sm.sin(alpha), 0],
                   [-sm.sin(alpha), sm.cos(alpha), 0],
                   [0, 0, 1]])

B_C_A


# In[7]:


B_C_N = B_C_A*A_C_N
B_C_N


# In[8]:


sm.trigsimp(B_C_N)


# In[9]:


import sympy.physics.mechanics as me


# In[10]:


N = me.ReferenceFrame('N')


# In[11]:


N.x, N.y, N.z


# In[12]:


A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')

N, A, B


# In[13]:


A_C_N


# In[14]:


N.orient_explicit(A, A_C_N)


# In[15]:


A.dcm(N)


# In[16]:


B.orient_axis(A, alpha, A.z)


# In[17]:


B.dcm(A)


# In[18]:


A.dcm(B)


# In[19]:


sm.trigsimp(B.dcm(A)*A.dcm(N))


# In[20]:


sm.trigsimp(B.dcm(N))


# In[21]:


sm.trigsimp(me.dot(B.x, N.x))


# In[22]:


psi = sm.symbols('psi')

A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')

B.orient_axis(A, psi, A.z)

B.dcm(A)


# In[23]:


theta = sm.symbols('theta')

C = me.ReferenceFrame('C')

C.orient_axis(B, theta, B.x)

C.dcm(B)


# In[24]:


phi = sm.symbols('varphi')

D = me.ReferenceFrame('D')

D.orient_axis(C, phi, C.y)

D.dcm(C)


# In[25]:


D.dcm(A)


# In[26]:


A = me.ReferenceFrame('A')
D = me.ReferenceFrame('D')

D.orient_body_fixed(A, (psi, theta, phi), 'zxy')

D.dcm(A)


# In[27]:


A = me.ReferenceFrame('A')
D = me.ReferenceFrame('D')

D.orient_body_fixed(A, (psi, theta, phi), 'zxz')

D.dcm(A)

