#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sympy as sm
import sympy.physics.mechanics as me
me.init_vprinting(use_latex='mathjax')


# In[2]:


m1, m2, m3 = sm.symbols("m1, m2, m3")
x1, x2, x3 = me.dynamicsymbols('x1, x2, x3')
y1, y2, y3 = me.dynamicsymbols('y1, y2, y3')
z1, z2, z3 = me.dynamicsymbols('z1, z2, z3')

A = me.ReferenceFrame('A')

r_O_So = (m1*(x1*A.x + y1*A.y + z1*A.z) +
          m2*(x2*A.x + y2*A.y + z2*A.z) +
          m3*(x3*A.x + y3*A.y + z3*A.z)) / (m1 + m2 + m3)
r_O_So


# In[3]:


r_O_So.xreplace({m2: 2*m1, m3: 3*m1}).simplify()


# In[4]:


v1, v2, v3 = sm.symbols('v1, v2, v3')
w1, w2, w3 = sm.symbols('w1, w2, w3')

A = me.ReferenceFrame('A')

v = v1*A.x + v2*A.y + v3*A.z
w = w1*A.x + w2*A.y + w3*A.z

Q = me.outer(v, w)
Q


# In[5]:


Q.to_matrix(A)


# In[6]:


me.outer(A.x, A.x)


# In[7]:


me.outer(A.x, A.x).to_matrix(A)


# In[8]:


me.outer(A.y, A.z)


# In[9]:


me.outer(A.y, A.z).to_matrix(A)


# In[10]:


theta = sm.symbols("theta")

A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')

B.orient_axis(A, theta, A.x)

P = 2*me.outer(B.x, B.x) + 3*me.outer(A.x, B.y) + 4*me.outer(B.z, A.z)
P


# In[11]:


P.express(A)


# In[12]:


P.to_matrix(A)


# In[13]:


P.express(B)


# In[14]:


P.to_matrix(B)


# In[15]:


U = me.outer(A.x, A.x) + me.outer(A.y, A.y) + me.outer(A.z, A.z)
U


# In[16]:


U.to_matrix(A)


# In[17]:


U.express(B).simplify()


# In[18]:


Ixx, Iyy, Izz = sm.symbols('I_{xx}, I_{yy}, I_{zz}')
Ixy, Iyz, Ixz = sm.symbols('I_{xy}, I_{yz}, I_{xz}')

I = me.inertia(A, Ixx, Iyy, Izz, ixy=Ixy, iyz=Iyz, izx=Ixz)
I


# In[19]:


I.to_matrix(A)


# In[20]:


sm.trigsimp(I.to_matrix(B))


# In[21]:


sm.trigsimp(B.dcm(A)*I.to_matrix(A)*A.dcm(B))


# In[22]:


Ixx, Iyy, Izz = sm.symbols('I_{xx}, I_{yy}, I_{zz}')
Ixy, Iyz, Ixz = sm.symbols('I_{xy}, I_{yz}, I_{xz}')
w1, w2, w3 = me.dynamicsymbols('omega1, omega2, omega3')

B = me.ReferenceFrame('B')

I = me.inertia(B, Ixx, Iyy, Izz, Ixy, Iyz, Ixz)

A_w_B = w1*B.x + w2*B.y + w3*B.z

I.dot(A_w_B)


# In[23]:


dx, dy, dz, m = sm.symbols('d_x, d_y, d_z, m')

N = me.ReferenceFrame('N')

r_O_Bo = dx*N.x + dy*N.y + dz*N.z

U = me.outer(N.x, N.x) + me.outer(N.y, N.y) + me.outer(N.z, N.z)

I_Bo_O = m*(me.dot(r_O_Bo, r_O_Bo)*U - me.outer(r_O_Bo, r_O_Bo))
I_Bo_O


# In[24]:


I_Bo_O.to_matrix(N)


# In[25]:


I = sm.Matrix([[1.0451, 0.0, -0.1123],
               [0.0, 2.403, 0.0],
               [-0.1123, 0.0, 1.8501]])
I


# In[26]:


I.eigenvects()

