#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sympy as sm
import sympy.physics.mechanics as me
sm.init_printing(use_latex='mathjax')


# In[2]:


N = me.ReferenceFrame('N')


# In[3]:


a, b, c, d, e, f = sm.symbols('a, b, c, d, e, f')


# In[4]:


v = a*N.x
v


# In[5]:


v.to_matrix(N)


# In[6]:


w = a*N.x + b*N.y + c*N.z
w


# In[7]:


w.to_matrix(N)


# In[8]:


x = d*N.x + e*N.y + f*N.z
x


# In[9]:


w + x


# In[10]:


y = 2*w
y


# In[11]:


z = -w
z


# In[12]:


N = me.ReferenceFrame('N')
l, theta = sm.symbols('l, theta')


# In[13]:


v1 = l*sm.cos(sm.pi/4)*N.x + l*sm.sin(sm.pi/4)*N.y
v1


# In[14]:


v2 = -10*N.y
v2


# In[15]:


v3 = -l*sm.sin(theta)*N.x + l*sm.cos(theta)*N.y
v3


# In[16]:


v1 + v2 - 5*v3


# In[17]:


N = me.ReferenceFrame('N')
w = a*N.x + b*N.y + c*N.z
x = d*N.x + e*N.y + f*N.z


# In[18]:


me.dot(w, x)


# In[19]:


w.dot(x)


# In[20]:


w.normalize()


# In[21]:


def normalize(vector):
    return vector/sm.sqrt(me.dot(vector, vector))

normalize(w)


# In[22]:


w.magnitude()


# In[23]:


w/w.magnitude()


# In[24]:


N = me.ReferenceFrame('N')
v1 = a*N.x + b*N.y + a*N.z
v2 = b*N.x + a*N.y + b*N.z


# In[25]:


sm.acos(v1.dot(v2) / (v1.magnitude()*v2.magnitude()))


# In[26]:


N = me.ReferenceFrame('N')
w = a*N.x + b*N.y + c*N.z
w


# In[27]:


x = d*N.x + e*N.y + f*N.z
x


# In[28]:


me.cross(w, x)


# In[29]:


w.cross(x)


# In[30]:


N = me.ReferenceFrame('N')

p1 = 23*N.x - 12* N.y
p2 = 16*N.x + 2*N.y - 4*N.z
p3 = N.x + 14*N.z

me.cross(p2 - p1, p3 - p1).magnitude() / 2


# In[31]:


N = me.ReferenceFrame('N')
A = me.ReferenceFrame('A')
a, b, theta = sm.symbols('a, b, theta')

v = a*A.x + b*N.y
v


# In[32]:


v + v


# In[33]:


A.orient_axis(N, theta, N.z)

v.express(N)


# In[34]:


v.express(A)


# In[35]:


q1, q2, q3, q4, q5 = sm.symbols('q1, q2, q3, q4, q5')
l1, l2, l3, l4 = sm.symbols('l1, l2, l3, l4')
N = me.ReferenceFrame('N')
A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
C = me.ReferenceFrame('C')


# In[36]:


A.orient_body_fixed(N, (q1, q2, 0), 'ZXZ')


# In[37]:


B.orient_axis(A, q3, A.x)


# In[38]:


C.orient_body_fixed(B, (q4, q5, 0), 'XZX')


# In[39]:


R_P1_P2 = l1*A.z
R_P2_P3 = l2*B.z
R_P3_P4 = l3*C.z - l4*C.y


# In[40]:


R_P1_P4 = R_P1_P2 + R_P2_P3 + R_P3_P4
R_P1_P4


# In[41]:


R_P1_P4.express(N)


# In[42]:


R_P1_P2.express(N)


# In[43]:


R_P1_P2.free_symbols(N)


# In[44]:


R_P1_P2.free_symbols(A)


# In[45]:


R_P1_P4.free_symbols(N)

