#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sympy as sm
import sympy.physics.mechanics as me
sm.init_printing(use_latex='mathjax')


# In[2]:


alpha, beta = sm.symbols('alpha, beta')
a, b, c, d, e, f = sm.symbols('a, b, c, d, e, f')

A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
C = me.ReferenceFrame('C')

B.orient_axis(A, alpha, A.x)
C.orient_axis(B, beta, B.y)

v = a*A.x + b*A.y + c*B.x + d*B.y + e*C.x + f*C.y
v


# In[3]:


dvdax = v.dot(A.x).diff(alpha)
dvdax


# In[4]:


dvday = v.dot(A.y).diff(alpha)
dvday


# In[5]:


dvdaz = v.dot(A.z).diff(alpha)
dvdaz


# In[6]:


dvda = dvdax*A.x + dvday*A.y + dvdaz*A.z
dvda


# In[7]:


dvdalpha = v.diff(alpha, A)
dvdalpha


# In[8]:


v.diff(alpha, A).simplify()


# In[9]:


v.diff(alpha, A).express(A).simplify()


# In[10]:


t = sm.symbols('t')
q_of = sm.Function('q')

q = q_of(t)
q


# In[11]:


q.diff(t)


# In[12]:


q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')
q1, q2, q3


# In[13]:


t = me.dynamicsymbols._t


# In[14]:


me.init_vprinting(use_latex='mathjax')
q1.diff(t), q2.diff(t, 2), q3.diff(t, 3)


# In[15]:


A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
B.orient_body_fixed(A, (q1, q2, q3), 'ZXZ')
v = q1*A.x + q2*A.y + t**2*A.z
v


# In[16]:


v.diff(t, A)


# In[17]:


v.dt(A)

