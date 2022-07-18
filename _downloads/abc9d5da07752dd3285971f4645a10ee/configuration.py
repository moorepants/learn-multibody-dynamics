#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sympy as sm
import sympy.physics.mechanics as me
me.init_vprinting(use_latex='mathjax')


# In[2]:


q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')
la, lb, lc, ln = sm.symbols('l_a, l_b, l_c, l_n')

N = me.ReferenceFrame('N')
A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
C = me.ReferenceFrame('C')

A.orient_axis(N, q1, N.z)
B.orient_axis(A, q2, A.z)
C.orient_axis(B, q3, B.z)

P1 = me.Point('P1')
P2 = me.Point('P2')
P3 = me.Point('P3')
P4 = me.Point('P4')


# In[3]:


P2.set_pos(P1, la*A.x)
P3.set_pos(P2, lb*B.x)
P4.set_pos(P3, lc*C.x)

P4.pos_from(P1)


# In[4]:


r_P1_P4 = ln*N.x


# In[5]:


loop = P4.pos_from(P1) - r_P1_P4
loop


# In[6]:


fhx = sm.trigsimp(loop.dot(N.x))
fhx


# In[7]:


fhy = sm.trigsimp(loop.dot(N.y))
fhy


# In[8]:


fh = sm.Matrix([fhx, fhy])
fh

