#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sympy as sm
import sympy.physics.mechanics as me
sm.init_printing(use_latex='mathjax')


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


alpha, beta = sm.symbols('alpha, beta')
a, b, c, d, e, f = sm.symbols('a, b, c, d, e, f')

A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
C = me.ReferenceFrame('C')

B.orient_axis(A, alpha, A.x)
C.orient_axis(B, beta, B.y)

v = a*A.x + b*A.y + c*B.x + d*B.y + e*C.x + f*C.y
v


# In[4]:


dvdalphaAx = v.dot(A.x).diff(alpha)
dvdalphaAx


# In[5]:


dvdalphaAy = v.dot(A.y).diff(alpha)
dvdalphaAy


# In[6]:


dvdalphaAz = v.dot(A.z).diff(alpha)
dvdalphaAz


# In[7]:


dvdalphaA = dvdalphaAx*A.x + dvdalphaAy*A.y + dvdalphaAz*A.z
dvdalphaA


# In[8]:


v.diff(alpha, A)


# In[9]:


dvdeBx = v.dot(B.x).diff(e)
dvdeBy = v.dot(B.y).diff(e)
dvdeBz = v.dot(B.z).diff(e)
dvdeBx*B.x + dvdeBy*B.y + dvdeBz*B.z


# In[10]:


v.diff(e, B).express(B)


# In[11]:


t = sm.symbols('t')
q_of = sm.Function('q')

q = q_of(t)
q


# In[12]:


q.diff(t)


# In[13]:


q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')
q1, q2, q3


# In[14]:


t = me.dynamicsymbols._t


# In[15]:


me.init_vprinting(use_latex='mathjax')
q1.diff(t), q2.diff(t, 2), q3.diff(t, 3)


# In[16]:


A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
B.orient_body_fixed(A, (q1, q2, q3), 'ZXZ')
v = q1*A.x + q2*A.y + t**2*A.z
v


# In[17]:


v.diff(t, A)


# In[18]:


v.dt(A)

