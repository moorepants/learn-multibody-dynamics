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


N = me.ReferenceFrame('N')


# In[4]:


a, b, c, d, e, f = sm.symbols('a, b, c, d, e, f')


# In[5]:


v = a*N.x
v


# In[6]:


v.to_matrix(N)


# In[7]:


w = a*N.x + b*N.y + c*N.z
w


# In[8]:


w.to_matrix(N)


# In[9]:


x = d*N.x + e*N.y + f*N.z
x


# In[10]:


w + x


# In[11]:


y = 2*w
y


# In[12]:


z = -w
z


# In[13]:


N = me.ReferenceFrame('N')
l, theta = sm.symbols('l, theta')


# In[14]:


v1 = l*sm.cos(sm.pi/4)*N.x + l*sm.sin(sm.pi/4)*N.y
v1


# In[15]:


v2 = -10*N.y
v2


# In[16]:


v3 = -l*sm.sin(theta)*N.x + l*sm.cos(theta)*N.y
v3


# In[17]:


v1 + v2 - 5*v3


# In[18]:


N = me.ReferenceFrame('N')
w = a*N.x + b*N.y + c*N.z
x = d*N.x + e*N.y + f*N.z


# In[19]:


me.dot(w, x)


# In[20]:


w.dot(x)


# In[21]:


w.normalize()


# In[22]:


def normalize(vector):
    return vector/sm.sqrt(me.dot(vector, vector))

normalize(w)


# In[23]:


w.magnitude()


# In[24]:


w/w.magnitude()


# In[25]:


N = me.ReferenceFrame('N')
v1 = a*N.x + b*N.y + a*N.z
v2 = b*N.x + a*N.y + b*N.z


# In[26]:


sm.acos(v1.dot(v2) / (v1.magnitude()*v2.magnitude()))


# In[27]:


N = me.ReferenceFrame('N')
w = a*N.x + b*N.y + c*N.z
w


# In[28]:


x = d*N.x + e*N.y + f*N.z
x


# In[29]:


me.cross(w, x)


# In[30]:


w.cross(x)


# In[31]:


N = me.ReferenceFrame('N')

p1 = 23*N.x - 12* N.y
p2 = 16*N.x + 2*N.y - 4*N.z
p3 = N.x + 14*N.z

me.cross(p2 - p1, p3 - p1).magnitude() / 2


# In[32]:


N = me.ReferenceFrame('N')
A = me.ReferenceFrame('A')

a, b, theta = sm.symbols('a, b, theta')

v = a*A.x + b*N.y
v


# In[33]:


v + v


# In[34]:


A.orient_axis(N, theta, N.z)

v.express(N)


# In[35]:


v.express(A)


# In[36]:


q1, q2, q3, q4, q5 = sm.symbols('q1, q2, q3, q4, q5')
l1, l2, l3, l4 = sm.symbols('l1, l2, l3, l4')
N = me.ReferenceFrame('N')
A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
C = me.ReferenceFrame('C')


# In[37]:


A.orient_body_fixed(N, (q1, q2, 0), 'ZXZ')


# In[38]:


B.orient_axis(A, q3, A.x)


# In[39]:


C.orient_body_fixed(B, (q4, q5, 0), 'XZX')


# In[40]:


R_P1_P2 = l1*A.z
R_P2_P3 = l2*B.z
R_P3_P4 = l3*C.z - l4*C.y


# In[41]:


R_P1_P4 = R_P1_P2 + R_P2_P3 + R_P3_P4
R_P1_P4


# In[42]:


R_P1_P4.express(N)


# In[43]:


R_P1_P2.express(N)


# In[44]:


R_P1_P2.free_symbols(N)


# In[45]:


R_P1_P2.free_symbols(A)


# In[46]:


R_P1_P4.free_symbols(N)

