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


# In[11]:


N = me.ReferenceFrame('N')


# In[12]:


N.x, N.y, N.z


# In[13]:


A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')

N, A, B


# In[14]:


A_C_N


# In[15]:


N.orient_explicit(A, A_C_N)


# In[16]:


A.dcm(N)


# In[17]:


N.dcm(A)


# In[18]:


beta = sm.symbols('beta')

D = me.ReferenceFrame('D')
F = me.ReferenceFrame('F')

F_C_D = sm.Matrix([[sm.cos(beta), 0, -sm.sin(beta)],
                   [0, 1, 0],
                   [sm.sin(beta), 0, sm.cos(beta)]])

F.orient_explicit(D, F_C_D.transpose())

F.dcm(D)


# In[19]:


B.orient_axis(A, alpha, A.z)


# In[20]:


B.dcm(A)


# In[21]:


A.dcm(B)


# In[22]:


sm.trigsimp(B.dcm(A)*A.dcm(N))


# In[23]:


sm.trigsimp(B.dcm(N))


# In[24]:


sm.trigsimp(me.dot(B.x, N.x))


# In[25]:


beta = sm.symbols('beta')

C = me.ReferenceFrame('C')
D = me.ReferenceFrame('D')
E = me.ReferenceFrame('E')

D.orient_axis(C, beta, -C.y)

D.dcm(C)


# In[26]:


E.orient_explicit(D, C.dcm(D))
E.dcm(D)


# In[27]:


psi = sm.symbols('psi')

A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')

B.orient_axis(A, psi, A.z)

B.dcm(A)


# In[28]:


theta = sm.symbols('theta')

C = me.ReferenceFrame('C')

C.orient_axis(B, theta, B.x)

C.dcm(B)


# In[29]:


phi = sm.symbols('varphi')

D = me.ReferenceFrame('D')

D.orient_axis(C, phi, C.y)

D.dcm(C)


# In[30]:


D.dcm(A)


# In[31]:


A = me.ReferenceFrame('A')
D = me.ReferenceFrame('D')

D.orient_body_fixed(A, (psi, theta, phi), 'zxy')

D.dcm(A)


# In[32]:


psi, theta, phi = sm.symbols('psi, theta, varphi')


# In[33]:


A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
C = me.ReferenceFrame('C')
D = me.ReferenceFrame('D')

B.orient_axis(A, psi, A.z)
C.orient_axis(B, theta, B.x)
D.orient_axis(C, phi, C.z)

D.dcm(A)


# In[34]:


A = me.ReferenceFrame('A')
D = me.ReferenceFrame('D')

D.orient_body_fixed(A, (psi, theta, phi), 'zxz')

D.dcm(A)


# In[35]:


N = me.ReferenceFrame('N')
A = me.ReferenceFrame('A')

q_0, qi, qj, qk = sm.symbols('q_0 q_i q_j q_k')
q = (q_0, qi, qj, qk)
A.orient_quaternion(N, q)
A.dcm(N)


# In[36]:


q = (sm.cos(theta/2), sm.sin(theta/2), 0, 0)
A.orient_quaternion(N, q)
sm.trigsimp(A.dcm(N))

