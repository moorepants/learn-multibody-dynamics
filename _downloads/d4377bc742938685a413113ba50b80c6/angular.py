#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sympy as sm
import sympy.physics.mechanics as me
me.init_vprinting(use_latex='mathjax')


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


cxx, cyy, czz = me.dynamicsymbols('c_{xx}, c_{yy}, c_{zz}')
cxy, cxz, cyx = me.dynamicsymbols('c_{xy}, c_{xz}, c_{yx}')
cyz, czx, czy = me.dynamicsymbols('c_{yz}, c_{zx}, c_{zy}')

B_C_A = sm.Matrix([[cxx, cxy, cxz],
                   [cyx, cyy, cyz],
                   [czx, czy, czz]])


# In[4]:


A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
B.orient_explicit(A, B_C_A.transpose())
B.dcm(A)


# In[5]:


B.x.express(A)


# In[6]:


B.y.express(A)


# In[7]:


B.z.express(A)


# In[8]:


B.y.express(A).dt(A)


# In[9]:


mnx = me.dot(B.y.express(A).dt(A), B.z)
mnx


# In[10]:


mny = me.dot(B.z.express(A).dt(A), B.x)
mny


# In[11]:


mnz = me.dot(B.x.express(A).dt(A), B.y)
mnz


# In[12]:


A_w_B = mnx*B.x + mny*B.y + mnz*B.z
A_w_B


# In[13]:


B_C_A = sm.Matrix([
    [ sm.sqrt(2)/4,  sm.sqrt(2)/2, sm.sqrt(6)/4],
    [-sm.sqrt(3)/2,          0,       sm.S(1)/2],
    [ sm.sqrt(2)/4, -sm.sqrt(2)/2, sm.sqrt(6)/4]
])
B_C_A


# In[14]:


B_C_A_dt = sm.Matrix([
    [-sm.sqrt(6)/2 - 3*sm.sqrt(2)/4, -sm.sqrt(6)/4 + 3*sm.sqrt(2)/2, -3*sm.sqrt(6)/4 + sm.sqrt(2)],
    [                      -1,                     -sm.S(1)/2,               -sm.sqrt(3)],
    [-sm.sqrt(6)/2 + 3*sm.sqrt(2)/4, -sm.sqrt(6)/4 + 3*sm.sqrt(2)/2,            3*sm.sqrt(6)/4]
])
B_C_A_dt


# In[15]:


mnx = (B_C_A[2, :]*B_C_A_dt[1, :].transpose())[0, 0]
mny = (B_C_A[0, :]*B_C_A_dt[2, :].transpose())[0, 0]
mnz = (B_C_A[1, :]*B_C_A_dt[0, :].transpose())[0, 0]

A_w_B = mnx*B.x + mny*B.y + mnz*B.z


# In[16]:


A_w_B.simplify()


# In[17]:


theta = me.dynamicsymbols('theta')

B_C_A = sm.Matrix([[sm.cos(theta), sm.sin(theta), 0],
                   [-sm.sin(theta), sm.cos(theta), 0],
                   [0, 0, 1]])

B_C_A


# In[18]:


A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
B.orient_explicit(A, B_C_A.transpose())

mnx = me.dot(B.y.express(A).dt(A), B.z)
mny = me.dot(B.z.express(A).dt(A), B.x)
mnz = me.dot(B.x.express(A).dt(A), B.y)

A_w_B = mnx*B.x + mny*B.y + mnz*B.z
A_w_B


# In[19]:


A_w_B.simplify()


# In[20]:


A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
B.orient_axis(A, theta, A.z)
B.ang_vel_in(A)


# In[21]:


theta = me.dynamicsymbols('theta')

A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
B.orient_axis(A, theta, A.x + A.y)
B.ang_vel_in(A)


# In[22]:


B.ang_vel_in(A).magnitude()


# In[23]:


psi, theta, phi = me.dynamicsymbols('psi, theta, varphi')

A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
B.orient_body_fixed(A, (psi, theta, phi), 'ZXZ')

mnx = me.dot(B.y.express(A).dt(A), B.z)
mny = me.dot(B.z.express(A).dt(A), B.x)
mnz = me.dot(B.x.express(A).dt(A), B.y)

A_w_B = mnx*B.x + mny*B.y + mnz*B.z
A_w_B.simplify()


# In[24]:


B.ang_vel_in(A)


# In[25]:


psi, theta, phi = me.dynamicsymbols('psi, theta, varphi')

N = me.ReferenceFrame('N')
T = me.ReferenceFrame('T')
T.orient_body_fixed(N, (psi, theta, phi), 'xyz')


# In[26]:


sm.trigsimp(T.dcm(N).xreplace({psi: 0}))


# In[27]:


sm.trigsimp(T.dcm(N).xreplace({psi: sm.pi}))


# In[28]:


sm.trigsimp(T.dcm(N).xreplace({theta: -sm.pi/2}))


# In[29]:


sm.trigsimp(T.dcm(N).xreplace({theta: sm.pi/2}))


# In[30]:


A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
B.orient_body_fixed(A, (psi, theta, 0), 'ZXZ')

u1, u2, u3 = me.dynamicsymbols('u1, u2, u3')

u = u1*B.x + u2*B.y + u3*B.z
u


# In[31]:


u.express(A)


# In[32]:


u.express(A).dt(A)


# In[33]:


u.dt(B)


# In[34]:


A_w_B = B.ang_vel_in(A)
A_w_B


# In[35]:


u.dt(B) + me.cross(A_w_B, u)


# In[36]:


u.express(A).dt(A).express(B).simplify()


# In[37]:


u.dt(A)


# In[38]:


u.dt(B) + me.cross(A_w_B, u)


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


A_w_B = B.ang_vel_in(A)
A_w_B


# In[41]:


B_w_C = C.ang_vel_in(B)
B_w_C


# In[42]:


C_w_D = D.ang_vel_in(C)
C_w_D


# In[43]:


A_w_D = A_w_B + B_w_C + C_w_D
A_w_D


# In[44]:


A2 = me.ReferenceFrame('A')
D2 = me.ReferenceFrame('D')
D2.orient_body_fixed(A2, (psi, theta, phi), 'YXY')
D2.ang_vel_in(A2).simplify()


# In[45]:


A_w_D.express(D)


# In[46]:


theta = me.dynamicsymbols('theta')

A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
B.orient_axis(A, theta, A.z)
B.ang_acc_in(A)


# In[47]:


B.ang_vel_in(A).dt(A)


# In[48]:


B.ang_vel_in(A).dt(B)


# In[49]:


psi, theta, phi = me.dynamicsymbols('psi, theta, varphi')

A = me.ReferenceFrame('A')
D = me.ReferenceFrame('D')
D.orient_body_fixed(A, (psi, theta, phi), 'YXY')

D.ang_acc_in(A).simplify()


# In[50]:


D.ang_vel_in(A).dt(A).simplify()


# In[51]:


D.ang_vel_in(A).dt(D).simplify()


# In[52]:


psi, theta, phi = me.dynamicsymbols('psi, theta, varphi')

A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
C = me.ReferenceFrame('C')
D = me.ReferenceFrame('D')

B.orient_axis(A, psi, A.y)
C.orient_axis(B, theta, B.x)
D.orient_axis(C, phi, C.y)


# In[53]:


A_alp_B = B.ang_acc_in(A)
A_alp_B


# In[54]:


B_alp_C = C.ang_acc_in(B)
B_alp_C


# In[55]:


C_alp_D = D.ang_acc_in(C)
C_alp_D


# In[56]:


A_alp_D = A_alp_B + B_alp_C + C_alp_D
A_alp_D.express(D).simplify()


# In[57]:


D.ang_vel_in(A).dt(A).express(D).simplify()

