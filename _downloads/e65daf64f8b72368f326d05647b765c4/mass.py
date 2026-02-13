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


m = sm.symbols('m')

m_total = m + m + m/2
m_total


# In[4]:


p, R, h = sm.symbols('rho, R, h')  # constants
r, z, theta = sm.symbols('r, z, theta')  # integration variables

sm.integrate(p*r, (r, 0, R/h*z), (theta, 0, 2*sm.pi), (z, 0, h))


# In[5]:


m1, m2, m3 = sm.symbols('m1, m2, m3')
x1, x2, x3 = me.dynamicsymbols('x1, x2, x3')
y1, y2, y3 = me.dynamicsymbols('y1, y2, y3')
z1, z2, z3 = me.dynamicsymbols('z1, z2, z3')

A = me.ReferenceFrame('A')

zeroth_moment = (m1 + m2 + m3)

first_moment = (m1*(x1*A.x + y1*A.y + z1*A.z) +
                m2*(x2*A.x + y2*A.y + z2*A.z) +
                m3*(x3*A.x + y3*A.y + z3*A.z))
first_moment


# In[6]:


r_O_So =  first_moment/zeroth_moment
r_O_So


# In[7]:


r_O_So.xreplace({m2: 2*m1, m3: 3*m1}).simplify()


# In[8]:


m, r, theta = sm.symbols('m, r, theta')
A = me.ReferenceFrame('A')


# In[9]:


r_O_m = (r + r*sm.sin(theta))*A.x + r*sm.cos(theta)*A.y
r_O_2m = (r + r*sm.sin(theta + sm.pi/7))*A.x + r*sm.cos(theta + sm.pi/7)*A.y
r_O_3m = (r + r*sm.sin(theta - sm.pi/6))*A.x + r*sm.cos(theta - sm.pi/6)*A.y


# In[10]:


Iyy = (m*me.dot(r_O_m.cross(A.y), r_O_m.cross(A.y)) +
       2*m*me.dot(r_O_2m.cross(A.y), r_O_2m.cross(A.y)) +
       3*m*me.dot(r_O_3m.cross(A.y), r_O_3m.cross(A.y)))
Iyy


# In[11]:


dIyydtheta = sm.simplify(Iyy.diff(theta))
dIyydtheta


# In[12]:


theta_sol = sm.nsolve((dIyydtheta/m/r**2).evalf(), theta, 0)
theta_sol


# In[13]:


import math

theta_sol*180/math.pi


# In[14]:


kyy = sm.sqrt(Iyy/m)
kyy


# In[15]:


sm.plot(kyy.xreplace({m: 1, r: 1}));


# In[16]:


kyy.xreplace({m: 1, r: 1, theta: theta_sol}).evalf()


# In[17]:


v1, v2, v3 = sm.symbols('v1, v2, v3')
w1, w2, w3 = sm.symbols('w1, w2, w3')

A = me.ReferenceFrame('A')

v = v1*A.x + v2*A.y + v3*A.z
w = w1*A.x + w2*A.y + w3*A.z

Q = me.outer(v, w)
Q


# In[18]:


Q.to_matrix(A)


# In[19]:


me.outer(A.x, A.x)


# In[20]:


me.outer(A.x, A.x).to_matrix(A)


# In[21]:


me.outer(A.y, A.z)


# In[22]:


me.outer(A.y, A.z).to_matrix(A)


# In[23]:


theta = sm.symbols("theta")

A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')

B.orient_axis(A, theta, A.x)

P = 2*me.outer(B.x, B.x) + 3*me.outer(A.x, B.y) + 4*me.outer(B.z, A.z)
P


# In[24]:


P.express(A)


# In[25]:


P.to_matrix(A)


# In[26]:


P.express(B)


# In[27]:


P.to_matrix(B)


# In[28]:


U = me.outer(A.x, A.x) + me.outer(A.y, A.y) + me.outer(A.z, A.z)
U


# In[29]:


U.to_matrix(A)


# In[30]:


U.express(B).simplify()


# In[31]:


Ixx, Iyy, Izz = sm.symbols('I_{xx}, I_{yy}, I_{zz}')
Ixy, Iyz, Ixz = sm.symbols('I_{xy}, I_{yz}, I_{xz}')

I = me.inertia(A, Ixx, Iyy, Izz, ixy=Ixy, iyz=Iyz, izx=Ixz)
I


# In[32]:


I.to_matrix(A)


# In[33]:


I.express(B).simplify()


# In[34]:


I.express(B).simplify().to_matrix(B)


# In[35]:


sm.simplify(B.dcm(A)*I.to_matrix(A)*A.dcm(B))


# In[36]:


N = me.ReferenceFrame('N')

I = (0.25*me.outer(N.x, N.x) +
     0.25*me.outer(N.y, N.y) +
     0.10*me.outer(N.z, N.z) -
     0.07*me.outer(N.x, N.z) -
     0.07*me.outer(N.z, N.x))
I


# In[37]:


H = me.ReferenceFrame('H')
H.orient_axis(N, sm.pi/2 - 68.0*sm.pi/180, N.y)


# In[38]:


I.dot(H.z).dot(H.z).evalf()


# In[39]:


I.to_matrix(N)


# In[40]:


I_H = (H.dcm(N) @ I.to_matrix(N) @ N.dcm(H)).evalf()
I_H


# In[41]:


I_H[2, 2]


# In[42]:


dx, dy, dz, m = sm.symbols('d_x, d_y, d_z, m')

N = me.ReferenceFrame('N')

r_O_Bo = dx*N.x + dy*N.y + dz*N.z

U = me.outer(N.x, N.x) + me.outer(N.y, N.y) + me.outer(N.z, N.z)

I_Bo_O = m*(me.dot(r_O_Bo, r_O_Bo)*U - me.outer(r_O_Bo, r_O_Bo))
I_Bo_O


# In[43]:


I_Bo_O.to_matrix(N)


# In[44]:


I = sm.Matrix([[1.0451, 0.0, -0.1123],
               [0.0, 2.403, 0.0],
               [-0.1123, 0.0, 1.8501]])
I


# In[45]:


ev1, ev2, ev3 = I.eigenvects()


# In[46]:


ev1[0]


# In[47]:


ev1[2][0]


# In[48]:


ev2[0]


# In[49]:


ev2[2][0]


# In[50]:


ev3[0]


# In[51]:


ev3[2][0]


# In[52]:


Ixx, Iyy, Izz = sm.symbols('I_{xx}, I_{yy}, I_{zz}')
Ixy, Iyz, Ixz = sm.symbols('I_{xy}, I_{yz}, I_{xz}')
w1, w2, w3 = me.dynamicsymbols('omega1, omega2, omega3')

B = me.ReferenceFrame('B')

I = me.inertia(B, Ixx, Iyy, Izz, Ixy, Iyz, Ixz)

A_w_B = w1*B.x + w2*B.y + w3*B.z

I.dot(A_w_B)


# In[53]:


I1, I2, I3 = sm.symbols('I_1, I_2, I_3')
w1, w2, w3 = me.dynamicsymbols('omega1, omega2, omega3')

B = me.ReferenceFrame('B')

I = me.inertia(B, I1, I2, I3)

A_w_B = w1*B.x + w2*B.y + w3*B.z

I.dot(A_w_B)

