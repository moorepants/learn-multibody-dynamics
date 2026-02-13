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


# In[4]:


P2.set_pos(P1, la*A.x)
P3.set_pos(P2, lb*B.x)
P4.set_pos(P3, lc*C.x)

P4.pos_from(P1)


# In[5]:


r_P1_P4 = ln*N.x


# In[6]:


loop = P4.pos_from(P1) - r_P1_P4
loop


# In[7]:


fhx = sm.trigsimp(loop.dot(N.x))
fhx


# In[8]:


fhy = sm.trigsimp(loop.dot(N.y))
fhy


# In[9]:


fh = sm.Matrix([fhx, fhy])
fh


# In[10]:


q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')
a, b = sm.symbols('a, b')

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

P2.set_pos(P1, b*A.x)
P3.set_pos(P2, 2*a*B.x)
P4.set_pos(P3, b*C.x)

P4.pos_from(P1)

r_P1_P4 = (2 - sm.S(1)/20)*b*N.x - 2*a*N.y

loop = P4.pos_from(P1) - r_P1_P4

fh_watts = sm.trigsimp(sm.Matrix([loop.dot(N.x), loop.dot(N.y)]))
fh_watts


# In[11]:


import math  # provides pi as a float

repl = {
    la: 1.0,
    lb: 4.0,
    lc: 3.0,
    ln: 5.0,
    q1: 30.0/180.0*math.pi,  # 30 degrees in radians
}
repl


# In[12]:


fh.xreplace(repl)


# In[13]:


q2_guess = -75.0/180.0*math.pi  # -75 degrees in radians
q3_guess = 100.0/180.0*math.pi  # 100 degrees in radians

sol = sm.nsolve(fh.xreplace(repl), (q2, q3), (q2_guess, q3_guess))
sol/math.pi*180.0  # to degrees


# In[14]:


repl = {
    a: 1.0,
    b: 4.0,
    q2: 3.0*math.pi/2.0 - 5.0/180.0*math.pi - q1,
}
repl


# In[15]:


fh_watts.xreplace(repl)


# In[16]:


q1_guess = 10.0/180.0*math.pi
q3_guess = 100.0/180.0*math.pi

sol = sm.nsolve(fh_watts.xreplace(repl), (q1, q3), (q1_guess, q3_guess))
sol/math.pi*180.0  # to degrees


# In[17]:


P1 = me.Point('P1')
P2 = me.Point('P2')
P3 = me.Point('P3')
P4 = me.Point('P4')
P2.set_pos(P1, la*A.x)
P3.set_pos(P2, lb*B.x)
P4.set_pos(P3, lc*C.x)


# In[18]:


P1.set_vel(N, 0)
P4.vel(N)


# In[19]:


sm.trigsimp(P4.vel(N).dot(N.x))


# In[20]:


sm.trigsimp(P4.vel(N).dot(N.y))


# In[21]:


t = me.dynamicsymbols._t
fhd = fh.diff(t)
fhd


# In[22]:


x = sm.Matrix([q2.diff(t), q3.diff(t)])
x


# In[23]:


A = fhd.jacobian(x)
A


# In[24]:


b = -fhd.xreplace({q2.diff(t): 0, q3.diff(t): 0})
b


# In[25]:


x_sol = sm.simplify(A.LUsolve(b))
x_sol


# In[26]:


P4.vel(N).free_dynamicsymbols(N)


# In[27]:


qd_dep_repl = {
  q2.diff(t): x_sol[0, 0],
  q3.diff(t): x_sol[1, 0],
}
qd_dep_repl


# In[28]:


P4.vel(N).xreplace(qd_dep_repl)


# In[29]:


P4.vel(N).xreplace(qd_dep_repl).free_dynamicsymbols(N)

