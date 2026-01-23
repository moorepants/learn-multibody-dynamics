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


x, y, theta = me.dynamicsymbols('x, y, theta')

N = me.ReferenceFrame('N')
A = me.ReferenceFrame('A')

A.orient_axis(N, theta, N.z)

O = me.Point('O')
P = me.Point('P')

P.set_pos(O, x*N.x + y*N.y)

O.set_vel(N, 0)

P.vel(N).express(A)


# In[4]:


fn = P.vel(N).dot(A.y)
fn


# In[5]:


t = me.dynamicsymbols._t

q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')
la, lb, lc, ln = sm.symbols('l_a, l_b, l_c, l_n')

fhx = la*sm.cos(q1) + lb*sm.cos(q1 + q2) + lc*sm.cos(q1 + q2 + q3) - ln
sm.trigsimp(fhx.diff(t))


# In[6]:


dfdx = fn.coeff(x.diff(t))
dfdy = fn.coeff(y.diff(t))
dfdth = fn.coeff(theta.diff(t))

dfdx, dfdy, dfdth


# In[7]:


dfdx.diff(y), dfdy.diff(x)


# In[8]:


dfdx.diff(theta), dfdth.diff(x)


# In[9]:


dfdy.diff(theta), dfdth.diff(y)


# In[10]:


fnx = fhx.diff(t)
dfdq1 = fnx.diff(q1)
dfdq2 = fnx.diff(q2)
dfdq3 = fnx.diff(q3)


# In[11]:


dfdq1.diff(q2) -  dfdq2.diff(q1)


# In[12]:


dfdq2.diff(q3) - dfdq3.diff(q2)


# In[13]:


dfdq3.diff(q1) - dfdq1.diff(q3)


# In[14]:


q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')

A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')

B.orient_body_fixed(A, (q1, q2, q3), 'ZXY')

A_w_B = B.ang_vel_in(A).simplify()
A_w_B


# In[15]:


u1, u2, u3 = me.dynamicsymbols('u1, u2, u3')

t = me.dynamicsymbols._t
qdot = sm.Matrix([q1.diff(t), q2.diff(t), q3.diff(t)])
u = sm.Matrix([u1, u2, u3])

A_w_B = A_w_B.xreplace(dict(zip(qdot, u)))
A_w_B


# In[16]:


Yk_plus_zk = qdot
Yk_plus_zk


# In[17]:


Yk = Yk_plus_zk.jacobian(qdot)
Yk


# In[18]:


qd_zero_repl = dict(zip(qdot, sm.zeros(3, 1)))
qd_zero_repl


# In[19]:


zk = Yk_plus_zk.xreplace(qd_zero_repl)
zk


# In[20]:


sm.Eq(qdot, Yk.LUsolve(u - zk))


# In[21]:


A_w_B = B.ang_vel_in(A).simplify()
A_w_B


# In[22]:


u1_expr = A_w_B.dot(B.x)
u2_expr = A_w_B.dot(B.y)
u3_expr = A_w_B.dot(B.z)

Yk_plus_zk = sm.Matrix([u1_expr, u2_expr, u3_expr])
Yk_plus_zk


# In[23]:


Yk = Yk_plus_zk.jacobian(qdot)
Yk


# In[24]:


zk = Yk_plus_zk.xreplace(qd_zero_repl)
zk


# In[25]:


sm.Eq(qdot, sm.trigsimp(Yk.LUsolve(u - zk)))


# In[26]:


A_w_B = B.ang_vel_in(A).express(A).simplify()
A_w_B


# In[27]:


u1_expr = A_w_B.dot(A.x)
u2_expr = A_w_B.dot(A.y)
u3_expr = A_w_B.dot(A.z)

Yk_plus_zk = sm.Matrix([u1_expr, u2_expr, u3_expr])
Yk_plus_zk


# In[28]:


Yk = Yk_plus_zk.jacobian(qdot)
Yk


# In[29]:


zk = Yk_plus_zk.xreplace(qd_zero_repl)
zk


# In[30]:


sm.Eq(qdot, sm.trigsimp(Yk.LUsolve(u - zk)))


# In[31]:


q1, q2, q3, q4, q5 = me.dynamicsymbols('q1, q2, q3, q4, q5')
l = sm.symbols('l')


# In[32]:


N = me.ReferenceFrame('N')
A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
C = me.ReferenceFrame('C')

A.orient_axis(N, q3, N.z)
B.orient_axis(A, q4, A.z)
C.orient_axis(A, q5, A.z)


# In[33]:


A.ang_vel_in(N)


# In[34]:


B.ang_vel_in(N)


# In[35]:


C.ang_vel_in(N)


# In[36]:


O = me.Point('O')
Ao = me.Point('A_o')
Bo = me.Point('B_o')
Co = me.Point('C_o')

Ao.set_pos(O, q1*N.x + q2*N.y)
Bo.set_pos(Ao, l/2*A.x)
Co.set_pos(Ao, -l/2*A.x)


# In[37]:


O.set_vel(N, 0)
Ao.vel(N)


# In[38]:


Bo.v2pt_theory(Ao, N, A)


# In[39]:


Co.v2pt_theory(Ao, N, A)


# In[40]:


fn = sm.Matrix([Bo.vel(N).dot(B.y),
                Co.vel(N).dot(C.y)])
fn = sm.trigsimp(fn)
fn


# In[41]:


u1, u2, u3, u4, u5 = me.dynamicsymbols('u1, u2, u3, u4, u5')

u_repl = {
    q1.diff(): u1,
    q2.diff(): u2,
    l*q3.diff()/2: u3,
    q4.diff(): u4,
    q5.diff(): u5
}

fn = fn.subs(u_repl)
fn


# In[42]:


us = sm.Matrix([u3, u4, u5])
ur = sm.Matrix([u1, u2])


# In[43]:


Ar = fn.jacobian(ur)
Ar


# In[44]:


As = fn.jacobian(us)
As


# In[45]:


brs = fn.xreplace(dict(zip([u1, u2, u3, u4, u5], [0, 0, 0, 0, 0])))
brs


# In[46]:


An = Ar.LUsolve(-As)
An = sm.simplify(An)
An


# In[47]:


bn = Ar.LUsolve(-brs)
bn


# In[48]:


sm.Eq(ur, An*us + bn)

