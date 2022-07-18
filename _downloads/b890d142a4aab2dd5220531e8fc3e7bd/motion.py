#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sympy as sm
import sympy.physics.mechanics as me
me.init_vprinting(use_latex='mathjax')


# In[2]:


x, y, theta = me.dynamicsymbols('x, y, theta')

N = me.ReferenceFrame('N')
A = me.ReferenceFrame('A')

A.orient_axis(N, theta, N.z)

O = me.Point('O')
P = me.Point('P')

P.set_pos(O, x*N.x + y*N.y)

O.set_vel(N, 0)

P.vel(N).express(A)


# In[3]:


fn = P.vel(N).dot(A.y)
fn


# In[4]:


t = me.dynamicsymbols._t

q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')
la, lb, lc, ln = sm.symbols('l_a, l_b, l_c, l_n')

fhx = la*sm.cos(q1) + lb*sm.cos(q1 + q2) + lc*sm.cos(q1 + q2 + q3) - ln
sm.trigsimp(fhx.diff(t))


# In[5]:


dfdx = fn.coeff(x.diff(t))
dfdy = fn.coeff(y.diff(t))
dfdth = fn.coeff(theta.diff(t))

dfdx, dfdy, dfdth


# In[6]:


dfdx.diff(y), dfdy.diff(x)


# In[7]:


dfdx.diff(theta), dfdth.diff(x)


# In[8]:


dfdy.diff(theta), dfdth.diff(y)


# In[9]:


q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')

A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')

B.orient_body_fixed(A, (q1, q2, q3), 'ZXY')

A_w_B = B.ang_vel_in(A).simplify()
A_w_B


# In[10]:


u1, u2, u3 = me.dynamicsymbols('u1, u2, u3')

t = me.dynamicsymbols._t
qdot = sm.Matrix([q1.diff(t), q2.diff(t), q3.diff(t)])
u = sm.Matrix([u1, u2, u3])

A_w_B = A_w_B.xreplace(dict(zip(qdot, u)))
A_w_B


# In[11]:


Yk_plus_zk = qdot
Yk_plus_zk


# In[12]:


Yk = Yk_plus_zk.jacobian(qdot)
Yk


# In[13]:


qd_zero_repl = dict(zip(qdot, sm.zeros(3, 1)))
qd_zero_repl


# In[14]:


zk = Yk_plus_zk.xreplace(qd_zero_repl)
zk


# In[15]:


sm.Eq(qdot, Yk.LUsolve(u - zk))


# In[16]:


A_w_B = B.ang_vel_in(A).simplify()
A_w_B


# In[17]:


u1_expr = A_w_B.dot(B.x)
u2_expr = A_w_B.dot(B.y)
u3_expr = A_w_B.dot(B.z)

Yk_plus_zk = sm.Matrix([u1_expr, u2_expr, u3_expr])
Yk_plus_zk


# In[18]:


Yk = Yk_plus_zk.jacobian(qdot)
Yk


# In[19]:


zk = Yk_plus_zk.xreplace(qd_zero_repl)
zk


# In[20]:


sm.Eq(qdot, sm.trigsimp(Yk.LUsolve(u - zk)))


# In[21]:


A_w_B = B.ang_vel_in(A).express(A).simplify()
A_w_B


# In[22]:


u1_expr = A_w_B.dot(A.x)
u2_expr = A_w_B.dot(A.y)
u3_expr = A_w_B.dot(A.z)

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


q1, q2, q3, q4, q5 = me.dynamicsymbols('q1, q2, q3, q4, q5')
l = sm.symbols('l')


# In[27]:


N = me.ReferenceFrame('N')
A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')
C = me.ReferenceFrame('C')

A.orient_axis(N, q3, N.z)
B.orient_axis(A, q4, A.z)
C.orient_axis(A, q5, A.z)


# In[28]:


A.ang_vel_in(N)


# In[29]:


B.ang_vel_in(N)


# In[30]:


C.ang_vel_in(N)


# In[31]:


O = me.Point('O')
Ao = me.Point('A_o')
Bo = me.Point('B_o')
Co = me.Point('C_o')

Ao.set_pos(O, q1*N.x + q2*N.y)
Bo.set_pos(Ao, l/2*A.x)
Co.set_pos(Ao, -l/2*A.x)


# In[32]:


O.set_vel(N, 0)
Ao.vel(N)


# In[33]:


Bo.v2pt_theory(Ao, N, A)


# In[34]:


Co.v2pt_theory(Ao, N, A)


# In[35]:


fn = sm.Matrix([Bo.vel(N).dot(B.y),
                Co.vel(N).dot(C.y)])
fn = sm.trigsimp(fn)
fn


# In[36]:


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


# In[37]:


us = sm.Matrix([u3, u4, u5])
ur = sm.Matrix([u1, u2])


# In[38]:


Ar = fn.jacobian(ur)
Ar


# In[39]:


As = -fn.jacobian(us)
As


# In[40]:


bs = -fn.xreplace(dict(zip([u1, u2, u3, u4, u5], [0, 0, 0, 0, 0])))
bs


# In[41]:


An = Ar.LUsolve(As)
An = sm.simplify(An)
An


# In[42]:


bn = Ar.LUsolve(bs)
bn


# In[43]:


sm.Eq(ur, An*us + bn)

