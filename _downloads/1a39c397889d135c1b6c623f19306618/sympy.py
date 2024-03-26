#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sympy as sm


# In[2]:


sm.init_printing(use_latex='mathjax')


# In[3]:


a = sm.symbols('a')
a


# In[4]:


type(a)


# In[5]:


b, t, omega, Omega = sm.symbols('b, t, omega, Omega')
b, t, omega, Omega


# In[6]:


pivot_angle, w2 = sm.symbols('alpha1, omega2')
pivot_angle, w2


# In[7]:


sm.symbols('q1:11')


# In[8]:


f = sm.Function('f')
f


# In[9]:


type(f)


# In[10]:


f(t)


# In[11]:


type(f(t))


# In[12]:


f(a, b, omega, t)


# In[13]:


x, y, z = sm.symbols('x, y, z')
sm.Function('H')(x, y, z)


# In[14]:


expr1 = a + b/omega**2
expr1


# In[15]:


type(expr1)


# In[16]:


sm.srepr(expr1)


# In[17]:


expr2 = f(t) + a*omega
expr2


# In[18]:


expr3 = a*sm.sin(omega) + sm.Abs(f(t))/sm.sqrt(b)
expr3


# In[19]:


expr4 = 5*sm.sin(12) + sm.Abs(-1001)/sm.sqrt(89.2)
expr4


# In[20]:


1/2*a


# In[21]:


sm.S(1)/2*a


# In[22]:


a/2


# In[23]:


expr5 = t*sm.sin(omega*f(t)) + f(t)/sm.sqrt(t)
expr5


# In[24]:


x, s, m = sm.symbols('x, sigma, mu')
sm.exp((x-m)**2/2/s**2)/sm.sqrt(2*sm.pi*s)


# In[25]:


sm.srepr(expr3)


# In[26]:


repr(expr3)


# In[27]:


print(expr3)


# In[28]:


sm.pprint(expr3)


# In[29]:


sm.latex(expr3)


# In[30]:


print(sm.latex(expr3))


# In[31]:


x, s, m = sm.symbols('x, sigma, mu')
print(sm.latex(sm.exp((x-m)**2/2/s**2)/sm.sqrt(2*sm.pi*s),
               mode='equation'))


# In[32]:


sm.diff(f(t), t)


# In[33]:


f(t).diff(t)


# In[34]:


expr3


# In[35]:


expr3.diff(b)


# In[36]:


expr3.diff(b, t)


# In[37]:


h = sm.Function('h')
sm.Abs(h(t)).diff(t)


# In[38]:


h = sm.Function('h', real=True)
sm.Abs(h(t)).diff(t)


# In[39]:


h = sm.Function('h', real=True, positive=True)
sm.Abs(h(t)).diff(t)


# In[40]:


expr5


# In[41]:


expr5.diff(t, omega)


# In[42]:


expr5.diff(t).diff(omega)


# In[43]:


repl = {omega: sm.pi/4, a: 2, f(t): -12, b: 25}


# In[44]:


expr3.xreplace(repl)


# In[45]:


expr3.evalf(n=10, subs=repl)


# In[46]:


type(expr3.evalf(n=10, subs=repl))


# In[47]:


expr3.evalf(n=80, subs=repl)


# In[48]:


float(expr3.evalf(n=300, subs=repl))


# In[49]:


type(float(expr3.evalf(n=300, subs=repl)))


# In[50]:


expr3


# In[51]:


eval_expr3 = sm.lambdify((omega, a, f(t), b), expr3)


# In[52]:


help(eval_expr3)


# In[53]:


eval_expr3(3.14/4, 2, -12, 25)


# In[54]:


type(eval_expr3(3.14/4, 2, -12, 25))


# In[55]:


eval_expr3(3.14/4, 2, -12, [25, 26, 27])


# In[56]:


type(eval_expr3(3.14/4, 2, -12, [25, 26, 27]))


# In[57]:


G, m1, m2, r = sm.symbols('G, m1, m2, r')
F = G*m1*m2/r**2
eval_F = sm.lambdify((G, m1, m2, r), F)
eval_F(6.67430E-11, 5.972E24, 80, 6371E3)


# In[58]:


mat1 = sm.Matrix([[a, 2*a], [b/omega, f(t)]])
mat1


# In[59]:


mat2 = sm.Matrix([[1, 2], [3, 4]])
mat2


# In[60]:


mat1.shape


# In[61]:


mat1[0, 1]


# In[62]:


mat1[0, 0:2]


# In[63]:


mat1[0:2, 1]


# In[64]:


mat1 + mat2


# In[65]:


mat1*mat2


# In[66]:


mat1@mat2


# In[67]:


sm.hadamard_product(mat1, mat2)


# In[68]:


mat3 = sm.Matrix([expr1, expr2, expr3, expr4, expr5])
mat3


# In[69]:


mat3.diff(a)


# In[70]:


mat3.diff(t)


# In[71]:


mat4 = sm.Matrix([a, b, omega, t])
mat4


# In[72]:


mat3.jacobian(mat4)


# In[73]:


def jacobian(v, x):
    """Returns the Jacobian of the vector function v with respect to the
    vector of variables x."""
    diffs = []
    for expr in v:
      for var in x:
         diffs.append(expr.diff(var))
    J_v_x = sm.Matrix(diffs).reshape(len(v), len(x))
    return J_v_x

jacobian(mat3, mat4)


# In[74]:


a1, a2 = sm.symbols('a1, a2')

exprs = sm.Matrix([
    [a1*sm.sin(f(t))*sm.cos(2*f(t)) + a2 + omega/sm.log(f(t), t) + 100],
    [a1*omega**2 + f(t)*a2 + omega + f(t)**3],
])
exprs


# In[75]:


A = exprs.jacobian([a1, a2])
A


# In[76]:


b = -exprs.xreplace({a1: 0, a2: 0})
b


# In[77]:


A.inv() @ b


# In[78]:


A.LUsolve(b)


# In[79]:


L1, L2, L3, L4, L5, L6, F1, F2 = sm.symbols('L1, L2, L3, L4, L5, L6, F1, F2')

exprs = sm.Matrix([
    -L1 + L2 - L3/sm.sqrt(2),
    L3/sm.sqrt(2) + L4 - F1,
    -L2 - L5/sm.sqrt(2),
    L5/sm.sqrt(2) - F2,
    L5/sm.sqrt(2) + L6,
    -L4 -L5/sm.sqrt(2),
])
exprs


# In[80]:


unknowns = sm.Matrix([L1, L2, L3, L4, L5, L6])

coef_mat = exprs.jacobian(unknowns)
rhs = exprs.xreplace(dict(zip(unknowns, [0]*6)))

sol = coef_mat.LUsolve(rhs)

sm.Eq(unknowns, sol)


# In[81]:


eval_sol = sm.lambdify((F1, F2), sol)
eval_sol(13, 32)


# In[82]:


a1, a2 = sm.symbols('a1, a2')
exprs = sm.Matrix([
    [a1*sm.sin(f(t))*sm.cos(2*f(t)) + a2 + omega/sm.log(f(t), t) + 100],
    [a1*omega**2 + f(t)*a2 + omega + f(t)**3],
])
A = exprs.jacobian([a1, a2])
b = -exprs.xreplace({a1: 0, a2: 0})
sol = A.LUsolve(b)


# In[83]:


sm.simplify(sol)


# In[84]:


trig_expr = sm.cos(omega)**2 + sm.sin(omega)**2
trig_expr


# In[85]:


sm.trigsimp(trig_expr)


# In[86]:


sm.count_ops(sol)


# In[87]:


substitutions, simplified = sm.cse(sol)


# In[88]:


substitutions[0]


# In[89]:


sm.Eq(*substitutions[0])


# In[90]:


sm.Eq(*substitutions[1])


# In[91]:


sm.Eq(*substitutions[2])


# In[92]:


sm.Eq(*substitutions[4])


# In[93]:


simplified[0]


# In[94]:


num_ops = sm.count_ops(simplified[0])
for sub in substitutions:
    num_ops += sm.count_ops(sub[1])
num_ops


# In[95]:


a, b, c, x, y, z = sm.symbols('a, b, c, x, y, z')
base_expr = a*sm.sin(x*x + b*sm.cos(x*y) + c*sm.sin(x*z))


# In[96]:


long_expr = base_expr.diff(x, 10)


# In[97]:


eval_long_expr = sm.lambdify((a, b, c, x, y, z), long_expr)
eval_long_expr_cse = sm.lambdify((a, b, c, x, y, z), long_expr, cse=True)


# In[98]:


get_ipython().run_cell_magic('timeit', '', 'eval_long_expr(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)\n')


# In[99]:


get_ipython().run_cell_magic('timeit', '', 'eval_long_expr_cse(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)\n')

