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


b, t, omega = sm.symbols('b, t, omega')
b, t, omega


# In[6]:


pivot_angle, w2 = sm.symbols('alpha1, omega2')
pivot_angle, w2


# In[7]:


f = sm.Function('f')
f


# In[8]:


type(f)


# In[9]:


f(t)


# In[10]:


type(f(t))


# In[11]:


f(a, b, omega, t)


# In[12]:


expr1 = a + b/omega**2
expr1


# In[13]:


type(expr1)


# In[14]:


sm.srepr(expr1)


# In[15]:


expr2 = f(t) + a*omega
expr2


# In[16]:


expr3 = a*sm.sin(omega) + sm.Abs(f(t))/sm.sqrt(b)
expr3


# In[17]:


expr4 = 5*sm.sin(12) + sm.Abs(-1001)/sm.sqrt(89.2)
expr4


# In[18]:


1/2*a


# In[19]:


sm.S(1)/2*a


# In[20]:


expr5 = t*sm.sin(omega*f(t)) + f(t)/sm.sqrt(t)
expr5


# In[21]:


sm.srepr(expr3)


# In[22]:


repr(expr3)


# In[23]:


print(expr3)


# In[24]:


sm.pprint(expr3)


# In[25]:


sm.latex(expr3)


# In[26]:


print(sm.latex(expr3))


# In[27]:


sm.diff(f(t), t)


# In[28]:


f(t).diff(t)


# In[29]:


expr3


# In[30]:


expr3.diff(b)


# In[31]:


expr3.diff(b, t)


# In[32]:


h = sm.Function('h')
sm.Abs(h(t)).diff(t)


# In[33]:


h = sm.Function('h', real=True)
sm.Abs(h(t)).diff(t)


# In[34]:


h = sm.Function('h', real=True, positive=True)
sm.Abs(h(t)).diff(t)


# In[35]:


expr5


# In[36]:


expr5.diff(t)


# In[37]:


repl = {omega: sm.pi/4, a: 2, f(t): -12, b: 25}


# In[38]:


expr3.xreplace(repl)


# In[39]:


expr3.evalf(n=31, subs=repl)


# In[40]:


type(expr3.evalf(n=31, subs=repl))


# In[41]:


expr3.evalf(n=300, subs=repl)


# In[42]:


float(expr3.evalf(n=300, subs=repl))


# In[43]:


type(float(expr3.evalf(n=300, subs=repl)))


# In[44]:


eval_expr3 = sm.lambdify((omega, a, f(t), b), expr3)


# In[45]:


help(eval_expr3)


# In[46]:


eval_expr3(3.14/4, 2, -12, 25)


# In[47]:


type(eval_expr3(3.14/4, 2, -12, 25))


# In[48]:


mat1 = sm.Matrix([[a, 2*a], [b/omega, f(t)]])
mat1


# In[49]:


mat2 = sm.Matrix([[1, 2], [3, 4]])
mat2


# In[50]:


mat1.shape


# In[51]:


mat1[0, 1]


# In[52]:


mat1[0, 0:2]


# In[53]:


mat1[0:2, 1]


# In[54]:


mat1 + mat2


# In[55]:


mat1*mat2


# In[56]:


mat1@mat2


# In[57]:


sm.hadamard_product(mat1, mat2)


# In[58]:


mat3 = sm.Matrix([expr1, expr2, expr3, expr4, expr5])
mat3


# In[59]:


mat3.diff(a)


# In[60]:


mat3.diff(t)


# In[61]:


mat4 = sm.Matrix([a, b, omega, t])
mat4


# In[62]:


mat3.jacobian(mat4)


# In[63]:


a1, a2 = sm.symbols('a1, a2')

exprs = sm.Matrix([
    [a1*sm.sin(f(t))*sm.cos(2*f(t)) + a2 + omega/sm.log(f(t), t) + 100],
    [a1*omega**2 + f(t)*a2 + omega + f(t)**3],
])
exprs


# In[64]:


A = exprs.jacobian([a1, a2])
A


# In[65]:


b = -exprs.xreplace({a1: 0, a2:0})
b


# In[66]:


A.LUsolve(b)


# In[67]:


sm.simplify(A.LUsolve(b))


# In[68]:


sm.trigsimp(sm.cos(omega)**2 + sm.sin(omega)**2)


# In[69]:


substitutions, simplified = sm.cse(A.LUsolve(b))


# In[70]:


substitutions[0]


# In[71]:


sm.Eq(*substitutions[0])


# In[72]:


sm.Eq(*substitutions[1])


# In[73]:


sm.Eq(*substitutions[2])


# In[74]:


sm.Eq(*substitutions[4])


# In[75]:


simplified[0]

