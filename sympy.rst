=====
SymPy
=====

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`sympy` or Jupyter notebook:
   :jupyter-download:notebook:`sympy`.

Import and setup
================

This is the standard import we will use in the first portion of the class. It
imports the core SymPy tools and the special tools from the mechanics package.
The init_printing() function is optional but it will increase the number of
outputs that print as symbolic math.

.. jupyter-execute::

   import sympy as sm
   sm.init_printing()

Symbols
=======

.. jupyter-execute::

   a = sm.symbols('a')
   a

.. jupyter-execute::

   type(a)

.. jupyter-execute::

   b, t, omega = sm.symbols('b, t, omega')
   b, t, omega

.. jupyter-execute::

   pivot_angle = sm.symbols('alpha1')
   pivot_angle

Undefined functions
===================

You can create arbitrary functions of variables. In this case we make a
function of t to represent time.

.. jupyter-execute::

   f = sm.Function('f')
   f

.. jupyter-execute::

   type(f)

.. jupyter-execute::

   f(t)

.. jupyter-execute::

   type(f(t))

.. jupyter-execute::

   f(a, b, omega, t)

Symbolic expressions
====================

.. jupyter-execute::

   expr1 = a + b/omega
   expr1

.. jupyter-execute::

   type(expr1)

.. jupyter-execute::

   sm.srepr(expr1)

https://docs.sympy.org/latest/tutorial/manipulation.html

.. jupyter-execute::

   expr2 = f(t) + a*omega
   expr2

https://docs.sympy.org/latest/modules/functions/index.html

.. jupyter-execute::

   expr3 = a*sm.sin(omega) + sm.Abs(f(t))/sm.sqrt(b)
   expr3

.. jupyter-execute::

   expr4 = 5*sm.sin(12) + sm.Abs(-1001)/sm.sqrt(89)
   expr4

.. jupyter-execute::

   expr5 = t*sm.sin(omega*f(t)) + f(t)/sm.sqrt(t)
   expr5

Printing
========

.. jupyter-execute::

   sm.srepr(expr3)

.. jupyter-execute::

   repr(expr3)

.. jupyter-execute::

   print(expr3)

.. jupyter-execute::

   sm.pprint(expr3)

.. jupyter-execute::

   sm.latex(expr3)

.. jupyter-execute::

   print(sm.latex(expr3))

Differentiating
===============

.. jupyter-execute::

   sm.diff(f(t), t)

.. jupyter-execute::

   f(t).diff(t)

.. jupyter-execute::

   expr3

.. jupyter-execute::

   expr3.diff(b)

.. jupyter-execute::

   expr3.diff(b, t)

.. jupyter-execute::

   expr5

.. jupyter-execute::

   expr5.diff(t)

Evaluating symbolic expressions
===============================

.. jupyter-execute::

   expr3.xreplace({omega: sm.pi/4, a: 2, f(t): -12, b: 25})

.. jupyter-execute::

   expr3.evalf(n=31, subs={omega: sm.pi/4, a: 2, f(t): -12, b: 25})

.. jupyter-execute::

   type(expr3.evalf(n=31, subs={omega: sm.pi/4, a: 2, f(t): -12, b: 25}))

.. jupyter-execute::

   eval_expr3 = sm.lambdify((omega, a, f(t), b), expr3)

.. jupyter-execute::

   help(eval_expr3)

.. jupyter-execute::

   eval_expr3(3.14/4, 2, -12, 25)

.. jupyter-execute::

   type(eval_expr3(3.14/4, 2, -12, 25))

Matrices
========

.. jupyter-execute::

   mat1 = sm.Matrix([[a, 2*a], [b/omega, f(t)]])
   mat1

.. jupyter-execute::

   mat2 = sm.Matrix([[1, 2], [3, 4]])
   mat2

.. jupyter-execute::

   mat1.shape

.. jupyter-execute::

   mat1[0, 1]

.. jupyter-execute::

   mat1 + mat2

.. jupyter-execute::

   mat1*mat2

.. jupyter-execute::

   mat3 = sm.Matrix([expr1, expr2, expr3, expr4, expr5])
   mat3

.. jupyter-execute::

   mat1.diff(a)

.. jupyter-execute::

   mat3.diff(t)

.. jupyter-execute::

   mat4 = sm.Matrix([a, b, omega, t])
   mat4

.. jupyter-execute::

   mat3.jacobian(mat4)

Solving Linear Systems
======================

.. jupyter-execute::

   exprs = sm.Matrix([
       [a*sm.sin(f(t))*sm.cos(2*f(t)) + b + omega/sm.log(f(t), t) + 100],
       [a*omega**2 + f(t)*b + omega + f(t)**3],
   ])
   exprs

.. jupyter-execute::

   A = exprs.jacobian([a, b])
   A

.. jupyter-execute::

   b = -exprs.xreplace({a: 0, b:0})
   b

.. jupyter-execute::

   A.LUsolve(b)

Simplification
==============

.. jupyter-execute::

   sm.simplify(A.LUsolve(b))

.. jupyter-execute::

   sm.trigsimp(sm.cos(omega)**2 + sm.sin(omega)**2)

.. jupyter-execute::

   substitutions, simplified = sm.cse(A.LUsolve(b))

.. jupyter-execute::

   substitutions

.. jupyter-execute::

   simplified

Learn more
==========

SymPy Tutorial

https://docs.sympy.org/latest/tutorial/index.html

Where to ask questions about SymPy:

- SymPy mailing list
- SymPy Gitter
- Stackoverflow
