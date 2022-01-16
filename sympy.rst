=====
SymPy
=====

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`sympy` or Jupyter notebook:
   :jupyter-download:notebook:`sympy`.

SymPy_ is an open source, collaboratively developed computer algebra system
(CAS) written in Python. We will use it extensively for manipulating symbolic
expressions and equations.  All of the mathematics needed to formulate the
equations of motion of multibody systems can be done with pencil and paper, but
the booking keeping becomes extremely tedious and error prone for systems with
even a small number of bodies. SymPy allows us to let a computer handle the
tedious aspects, e.g. differentiation, solving linear systems of equations, and
reduce the errors one would encounter with pencil and paper. This chapter
introduces SymPy and the primary core SymPy features needed for developing the
equations of motion of multibody systems.

.. _SymPy: https://www.sympy.org

Import and setup
================

We will consistently import SymPy as follows:

.. jupyter-execute::

   import sympy as sm

Since SymPy lets you work with mathematical symbols it's nice to view SymPy
objects in a format that is similar to the math in a textbook. Executing
``init_printing()`` at the beginning of your Jupyter notebook will ensure that
SymPy objects render as typeset mathematics.

.. jupyter-execute::

   sm.init_printing(use_latex='mathjax')

Symbols
=======

Symbols are created with the ``symbols()`` function. A symbol :math:`a` is
created like so:

.. jupyter-execute::

   a = sm.symbols('a')
   a

This symbol object is of the ``Symbol`` type:

.. jupyter-execute::

   type(a)

Multiple symbols can be created with one call to ``symbols()`` and SymPy
recognizes common Greek symbols by spelled out name.

.. jupyter-execute::

   b, t, omega = sm.symbols('b, t, omega')
   b, t, omega

Note that the argument provided to ``symbols()`` does not need to match the
Python variable name it is assigned to. Using more verbose Python variable
names can make code easier to read and understand, especially if there are many
mathematical variables.

.. jupyter-execute::

   pivot_angle, w2 = sm.symbols('alpha1, omega2')
   pivot_angle, w2

Undefined functions
===================

We will also work with undefined mathematical functions in addition to symbols.
These will play an important role in setting up differential equations, where
we typically don't know the function, but only its derivative(s). You can
create arbitrary functions of variables. In this case we make a function of
:math:`t`. First create the function name:

.. jupyter-execute::

   f = sm.Function('f')
   f

This is of a type ``UndefinedFunction``.

.. jupyter-execute::

   type(f)

Now we can create functions of one or more variables like so:

.. jupyter-execute::

   f(t)

Due to SymPy's internal implementations, the type of this object is not defined
as expected:

.. jupyter-execute::

   type(f(t))

The same ``UndefinedFunction`` can be used to create multivariate functions:

.. jupyter-execute::

   f(a, b, omega, t)

Symbolic expressions
====================

Now that we have mathematical variables and functions available, they can be
used to construct mathematical expressions. The most basic way to construct
expressions is with the standard Python operators ``+ - * / **``. For example:

.. jupyter-execute::

   expr1 = a + b/omega**2
   expr1

An expression will have the type ``Add, Mul, or Pow``:

.. jupyter-execute::

   type(expr1)

This is because SymPy stores expressions as a tree_. You can inspect this
internal representation by using the ``srepr()`` function:

.. _tree: https://en.wikipedia.org/wiki/Tree_(graph_theory)

.. jupyter-execute::

   sm.srepr(expr1)

This representation is SymPy's "true" representation of the symbolic
expression. SymPy can display this representation in many other
representations, for example the typeset mathematical expression you have
already seen is one of those representations. This is important to know,
because sometimes the expressions is displayed to you in a way that may be
confusing and checking the ``srepr()`` version can help clear up
misunderstandings. See the `manipulation section`_ of the SymPy tutorial for
more information on this.

.. _manipulation section: https://docs.sympy.org/latest/tutorial/manipulation.html

Undefined functions can also be used to build up expressions:

.. jupyter-execute::

   expr2 = f(t) + a*omega
   expr2

SymPy has a large number of elementary and special functions. See the SymPy
`documentation on functions`_ for more information. For example, here is an
expression that uses ``sin()``, ``Abs()``, and ``sqrt()``:

.. _documentation on functions: https://docs.sympy.org/latest/modules/functions/index.html

.. jupyter-execute::

   expr3 = a*sm.sin(omega) + sm.Abs(f(t))/sm.sqrt(b)
   expr3

Note that Python integers and floats can also be used when constructing
expressions:

.. jupyter-execute::

   expr4 = 5*sm.sin(12) + sm.Abs(-1001)/sm.sqrt(89.2)
   expr4

.. jupyter-execute::

   expr5 = t*sm.sin(omega*f(t)) + f(t)/sm.sqrt(t)
   expr5

Printing
========

We introduced the ``srepr()`` form of SymPy expressions above and mentioned
that expressions can have different representations. For the following
``srepr()`` form:

.. jupyter-execute::

   sm.srepr(expr3)

There is also a standard representation accessed with the ``repr()`` function:

.. jupyter-execute::

   repr(expr3)

This form matches what you typically would type to create the function and it
returns a string. ``print()`` will display that string:

.. jupyter-execute::

   print(expr3)

SymPy also has a "pretty printer" that makes use of unicode symbols to provide
a form that more closely resembles typeset math:

.. jupyter-execute::

   sm.pprint(expr3)

Lastly, the following lines show how SymPy expressions can be represented as
:math:`\LaTeX` code. The double backslashes are present because double
backslashes represent the escape character in Python strings.

.. jupyter-execute::

   sm.latex(expr3)

.. jupyter-execute::

   print(sm.latex(expr3))

.. warning::

   When you are workign with long ezpressions, which will be the case in this
   course, there is no need to print them to the screen. In fact, printing them
   to the screen make take a long time and fill your entire notebook with an
   unreadable mess.

Differentiating
===============

One of the most tedious tasks in formulating equations of motion is the
differentiation of complex trigonometric expressions. SymPy can calculate
derivatives effortlessly. The ``diff()`` SymPy function takes an undefined
function or an expression and differentiates it with respect to the symbol
provided as the second argument:

.. jupyter-execute::

   sm.diff(f(t), t)

All functions and expressions also have a ``.diff()`` method which can be used
like so (many SymPy functions exist as standalone functions and methods):

.. jupyter-execute::

   f(t).diff(t)

``expr3`` is a more complicated expression:

.. jupyter-execute::

   expr3

It can be differentiated, for example, with respect to :math:`b`:

.. jupyter-execute::

   expr3.diff(b)

You can also calculate partial derivatives with respect to successive
variables, as in the following:

.. math::

   \frac{\partial h(a, \omega, t, b)}{\partial t \partial b}

.. jupyter-execute::

   expr3.diff(b, t)

Note that the answer includes real and imaginary components and the signum
function.

.. warning::

   SymPy assumes all symbols are complex valued unless told otherwise. You can
   attach assumptions to symbols to force them to be real, positive, negatives,
   etc. For example, compare these three outputs:

   .. jupyter-execute::

      h = sm.Function('h')
      sm.Abs(h(t)).diff(t)

   .. jupyter-execute::

      h = sm.Function('h', real=True)
      sm.Abs(h(t)).diff(t)

   .. jupyter-execute::

      h = sm.Function('h', real=True, positive=True)
      sm.Abs(h(t)).diff(t)

   Sometimes you may need to add assumptions to variables, but in general it
   will not be necessary.

Lastly, a typical type of derivative you may encounter:

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

This section only scratches the surface of what SymPy can do. The presented
concepts are the basic ones needed for this course, but getting more familiar
with SymPy and what it can do will help. I recommend doing the `SymPy
Tutorial`_. The "Gotchas" section is particularly helpful for common mistakes
when using SymPy. The tutorial is part of the SymPy documentation
https://docs.sympy.org, where you will find general information on SymPy.

.. _SymPy Tutorial: https://docs.sympy.org/latest/tutorial/index.html

If you want to ask a question about using SymPy (or search to see if someone
else has asked your question), you can do so at the following places:

- `SymPy mailing list <https://groups.google.com/g/sympy>`_: Ask questions via
  email.
- `SymPy Gitter <https://gitter.im/sympy/sympy>`_: Ask questions in a live
  chat.
- `Stackoverflow
  <https://stackoverflow.com/questions/tagged/sympy?tab=Votes>`_: Ask and
  search questions on the most popular coding Q&A website.
