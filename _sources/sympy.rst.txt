=====
SymPy
=====

.. note::

   You can download this example as a Python script:
   :jupyter-download-script:`sympy` or Jupyter notebook:
   :jupyter-download-notebook:`sympy`.

Learning Objectives
===================

After completing this chapter readers will be able to:

- Write mathematical expressions with symbols and functions using SymPy.
- Print different forms of expressions and equations with SymPy.
- Differentiate mathematical expressions using SymPy.
- Evaluate mathematical expressions using SymPy.
- Create matrices and do linear algebra using SymPy.
- Solve a linear system of equations with SymPy.
- Simplify mathematical expressions with SymPy.

Introduction
============

SymPy_ is an open source, collaboratively developed `computer algebra system`_
(CAS) written in Python. It will be used extensively for manipulating symbolic
expressions and equations. All of the mathematics needed to formulate the
equations of motion of multibody systems can be done with pencil and paper, but
the bookkeeping becomes extremely tedious and error prone for systems with even
a small number of bodies. SymPy lets a computer handle the tedious aspects
(e.g. differentiation or solving linear systems of equations) and reduces the
errors one would encounter with pencil and paper. This chapter introduces SymPy
and the primary SymPy features we will be using.

.. _SymPy: https://www.sympy.org
.. _computer algebra system: https://en.wikipedia.org/wiki/Computer_algebra_system

Import and Setup
================

I will import SymPy as follows throughout this book:

.. jupyter-execute::

   import sympy as sm

Since SymPy works with mathematical symbols it's nice to view SymPy objects in
a format that is similar to the math in a textbook. Executing
:external:py:func:`~sympy.interactive.printing.init_printing` at the beginning
of your Jupyter Notebook will ensure that SymPy objects render as typeset
mathematics. I use the ``use_latex='mathjax'`` argument here to disable math
png image generation, but that keyword argument is not necessary.

.. jupyter-execute::

   sm.init_printing(use_latex='mathjax')

Symbols
=======

Mathematical symbols are created with the
:external:py:func:`~sympy.core.symbol.symbols` function. A symbol :math:`a` is
created like so:

.. jupyter-execute::

   a = sm.symbols('a')
   a

This symbol object is of the :external:py:class:`~sympy.core.symbol.Symbol` type:

.. jupyter-execute::

   type(a)

Multiple symbols can be created with one call to ``symbols()`` and SymPy
recognizes common Greek symbols by their spelled-out name.

.. jupyter-execute::

   b, t, omega, Omega = sm.symbols('b, t, omega, Omega')
   b, t, omega, Omega

Note that the argument provided to ``symbols()`` does not need to match the
Python variable name it is assigned to. Using more verbose Python variable
names may make code easier to read and understand, especially if there are many
mathematical variables that you need to keep track of. Note that the subscripts
are recognized too.

.. jupyter-execute::

   pivot_angle, w2 = sm.symbols('alpha1, omega2')
   pivot_angle, w2

.. admonition:: Exercise

   Review the SymPy documentation and create symbols :math:`q_1, q_2, \ldots,
   q_{10}` with a very succint call to
   :external:py:func:`~sympy.core.symbol.symbols`.

.. admonition:: Solution
   :class: dropdown

   .. jupyter-execute::

      sm.symbols('q1:11')

Undefined Functions
===================

You will also work with undefined mathematical functions in addition to symbols.
These will play an important role in setting up differential equations, where
you typically don't know the function, but only its derivative(s). You can
create arbitrary functions of variables. In this case, you make a function of
:math:`t`. First create the function name:

.. jupyter-execute::

   f = sm.Function('f')
   f

This is of a type ``sympy.core.function.UndefinedFunction``.

.. jupyter-execute::

   type(f)

Now you can create functions of one or more variables like so:

.. jupyter-execute::

   f(t)

.. warning::

   Due to SymPy's internal implementations, the type of a function with its
   argument is not defined as expected:

   .. jupyter-execute::

      type(f(t))

   This can be confusing if you are checking types.

The same ``UndefinedFunction`` can be used to create multivariate functions:

.. jupyter-execute::

   f(a, b, omega, t)

.. admonition:: Exercise

   Create a function :math:`H(x, y, z)`.


.. admonition:: Solution
   :class: dropdown

   .. jupyter-execute::

      x, y, z = sm.symbols('x, y, z')
      sm.Function('H')(x, y, z)

Symbolic Expressions
====================

Now that you have mathematical variables and functions available, they can be
used to construct mathematical expressions. The most basic way to construct
expressions is with the standard Python operators ``+``, ``-``, ``*``, ``/``,
and ``**``. For example:

.. jupyter-execute::

   expr1 = a + b/omega**2
   expr1

An expression will have the type ``Add, Mul, or Pow``:

.. jupyter-execute::

   type(expr1)

This is because SymPy stores expressions behind the scenes as a tree_. You can
inspect this internal representation by using the
:external:py:func:`~sympy.printing.repr.srepr` function:

.. _tree: https://en.wikipedia.org/wiki/Tree_(graph_theory)

.. jupyter-execute::

   sm.srepr(expr1)

This is a visual representation of the tree:

.. Use ``print(sm.dotprint(expr1))`` to get the following code.

.. graphviz::
   :align: center

   digraph{

   # Graph style
   "ordering"="out"
   "rankdir"="TD"

   #########
   # Nodes #
   #########

   "Add(Symbol('a'), Mul(Symbol('b'), Pow(Symbol('omega'), Integer(-2))))_()" ["color"="black", "label"="Add", "shape"="ellipse"];
   "Symbol('a')_(0,)" ["color"="black", "label"="a", "shape"="ellipse"];
   "Mul(Symbol('b'), Pow(Symbol('omega'), Integer(-2)))_(1,)" ["color"="black", "label"="Mul", "shape"="ellipse"];
   "Symbol('b')_(1, 0)" ["color"="black", "label"="b", "shape"="ellipse"];
   "Pow(Symbol('omega'), Integer(-2))_(1, 1)" ["color"="black", "label"="Pow", "shape"="ellipse"];
   "Symbol('omega')_(1, 1, 0)" ["color"="black", "label"="omega", "shape"="ellipse"];
   "Integer(-2)_(1, 1, 1)" ["color"="black", "label"="-2", "shape"="ellipse"];

   #########
   # Edges #
   #########

   "Add(Symbol('a'), Mul(Symbol('b'), Pow(Symbol('omega'), Integer(-2))))_()" -> "Symbol('a')_(0,)";
   "Add(Symbol('a'), Mul(Symbol('b'), Pow(Symbol('omega'), Integer(-2))))_()" -> "Mul(Symbol('b'), Pow(Symbol('omega'), Integer(-2)))_(1,)";
   "Mul(Symbol('b'), Pow(Symbol('omega'), Integer(-2)))_(1,)" -> "Symbol('b')_(1, 0)";
   "Mul(Symbol('b'), Pow(Symbol('omega'), Integer(-2)))_(1,)" -> "Pow(Symbol('omega'), Integer(-2))_(1, 1)";
   "Pow(Symbol('omega'), Integer(-2))_(1, 1)" -> "Symbol('omega')_(1, 1, 0)";
   "Pow(Symbol('omega'), Integer(-2))_(1, 1)" -> "Integer(-2)_(1, 1, 1)";
   }

This representation is SymPy's "true" representation of the symbolic
expression. SymPy can display this expression in many other representations,
for example the typeset mathematical expression you have already seen is one of
those representations. This is important to know, because sometimes the
expressions are displayed to you in a way that may be confusing and checking
the ``srepr()`` version can help clear up misunderstandings. See the
`manipulation section`_ of the SymPy tutorial for more information on this.

.. _manipulation section: https://docs.sympy.org/latest/tutorial/manipulation.html

Undefined functions can also be used in expressions just like symbols:

.. jupyter-execute::

   expr2 = f(t) + a*omega
   expr2

SymPy has a large number of elementary and special functions. See the SymPy
`documentation on functions`_ for more information. For example, here is an
expression that uses
:external:py:class:`~sympy.functions.elementary.trigonometric.sin`,
:external:py:class:`~sympy.functions.elementary.complexes.Abs`, and
:external:py:func:`~sympy.functions.elementary.miscellaneous.sqrt`:

.. _documentation on functions: https://docs.sympy.org/latest/modules/functions/index.html

.. jupyter-execute::

   expr3 = a*sm.sin(omega) + sm.Abs(f(t))/sm.sqrt(b)
   expr3

Note that Python integers and floats can also be used when constructing
expressions:

.. jupyter-execute::

   expr4 = 5*sm.sin(12) + sm.Abs(-1001)/sm.sqrt(89.2)
   expr4

.. warning::

   Be careful with numbers, as SymPy may not intepret them as expected. For
   example:

   .. jupyter-execute::

      1/2*a

   Python does the division before it is multiplied by ``a``, thus a floating
   point value is created. To fix this you can use the ``S()`` function to
   "sympify" numbers:

   .. jupyter-execute::

      sm.S(1)/2*a

   Or you can ensure the symbol comes first in the division operation:

   .. jupyter-execute::

      a/2

Lastly, an expression of ``t``:

.. jupyter-execute::

   expr5 = t*sm.sin(omega*f(t)) + f(t)/sm.sqrt(t)
   expr5

.. admonition:: Exercise

   Create an expression for the normal distribution function:

   .. math::

      \frac{1}{\sqrt{2\pi\sigma}}e^{\frac{(x-\mu)^2}{2\sigma^2}}

.. admonition:: Solution
   :class: dropdown

   .. jupyter-execute::

      x, s, m = sm.symbols('x, sigma, mu')
      sm.exp((x-m)**2/2/s**2)/sm.sqrt(2*sm.pi*s)

   Notice that SymPy does some minor manipulation of the expression, but it is
   equivalent to the form shown in the prompt.

Printing
========

I introduced the ``srepr()`` form of SymPy expressions above and mentioned that
expressions can have different representations. For the following ``srepr()``
form:

.. jupyter-execute::

   sm.srepr(expr3)

There is also a standard representation accessed with the ``repr()`` function:

.. jupyter-execute::

   repr(expr3)

This form matches what you typically would type to create the expression and it
returns a string. The ``print()`` function will display that string:

.. jupyter-execute::

   print(expr3)

SymPy also has a "pretty printer" (:external:py:func:`pprint()
<sympy.printing.pretty.pretty.pretty_print>`) that makes use of unicode symbols
to provide a form that more closely resembles typeset math:

.. jupyter-execute::

   sm.pprint(expr3)

Lastly, the following lines show how SymPy expressions can be represented as
LaTeX code using :external:py:func:`sympy.printing.latex.latex`. The double
backslashes are present because double backslashes represent the escape
character in Python strings.

.. jupyter-execute::

   sm.latex(expr3)

.. jupyter-execute::

   print(sm.latex(expr3))

.. warning::

   When you are working with long expressions, which will be the case in this
   course, there is no need to print them to the screen. In fact, printing them
   to the screen make take a long time and fill your entire notebook with an
   unreadable mess.

.. admonition:: Exercise

   Print the normal distribution expression

   .. math::

      \frac{1}{\sqrt{2\pi\sigma}}e^{\frac{(x-\mu)^2}{2\sigma^2}}

   as a LaTeX string inside an equation environment.

.. admonition:: Solution
   :class: dropdown

   .. jupyter-execute::

      x, s, m = sm.symbols('x, sigma, mu')
      print(sm.latex(sm.exp((x-m)**2/2/s**2)/sm.sqrt(2*sm.pi*s),
                     mode='equation'))

Differentiating
===============

One of the most tedious tasks in formulating equations of motion is
differentiating complex trigonometric expressions. SymPy can calculate
derivatives effortlessly. The :external:py:func:`~sympy.core.function.diff`
SymPy function takes an undefined function or an expression and differentiates
it with respect to the symbol provided as the second argument:

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
variables. If you want to first differentiate with respect to :math:`b` and
then with respect to :math:`t` as in the following operation:

.. math::

   \frac{\partial^2 h(a, \omega, t, b)}{\partial t \partial b}

where:

.. math::

   h(a, \omega, t, b) = \displaystyle a \sin{\left(\omega \right)} + \frac{\left|{f{\left(t \right)}}\right|}{\sqrt{b}}

then you can use successive arguments to ``.diff()``:

.. jupyter-execute::

   expr3.diff(b, t)

Note that the answer includes real and imaginary components and the `signum
function`_.

.. _signum function: https://en.wikipedia.org/wiki/Sign_function

.. warning::

   SymPy assumes all symbols are complex-valued unless told otherwise. You can
   attach assumptions to symbols to force them to be real, positive, negative,
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
   will not be necessary. Read more about `assumptions in SymPy's guide
   <https://docs.sympy.org/latest/guides/assumptions.html>`_.

.. admonition:: Exercise

   Differentiate ``expr5`` above using this operator:

   .. math::

      \frac{\partial^2}{\partial \omega \partial t}

.. admonition:: Solution
   :class: dropdown

   First show ``expr5``:

   .. jupyter-execute::

      expr5

   The twice partial derivative is:

   .. jupyter-execute::

      expr5.diff(t, omega)

   or you can chain ``.diff()`` calls:

   .. jupyter-execute::

      expr5.diff(t).diff(omega)

Evaluating Symbolic Expressions
===============================

SymPy expressions can be evaluated numerically in several ways. The
:external:py:meth:`~sympy.core.basic.Basic.xreplace` method allows substitution
of exact symbols or sub-expressions. First create a dictionary that maps
symbols, functions or sub-expressions to the replacements:

.. jupyter-execute::

   repl = {omega: sm.pi/4, a: 2, f(t): -12, b: 25}

This dictionary can then be passed to ``.xreplace()``:

.. jupyter-execute::

   expr3.xreplace(repl)

Notice how the square root and fraction do not automatically reduce to their
decimal equivalents. To do so, you must use the
:external:py:meth:`~sympy.core.evalf.EvalfMixin.evalf` method. This method will
evaluate an expression to an arbitrary number of decimal points.  You provide
the number of decimal places and the substitution dictionary to evaluate:

.. jupyter-execute::

   expr3.evalf(n=10, subs=repl)

.. jupyter-execute::

   type(expr3.evalf(n=10, subs=repl))

Note that this is a SymPy :external:py:class:`~sympy.core.numbers.Float`
object, which is a special object that can have an arbitrary number of decimal
places, for example here is the expression evaluated to 80 decimal places:

.. jupyter-execute::

   expr3.evalf(n=80, subs=repl)

To convert this to Python floating point number, use ``float()``:

.. jupyter-execute::

   float(expr3.evalf(n=300, subs=repl))

.. jupyter-execute::

   type(float(expr3.evalf(n=300, subs=repl)))

This value is a `machine precision`_ floating point value and can be used with
standard Python functions that operate on floating point numbers.

.. _machine precision: https://en.wikipedia.org/wiki/Machine_epsilon

To obtain machine precision floating point numbers directly and with more
flexibility, it is better to use the
:external:py:func:`~sympy.utilities.lambdify.lambdify` function to convert the
expression to a Python function. When using ``lambdify()``, all symbols and
functions should be converted to numbers so first identify what symbols and
functions make up the expression.

.. jupyter-execute::

   expr3

:math:`\omega, a, f(t)`, and :math:`b` are all present in the expression. The
first argument of ``lambdify()`` should be a sequence of all these symbols and
functions and the second argument should be the expression.

.. jupyter-execute::

   eval_expr3 = sm.lambdify((omega, a, f(t), b), expr3)

``lambdify()`` generates a Python function and, in this case, we store that
function in the variable ``eval_expr3``. You can see what the inputs and
outputs of the function are with ``help()``:

.. jupyter-execute::

   help(eval_expr3)

This function operates on and returns floating point values, for example:

.. jupyter-execute::

   eval_expr3(3.14/4, 2, -12, 25)

The type of lambdify's return values will be NumPy_ floats.

.. jupyter-execute::

   type(eval_expr3(3.14/4, 2, -12, 25))

.. _NumPy: https://www.numpy.org

These floats are interoperable with Python floats for single values (unlike
SymPy Floats) but also support arrays of floats. For example:

.. jupyter-execute::

   eval_expr3(3.14/4, 2, -12, [25, 26, 27])

.. jupyter-execute::

   type(eval_expr3(3.14/4, 2, -12, [25, 26, 27]))

More on NumPy arrays of floats will be introduced in a later chapter.

.. warning:: Python and NumPy floats can be mixed, but avoid mixing SymPy
   Floats with either.

.. note::

   This distinction between SymPy ``Float`` objects and regular Python and
   NumPy ``float`` objects is important. In this case, the Python float and the
   NumPy float are equivalent. The later will compute much faster because
   arbitrary precision is not required. In this book, you will almost always
   want to convert SymPy expressions to machine precision floating point
   numbers, so use ``lambdify()``.

.. admonition:: Exercise

   Create a symbolic expression representing `Newton's Law of Universal
   Gravitation
   <https://en.wikipedia.org/wiki/Newton's_law_of_universal_gravitation>`_. Use
   ``lambdify()`` to evaluate the expression for two mass of 5.972E24 kg and 80
   kg at a distance of 6371 km apart to find the gravitational force in
   Newtons.

.. admonition:: Solution
   :class: dropdown

   .. jupyter-execute::

      G, m1, m2, r = sm.symbols('G, m1, m2, r')
      F = G*m1*m2/r**2
      eval_F = sm.lambdify((G, m1, m2, r), F)
      eval_F(6.67430E-11, 5.972E24, 80, 6371E3)

Matrices
========

SymPy supports matrices of expressions and linear algebra. Many of the
operations needed in multibody dynamics are more succinctly formulated with
matrices and linear algebra. Matrices can be created by passing nested lists to
the :external:py:class:`Matrix() <sympy.matrices.dense.MutableDenseMatrix>`
object. For example:

.. jupyter-execute::

   mat1 = sm.Matrix([[a, 2*a], [b/omega, f(t)]])
   mat1

.. jupyter-execute::

   mat2 = sm.Matrix([[1, 2], [3, 4]])
   mat2

All matrices are two dimensional and the number of rows and columns, in that
order, are stored in the ``.shape`` attribute.

.. jupyter-execute::

   mat1.shape

Individual elements of the matrix can be extracted with the bracket notation
taking the row and column indices (remember Python indexes from 0):

.. jupyter-execute::

   mat1[0, 1]

The slice notation can extract rows or columns:

.. jupyter-execute::

   mat1[0, 0:2]

.. jupyter-execute::

   mat1[0:2, 1]

Matrix algebra can be performed. Matrices can be added:

.. jupyter-execute::

   mat1 + mat2

Both the ``*`` and the ``@`` operator perform matrix multiplication:

.. jupyter-execute::

   mat1*mat2

.. jupyter-execute::

   mat1@mat2

Element-by-element multiplication requires the ``sympy.hadamard_product()``
function:

.. jupyter-execute::

   sm.hadamard_product(mat1, mat2)

Note that NumPy uses ``*`` for element-by-element multiplication and ``@`` for matrix multiplication,
so to avoid possible confusion, use ``@`` for SymPy matrix multiplication.

Differentiation operates on each element of the matrix:

.. jupyter-execute::

   mat3 = sm.Matrix([expr1, expr2, expr3, expr4, expr5])
   mat3

.. jupyter-execute::

   mat3.diff(a)

.. jupyter-execute::

   mat3.diff(t)

If you have column vectors :math:`\bar{v}` and :math:`\bar{u}`, the
:math:`(i,j)` entries of the Jacobian of :math:`\bar{v}` with respect to the
entries in vector :math:`\bar{u}` are found with :math:`\mathbf{J}_{ij} =
\frac{\partial v_i}{\partial u_j}`.  The Jacobian_ matrix of vector (column
matrix) can be formed with the
:external:py:meth:`~sympy.matrices.matrices.MatrixCalculus.jacobian` method.
This calculates the partial derivatives of each element in the vector with
respect to a vector (or sequence) of variables.

.. jupyter-execute::

   mat4 = sm.Matrix([a, b, omega, t])
   mat4

.. jupyter-execute::

   mat3.jacobian(mat4)

.. _Jacobian: https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant

.. admonition:: Exercise

   Write your own function that produces a Jacobian given a column matrix of
   expressions. It should look like::

      def jacobian(v, x):
          """Returns the Jacobian of the vector function v with respect to the
          vector of variables x."""
          # fill in your code here
          return J_v_x

   Show that it gives the same solution as the above ``.jacobian()`` method. Do
   not use the ``.jacobian()`` method in your function.

.. admonition:: Solution
   :class: dropdown

   .. jupyter-execute::

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

Solving Linear Systems
======================

You'll need to solve `linear systems of equations`_ often in this book. SymPy
offers `a number of ways to do this`_, but the best way to do so if you know a
set of equations are linear in specific variables is the method described
below. First, you should confirm you have equations of this form:

.. math::

   a_{11} x_1 + a_{12} x_2 + \ldots + a_{1n} x_n + b_1 = 0 \\
   a_{21} x_1 + a_{22} x_2 + \ldots + a_{2n} x_n + b_2 = 0 \\
   \vdots \\
   a_{n1} x_1 + a_{n2} x_2 + \ldots + a_{nn} x_n + b_n = 0

.. _linear systems of equations: https://en.wikipedia.org/wiki/System_of_linear_equations
.. _a number of ways to do this: https://docs.sympy.org/latest/guides/solving/index.html

These equations can be put into matrix form:

.. math::

   \mathbf{A}\bar{x} = \bar{b}

where:

.. math::

   \mathbf{A} =
   \begin{bmatrix}
     a_{11} & a_{12} & \ldots & a_{1n} \\
     a_{21} & a_{22} & \ldots & a_{2n} \\
     \ldots & \ldots & \ldots & \ldots \\
     a_{n1} & a_{n2} & \ldots & a_{nn}
   \end{bmatrix},
   \bar{x} =
   \begin{bmatrix}
     x_1 \\
     x_2 \\
     \ldots \\
     x_n
   \end{bmatrix},
   \bar{b} =
   \begin{bmatrix}
     -b_1 \\
     -b_2 \\
     \ldots \\
     -b_n
   \end{bmatrix}

:math:`\bar{x}`, the solution, is found with matrix inversion (if the matrix is
invertible):

.. math::

   \bar{x} = \mathbf{A}^{-1}\bar{b}

Taking the inverse is not computationally efficient and potentially numerically
inaccurate, so some form of `Gaussian elmination`_ should be used to solve the
system.

.. _Gaussian elmination: https://en.wikipedia.org/wiki/Gaussian_elimination

To solve with SymPy, start with a column matrix of linear expressions:

.. jupyter-execute::

   a1, a2 = sm.symbols('a1, a2')

   exprs = sm.Matrix([
       [a1*sm.sin(f(t))*sm.cos(2*f(t)) + a2 + omega/sm.log(f(t), t) + 100],
       [a1*omega**2 + f(t)*a2 + omega + f(t)**3],
   ])
   exprs

Since we know these two expressions are linear in the :math:`a_1` and
:math:`a_2` variables, the partial derivatives with respect to those two
variables will return the linear coefficients. The :math:`\mathbf{A}` matrix
can be formed in one step with the ``.jacobian()`` method:

.. jupyter-execute::

   A = exprs.jacobian([a1, a2])
   A

The :math:`\bar{b}` vector can be formed by setting :math:`a_1=a_2=0`, leaving the
terms that are not linear in :math:`a_1` and :math:`a_2`.

.. jupyter-execute::

   b = -exprs.xreplace({a1: 0, a2: 0})
   b

The :external:py:meth:`~sympy.matrices.matrices.MatrixBase.inv` method can
compute the inverse of A to find the solution:

.. jupyter-execute::

   A.inv() @ b

But it is best to use the
:external:py:meth:`~sympy.matrices.matrices.MatrixBase.LUsolve` method to
perform an `LU decomposition`_ Gaussian-Elimination to solve the system,
especially as the dimension of :math:`\mathbf{A}` grows:

.. jupyter-execute::

   A.LUsolve(b)

.. _LU decomposition: https://en.wikipedia.org/wiki/LU_decomposition

.. warning::

   This method of solving symbolic linear systems is fast, but it can give
   incorrect answers for:

   1. expressions that are not acutally linear in the variables the Jacobian is
      taken with respect to
   2. :math:`\mathbf{A}` matrix entries that would evaluate to zero if
      simplified or specific numerical values are provided

   So only use this method if you are sure your equations are linear and if
   your :math:`\mathbf{A}` matrix is made up of complex expressions, watch out
   for ``nan`` results after lambdifying.
   :external:py:func:`~sympy.solvers.solvers.solve` and
   :external:py:func:`~sympy.solvers.solveset.linsolve` can also solve linear
   systems and they check for linearity and properties of the A matrix.  The
   cost is that they can be extremely slow for large expressions (which we will
   have in this book).

.. admonition:: Exercise

   Solve the following equations for all of the :math:`L`'s and then use
   ``lambdify()`` to evaluate the solution for :math:`F_1=13` and
   :math:`F_2=32`.

   .. math::

      -L_1 + L_2 - L_3/\sqrt{2} = & 0 \\
      L_3/\sqrt{2} + L_4 = &  F_1 \\
      -L_2 - L_5/\sqrt{2} = &  0 \\
      L_5/\sqrt{2} = & F_2 \\
      L_5/\sqrt{2} + L_6 = &  0 \\
      -L_4 -L_5/\sqrt{2} = &  0

.. admonition:: Solution
   :class: dropdown

   .. jupyter-execute::

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

   .. jupyter-execute::

      unknowns = sm.Matrix([L1, L2, L3, L4, L5, L6])

      coef_mat = exprs.jacobian(unknowns)
      rhs = exprs.xreplace(dict(zip(unknowns, [0]*6)))

      sol = coef_mat.LUsolve(rhs)

      sm.Eq(unknowns, sol)

   .. jupyter-execute::

      eval_sol = sm.lambdify((F1, F2), sol)
      eval_sol(13, 32)

Simplification
==============

The above result from
:external:py:meth:`~sympy.matrices.matrices.MatrixBase.LUsolve` is a bit
complicated. Reproduced here:

.. jupyter-execute::

   a1, a2 = sm.symbols('a1, a2')
   exprs = sm.Matrix([
       [a1*sm.sin(f(t))*sm.cos(2*f(t)) + a2 + omega/sm.log(f(t), t) + 100],
       [a1*omega**2 + f(t)*a2 + omega + f(t)**3],
   ])
   A = exprs.jacobian([a1, a2])
   b = -exprs.xreplace({a1: 0, a2: 0})
   sol = A.LUsolve(b)

SymPy has some functionality for automatically simplifying symbolic
expressions. The function :external:py:func:`~sympy.simplify.simplify.simplify`
will attempt to find a simpler version:

.. jupyter-execute::

   sm.simplify(sol)

But you'll have the best luck at simplifying if you use simplification functions that
target the type of expression you have. The
:external:py:func:`~sympy.simplify.trigsimp.trigsimp` function only attempts
trigonometric simplifications, for example:

.. jupyter-execute::

   trig_expr = sm.cos(omega)**2 + sm.sin(omega)**2
   trig_expr

.. jupyter-execute::

   sm.trigsimp(trig_expr)

.. warning::

   Only attempt simplification on expressions that are several lines of text.
   Larger expressions become increasingly computationally intensive to simplify
   and there is generally no need to do so.

As mentioned earlier, SymPy represents expressions as trees. Symbolic
expressions can also be represented as `directed acyclic graphs`_ that contain
only one node for each unique expression (unlike SymPy's trees which may have
the same expression in more than one node). These unique expressions, or
"common subexpressions", can be found with the
:external:py:func:`~sympy.simplify.cse_main.cse` function. This function will
provide a simpler form of the equations that minimizes the number of operations
to compute the answer. We can count the number of basic operations (additions,
multiplies, etc.) using :external:py:func:`~sympy.core.function.count_ops`:

.. jupyter-execute::

   sm.count_ops(sol)

.. _Directed acyclic graphs: https://en.wikipedia.org/wiki/Directed_acyclic_graph

We can simplify with ``cse()``:

.. jupyter-execute::

   substitutions, simplified = sm.cse(sol)

The ``substitutions`` variable contains a list of tuples, where each tuple has
a new intermediate variable and the sub-expression it is equal to.

.. jupyter-execute::

   substitutions[0]

The :external:py:class:`Eq() <sympy.core.relational.Equality>` class with tuple
unpacking (``*``) can be used to display these tuples as equations:

.. jupyter-execute::

   sm.Eq(*substitutions[0])

.. jupyter-execute::

   sm.Eq(*substitutions[1])

.. jupyter-execute::

   sm.Eq(*substitutions[2])

.. jupyter-execute::

   sm.Eq(*substitutions[4])

The ``simplified`` variable contains the simplified expression, made up of the
intermediate variables.

.. jupyter-execute::

   simplified[0]

We can count the number of operations of the simplified version:

.. jupyter-execute::

   num_ops = sm.count_ops(simplified[0])
   for sub in substitutions:
       num_ops += sm.count_ops(sub[1])
   num_ops

.. admonition:: Exercise

   :external:py:func:`~sympy.utilities.lambdify.lambdify` has an optional
   argument ``cse=True|False`` that applies common subexpression elimination
   internally to simplify the number of operations. Differentiate the
   ``base_expr`` with respect to ``x`` 10 times to generate a very long
   expression. Create two functions using ``lambdify()``, one with ``cse=True``
   and one with ``cse=False``. Compare how long it takes to numerically
   evaluate the resulting functions using the ``%timeit`` magic.

   .. jupyter-execute::

      a, b, c, x, y, z = sm.symbols('a, b, c, x, y, z')
      base_expr = a*sm.sin(x*x + b*sm.cos(x*y) + c*sm.sin(x*z))

.. admonition:: Solution
   :class: dropdown

   Differentiate 10 times:

   .. jupyter-execute::

      long_expr = base_expr.diff(x, 10)

   Create the numerical functions:

   .. jupyter-execute::

      eval_long_expr = sm.lambdify((a, b, c, x, y, z), long_expr)
      eval_long_expr_cse = sm.lambdify((a, b, c, x, y, z), long_expr, cse=True)

   Now time each function:

   .. jupyter-execute::

      %%timeit
      eval_long_expr(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)

   .. jupyter-execute::

      %%timeit
      eval_long_expr_cse(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)

Learn more
==========

This section only scratches the surface of what SymPy can do. The presented
concepts are the basic ones needed for this book, but getting more familiar
with SymPy and what it can do will help. I recommend doing the `SymPy
Tutorial`_. The "Gotchas" section is particularly helpful for common mistakes
when using SymPy. The tutorial is part of the SymPy documentation
https://docs.sympy.org, where you will find general information on SymPy.

.. _SymPy Tutorial: https://docs.sympy.org/latest/tutorial/index.html

The tutorial is also available on video:

.. raw:: html

   <iframe width="560" height="315"
   src="https://www.youtube.com/embed/AqnpuGbM6-Q" title="YouTube video player"
   frameborder="0" allow="accelerometer; autoplay; clipboard-write;
   encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

If you want to ask a question about using SymPy (or search to see if someone
else has asked your question), you can do so at the following places:

- `SymPy mailing list <https://groups.google.com/g/sympy>`_: Ask questions via
  email.
- `SymPy Github Discussions <https://github.com/sympy/sympy/discussions>`_: Ask
  questions via Github.
- `Stackoverflow
  <https://stackoverflow.com/questions/tagged/sympy?tab=Votes>`_: Ask and
  search questions on the most popular coding Q&A website.
