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
the book keeping becomes extremely tedious and error prone for systems with
even a small number of bodies. SymPy allows us to let a computer handle the
tedious aspects, e.g. differentiation, solving linear systems of equations, and
reduce the errors one would encounter with pencil and paper. This chapter
introduces SymPy and the primary core SymPy features needed we will be using.

.. _SymPy: https://www.sympy.org

Import and setup
================

We will consistently import SymPy as follows:

.. jupyter-execute::

   import sympy as sm

Since SymPy lets you work with mathematical symbols it's nice to view SymPy
objects in a format that is similar to the math in a textbook. Executing
:external:py:func:`~sympy.interactive.printing.init_printing` at the beginning
of your Jupyter Notebook will ensure that SymPy objects render as typeset
mathematics. I use the ``use_latex='mathjax'`` argument here to disable math
image generation.

.. jupyter-execute::

   sm.init_printing(use_latex='mathjax')

Symbols
=======

Symbols are created with the :external:py:func:`~sympy.core.symbol.symbols`
function. A symbol :math:`a` is created like so:

.. jupyter-execute::

   a = sm.symbols('a')
   a

This symbol object is of the :external:py:class:`~sympy.core.symbol.Symbol` type:

.. jupyter-execute::

   type(a)

Multiple symbols can be created with one call to ``symbols()`` and SymPy
recognizes common Greek symbols by spelled out name.

.. jupyter-execute::

   b, t, omega = sm.symbols('b, t, omega')
   b, t, omega

Note that the argument provided to ``symbols()`` does not need to match the
Python variable name it is assigned to. Using more verbose Python variable
names may make code easier to read and understand, especially if there are many
mathematical variables that you forget the meaning of. Note that the subscripts
are recognized too.

.. jupyter-execute::

   pivot_angle, w2 = sm.symbols('alpha1, omega2')
   pivot_angle, w2

Undefined functions
===================

We will also work with undefined mathematical functions in addition to symbols.
These will play an important role in setting up differential equations, where
we typically don't know the function, but only its derivative(s). You can
create arbitrary functions of variables. In this case, we make a function of
:math:`t`. First create the function name:

.. jupyter-execute::

   f = sm.Function('f')
   f

This is of a type :external:py:class:`~sympy.core.function.UndefinedFunction`.

.. jupyter-execute::

   type(f)

Now we can create functions of one or more variables like so:

.. jupyter-execute::

   f(t)

Due to SymPy's internal implementations, the type of this object is not defined
as expected:

.. jupyter-execute::

   type(f(t))

so be aware of that.

The same ``UndefinedFunction`` can be used to create multivariate functions:

.. jupyter-execute::

   f(a, b, omega, t)

Symbolic expressions
====================

Now that we have mathematical variables and functions available, they can be
used to construct mathematical expressions. The most basic way to construct
expressions is with the standard Python operators ``+``, ``-``, ``*``, ``/``,
and ``**``. For example:

.. jupyter-execute::

   expr1 = a + b/omega**2
   expr1

An expression will have the type ``Add, Mul, or Pow``:

.. jupyter-execute::

   type(expr1)

This is because SymPy stores expressions as a tree_. You can inspect this
internal representation by using the
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

Undefined functions can also be used in expressions:

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

Lastly, an expression of ``t``:

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

SymPy also has a "pretty printer" (:external:py:func:`pprint()
<sympy.printing.pretty.pretty.pretty_print>`) that makes use of unicode symbols
to provide a form that more closely resembles typeset math:

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

   When you are working with long expressions, which will be the case in this
   course, there is no need to print them to the screen. In fact, printing them
   to the screen make take a long time and fill your entire notebook with an
   unreadable mess.

Differentiating
===============

One of the most tedious tasks in formulating equations of motion is the
differentiation of complex trigonometric expressions. SymPy can calculate
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
variables, as in the following:

.. math::

   \frac{\partial^2 h(a, \omega, t, b)}{\partial t \partial b}

.. jupyter-execute::

   expr3.diff(b, t)

Note that the answer includes real and imaginary components and the `signum
function`_.

.. _signum function: https://en.wikipedia.org/wiki/Sign_function

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

   expr3.evalf(n=31, subs=repl)

.. jupyter-execute::

   type(expr3.evalf(n=31, subs=repl))

Note that this is a SymPy :external:py:class:`~sympy.core.numbers.Float` object, which is a special object that can
have an arbitrary number of decimal places, for example here is the expression
evaluated to 300 decimal places:

.. jupyter-execute::

   expr3.evalf(n=300, subs=repl)

To convert this to Python floating point number, use ``float()``:

.. jupyter-execute::

   float(expr3.evalf(n=300, subs=repl))

.. jupyter-execute::

   type(float(expr3.evalf(n=300, subs=repl)))

This value is a machine precision floating point value and can be used with
standard Python functions that operating on floating point numbers.

To obtain machine precisions floating point numbers directly, it is better to
use the :external:py:func:`~sympy.utilities.lambdify.lambdify` function to convert the expression into a Python
function:

.. jupyter-execute::

   eval_expr3 = sm.lambdify((omega, a, f(t), b), expr3)

.. jupyter-execute::

   help(eval_expr3)

Now you have a function that operates on and returns floating point values:

.. jupyter-execute::

   eval_expr3(3.14/4, 2, -12, 25)

.. jupyter-execute::

   type(eval_expr3(3.14/4, 2, -12, 25))

This distinction between SymPy ``Float`` objects and regular Python and NumPy
``float`` objects is important. In this case, the Python float and the NumPy
float are equivalent. The later will compute much faster because arbitrary
precision is not required.

.. note::

   In these materials, you will almost always want to convert SymPy expressions
   into machine precision floating point numbers, so use ``lambdify()`` almost
   exclusively.

Matrices
========

SymPy supports matrices of expressions and linear algebra. Many of the
operations needed in multibody dynamics are more succinctly formulated with
matrices and linear algebra. Matrices can be created by passing nested lists to
the :external:py:class:`Matrix() <sympy.matrices.dense.MutableDenseMatrix>` object. For example:

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

Element-by-element multiplication requires the
:external:py:func:`~sympy.matrices.expressions.hadamard.hadamard_product`
function:

.. jupyter-execute::

   sm.hadamard_product(mat1, mat2)

Differentiation operates on each element of the matrix:

.. jupyter-execute::

   mat3 = sm.Matrix([expr1, expr2, expr3, expr4, expr5])
   mat3

.. jupyter-execute::

   mat3.diff(a)

.. jupyter-execute::

   mat3.diff(t)

The Jacobian_ matrix of vector (column matrix) can be formed with the
:external:py:meth:`~sympy.matrices.matrices.DenseMatrix.jacobian` method. This
calculates the partial derivatives of each element in the vector with respect
to a vector (or sequence) of variables.

.. jupyter-execute::

   mat4 = sm.Matrix([a, b, omega, t])
   mat4

.. jupyter-execute::

   mat3.jacobian(mat4)

.. _Jacobian: https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant

.. _sec-solving-linear-systems:
Solving Linear Systems
======================

You'll need to solve linear systems of equations often in this course. SymPy
offers a number of ways to do this, but the best way to do so if you know a set
of equations are linear in specific variables is the method described below.
First, you should know you have equations of this form:

.. math::

   a_{11} x_1 + a_{12} x_2 + \ldots + a_{1n} x_n + b_1 = 0 \\
   a_{21} x_1 + a_{22} x_2 + \ldots + a_{2n} x_n + b_2 = 0 \\
   \ldots \\
   a_{n1} x_1 + a_{n2} x_2 + \ldots + a_{nn} x_n + b_n = 0

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
   \end{bmatrix}

   \bar{x} =
   \begin{bmatrix}
     x_1 \\
     x_2 \\
     \ldots \\
     x_n
   \end{bmatrix}

   \bar{b} =
   \begin{bmatrix}
     -b_1 \\
     -b_2 \\
     \ldots \\
     -b_n
   \end{bmatrix}

Finally, :math:`\bar{x}` is found with matrix inversion (if the matrix is
invertible):

.. math::

   \bar{x} = \mathbf{A}^{-1}\bar{b}

Taking the inverse is not computationally efficient, so some form of `Gaussian
elmination`_ should be used to solve the system.

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

The :math:`\bar{b}` vector can be formed by setting :math:`a=b=0`, leaving the
terms that are not linear in :math:`a_1` and :math:`a_2`.

.. jupyter-execute::

   b = -exprs.xreplace({a1: 0, a2:0})
   b

Lastly, the ``LUsolve()`` method performs Gaussian-Elimination to solve the
system:

.. jupyter-execute::

   A.LUsolve(b)

Simplification
==============

The above result from
:external:py:meth:`~sympy.matrices.matrices.MutableDenseMatrix.LUsolve` is a bit
complicated. SymPy has some functionality for automatically simplifying
symbolic expressions. The function
:external:py:func:`~sympy.simplify.simplify.simplify` will attempt to find a
simpler version:

.. jupyter-execute::

   sm.simplify(A.LUsolve(b))

But you'll have the best luck at simplifying if you use specific functions that
target what type of expression you may have. The
:external:py:func:`~sympy.simplify.trigsimp.trigsimp` function only attempts
trigonometric simplifications, for example:

.. jupyter-execute::

   sm.trigsimp(sm.cos(omega)**2 + sm.sin(omega)**2)

.. warning::

   Only attempt simplification on expressions that are several lines of text.
   Larger expressions become increasingly computationally intensive to simplify
   and there is generally no need to do so in these materials.

As mentioned earlier, SymPy represents expressions as graphs (trees). Symbolic
expressions can also be represented as `directed acyclic graphs`_ that contain
only one node for each unique expression (unlike SymPy's trees which may have
repeated expressions in nodes). These unique expressions, or "common
sub-expressions", can be found with the
:external:py:func:`~sympy.simplify.cse_main.cse` function. This function will
provide a simpler form of the equations that minimizes the number of operations
to compute the answer.

.. _Directed acyclic graphs: https://en.wikipedia.org/wiki/Directed_acyclic_graph

.. jupyter-execute::

   substitutions, simplified = sm.cse(A.LUsolve(b))

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

Learn more
==========

This section only scratches the surface of what SymPy can do. The presented
concepts are the basic ones needed for this course, but getting more familiar
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
- `SymPy Gitter <https://gitter.im/sympy/sympy>`_: Ask questions in a live
  chat.
- `Stackoverflow
  <https://stackoverflow.com/questions/tagged/sympy?tab=Votes>`_: Ask and
  search questions on the most popular coding Q&A website.
