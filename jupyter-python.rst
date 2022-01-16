Introduction
============

The following is a brief introduction to Python and the IPython
Notebook. There is much more to learn than what is covered here. This is
just enough to get you started with the tutorial.

The IPython Notebook
====================

IPython consists of two processes: a kernel and a frontend. The kernel
executes code while the frontend provides an interface for the user
enter their code. The IPython notebook is a frontend for Python which
provides an interactive web-environment for executing code and
displaying rich media.

Using the Notebook
~~~~~~~~~~~~~~~~~~

To start an IPython notebook session open a command prompt, navigate to
a desired working directory then issue the following command:

``ipython notebook``

A new window will open in your web browser where you can open an
existing notebook or start a new one. Notebooks are organized with
cells. You may have a code cell for inputing commands followed by its
result cell which contains the output of the code. You may also have a
text cell, such as this one you’re reading right now. Cell type can be
changed in the above dropdown menu.

There is the menubar above for navigating a notebook but you will find
the following shortcuts helpful:

-  ``Enter`` : Create a new line with in cell
-  ``Shift + Enter`` : Execute cell and advance to next cell
-  ``Ctrl + Enter`` : Execute cell in place (do not advance to the next
   cell)
-  Press ``esc`` (command mode) then ``h`` to display keyboard shortcuts

At times you might run code that gets stuck in an infinite loop or
simply you want to clear all your workspace variables and start over. To
solve each of these problems you can click on the menu:

``Kernel -> Interrupt``

``Kernel -> Restart``

Magic Commands
~~~~~~~~~~~~~~

These are commands to control IPython itself.

.. jupyter-execute::

    # list available magic commands
    %lsmagic

Need Help?
~~~~~~~~~~

In case you’re lost help isn’t far. The following commands should
provide assistance.

.. jupyter-execute::

    # Displays an overview of IPython's features
    ?

.. jupyter-execute::

    # A Quick Reference for IPython
    %quickref

.. jupyter-execute::

    # For details about an object
    # object_name?
    round??

.. jupyter-execute::

    # list functions available in your workspace
    dir()

.. jupyter-execute::

    # list variables available in your workspace
    %whos

Python
======

Basic Data Types
~~~~~~~~~~~~~~~~

.. jupyter-execute::

    a = 5
    b = 5.0
    c = float(5)
    d = 'dee'
    e = 'e'
    
    type(a), type(b), type(c), type(d), type(e)

Data Structures
~~~~~~~~~~~~~~~

Python offers several builtin data structures for arranging data in
memory. We will be making use of lists, dictionaries, tuples during this
tutorial.

Lists
^^^^^

A list is a versatile container that holds objects in the order given.
Lists are typically used to group similar items but may contain
heterogenous data types.

.. jupyter-execute::

    empty_list = []
    
    string_list = ['lions', 'tigers', 'bears', 'sharks', 'hamsters']
    
    int_list = [0, 1, 2, 3, 4]
    int_list2 = range(5,10)
    
    list_from_variables = [a,b,c,d,e]
    
    list_of_lists = [empty_list,
                     string_list,
                     list_from_variables,
                     int_list,
                     int_list2]
    
    print(list_of_lists)

Elements of a list are accessible by their index.

.. jupyter-execute::

    print(string_list[0])
    print(string_list[1:4])
    print(int_list[::2])  # get every 2nd element
    print(list_of_lists[1][4])  # get a nested item

List are mutable, meaning after a list is created we can change, add, or
remove elements.

.. jupyter-execute::

    int_list[2] = 222
    
    int_list.append(5)
    
    string_list.remove('lions')
    
    list_from_variables.extend(int_list)
    
    print(int_list)
    print(string_list)
    print(list_from_variables)

Tuples
^^^^^^

Tuples share similarites with lists. A tuple is good for organizing
related data that may be of different types. Notice they are defined
with parenthesis, ``()``, rather than brackets.

.. jupyter-execute::

    joe_blow = (32, 'tall', 'likes hats')
    
    print(joe_blow[1])

Unlike lists, tuples are immutable. They cannot be changed once defined.

.. jupyter-execute::

    # this won't work
    #joe_blow.append('married')
    
    # neither will this
    #joe_blow[2] = 'not really a fan of hats'

In python a function can return multiple values. These ouputs are packed
into a tuple. Tuple unpacking assigns individual elements of a tuple to
separate variables.

.. jupyter-execute::

    pets = ('elephant', 'cow', 'rock')
    
    pet1, pet2, pet3 = pets

A peculiar thing about tuples in python is defining a single element
tuple. Note the trialing comma. This is necessary for python to know you
want a one-element tuple.

.. jupyter-execute::

    (pet1,)

Dictionaries
^^^^^^^^^^^^

A dictionary is an unordered set of *key:value* pairs. Much like a
language dictionary where you look up a *word* and get its *definition*
in a python dictionary you look up a *key* and get its *value*.

.. jupyter-execute::

    # numbers or strings may be used as keys
    dictionary0 = {'key1':'value1', 'key2':'value2', 'key3':'value3'}
    dictionary1 = {1:'value1', 2:'value2', 3:'value3'}
    
    cylinder = {'mass':50, 'base':10, 'height':100}
    
    print(dictionary0)
    print(dictionary1.keys())
    print(cylinder['mass'])

The zip function is a convenient function to help generate a dictionary.
It takes sequence objects and combines them into a list of tuples. We
can subsequently use the list of two element tuples to create a
dictionary.

.. jupyter-execute::

    keys = ['mass01', 'inertia01', 'mass02', 'inertia02']
    values = [10, 1, 50, 5]
    
    dict(zip(keys, values))

Functions
~~~~~~~~~

Python does not use braces, ``{}``, or ``end`` statements to seperate
blocks of code. Rather, code blocks are initialized with colon, ``:``,
and defined by their indentation. It is convention to use four spaces
for each level of indentation.

.. jupyter-execute::

    def abs_value(A):
        if A < 0:
            A = -A 
        return A
    
    abs_value(-100)

.. jupyter-execute::

    def long_div(dividend, divisor):
        quotient = dividend // divisor # // : floor division
        remainder = dividend % divisor # % : modulo
        return quotient, remainder
    
    
    a = 430
    b = 25
    
    # an example of tuple unpacking
    quo, rem = long_div(a, b)
    
    print('%d divided %d is %d remainder %d' % (a, b, quo, rem))

Modules
~~~~~~~

Modules add additional functionality not present in the default
installation of python. Throughout this tutorial we will either import
an entire module or import specific functions from a module.

.. jupyter-execute::

    # import object from sympy into the current namespace
    from numpy import array
    
    # import multiple objects from sympy
    from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point

Objects from these modules are now available in your namespace.

.. jupyter-execute::

    # from numpy
    arr = array([1,2,3,4,5])
    
    # from sympy.physics.mechanics
    inertial_frame = ReferenceFrame('I')

The SciPy Stack
===============

SciPy is a collection of open source software that brings scientific
computing to Python.

-  IPython : Enhanced Interactive Console
-  Numpy : N-dimensional array package
-  Matplotlib : a fully featured 2D plotting package
-  Scipy : a collection of toolboxes for signal processing,
   optimization, statisitics, etc.
-  Sympy : symbolic mathematics and computer algebra
-  and much more…

Numpy
-----

Setup
~~~~~

.. jupyter-execute::

    from numpy import random, linspace, zeros, arange

Creating Arrays
~~~~~~~~~~~~~~~

.. jupyter-execute::

    # An array from a list
    print(array([5, 12, -2, 9.3, 7]))
    
    # random values
    print(random.random((5)))
    
    # linearly spaced values
    # 5 values between 0 and 10 inclusive
    print(linspace(0,10,5))
    
    # range of values with a defined stepsize
    # start at 0 and increase by 3 
    print(arange(0,14,3))

Accessing Array Elements
~~~~~~~~~~~~~~~~~~~~~~~~

.. jupyter-execute::

    P = random.random((3,5))
    
    # individual element
    print(P[0,3])
    
    # entire row
    print(P[2])
    
    # entire column
    print(P[:,4])
    
    # every third element
    print(P[::3])

Operations on Arrays
~~~~~~~~~~~~~~~~~~~~

.. jupyter-execute::

    # mathematical operations are always elementwise
    x = arange(5)
    
    print(x)
    
    # the double asterisk represents exponentiation
    print(x + x**2)

Matplotlib
----------

Matplotlib provides an API similar to MATLAB’s

Use the magic command ``%matplotlib inline`` to work with matplotlib
interactively. The ``inline`` argument allows for plots to be embedded
within the IPython notebook

.. jupyter-execute::

    %matplotlib inline

.. jupyter-execute::

    from matplotlib.pyplot import plot, subplot, xlabel, ylabel, legend, tight_layout
    from numpy import sin, cos, pi

Examples
~~~~~~~~

.. jupyter-execute::

    x = arange(-pi,pi,0.1)
    y1 = 2*sin(x)
    y2 = x + cos(4*x)
    
    plot(x, y1, 'r', x, y2, '--b')
    plot(x[::5], y1[::5], 'og') # plot every 5th point
    xlabel('x axis')
    ylabel('y axis')
    legend(['red', 'blue', 'green'], loc='upper left')

.. jupyter-execute::

    x = linspace(-100,100)
    
    for i in range(1,5):
        subplot(2,2,i)
        plot(x, x**i)
        
    tight_layout() # this prevents the axis labels from overlapping

Exercise
~~~~~~~~

Use the provided function to create 3 different sine waves at various
frequencies. Plot the 3 functions with labeled axis.

.. jupyter-execute::

    def three_sine_waves(t, A, B, C):
        """
        t : (type: array) an monotonically increasing array of time values
        A,B,C : (type: float) frequency of sine waves
        """
        y1 = sin(A*t)
        y2 = sin(B*t)
        y3 = sin(C*t)
        
        return y1, y2, y3

.. jupyter-execute::

    #t = 
    #y1, y2, y3 = three_sine_waves()

Scipy odeint
------------

Scipy provides a routine for integrating first order ordinary
differential equations.

Setup
~~~~~

.. jupyter-execute::

    from scipy.integrate import odeint

Examples
~~~~~~~~

.. jupyter-execute::

    def dy(y,x):
        return x
    
    y0 = 0.0
    x = linspace(-5.0, 5.0, 1000)
    
    y = odeint(dy,y0,x)
    
    plot(x,y)

.. jupyter-execute::

    def dy(y,t,coeff):
        A = coeff['A']
        B = coeff['B']
        C = coeff['C']
        D = coeff['D']
        return A*t**3 + B*t**2 + C*t + D
    
    y0 = 2.0
    t = linspace(-5.0, 3.0, 1000)
    sys = {'A' : 0.25,
           'B' : 0.75,
           'C' : -1.5,
           'D' : -2.0}
    
    y = odeint(dy, y0, t, args=(sys,))
    
    plot(t,y)

SymPy
-----

Setup
~~~~~

.. jupyter-execute::

    from sympy import *
    
    # This is for prettier printing of equations
    interactive.printing.init_printing()

Creating Symbolic Variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. jupyter-execute::

    a = symbols('a')
    b = symbols('b')
    gravity, mass, spring_const, time = symbols('g, m, k, t')

Expressions Using Symbolic Variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. jupyter-execute::

    a**2 + b/pi

.. jupyter-execute::

    simplify(4*(a*a*a)*(b+b+b))

.. jupyter-execute::

    diff(-gravity*time**2/2, time)

.. jupyter-execute::

    # indefinte integral
    
    integrate(-gravity,time)

.. jupyter-execute::

    # definite integral
    
    v0 = 5
    t1 = 0
    t2 = .35
    position = integrate(-gravity*time + v0,(time,t1,t2))
    position.subs(gravity, 9.81)

Additional Resources
====================

[1] http://docs.python.org/2/tutorial

[2]
http://nbviewer.ipython.org/github/ipython/ipython/blob/master/examples/notebooks/Cell%20Magics.ipynb

[3] http://www.scipy.org/index.html
