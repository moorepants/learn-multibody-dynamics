=======
Vectors
=======

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`vectors` or Jupyter Notebook:
   :jupyter-download:notebook:`vectors`.

.. jupyter-execute::

   import sympy as sm
   sm.init_printing(use_latex='mathjax')

What is a vector?
=================

Vectors have three characteristics:

1. magnitude
2. orientation
3. sense

The direction the vector points is derived from both the orientation and the
sense.

In this text we will distinguish scalar variables, e.g. :math:`v`, from vectors
by including a bar over the top of the symbol, e.g. :math:`\bar{v}`. Vectors
are equal when all three characteristics are the same.

Vectors have these mathematical properites:

- scalar mutiplicative: :math:`\lambda\bar{a}` where :math:`\lambda` can only
  change the magnitude and the sense of the vector
- commutative: :math:`\bar{a} + \bar{b} = \bar{b} + \bar{a}`
- distributive: :math:`\lambda(\bar{a} + \bar{b}) = \lambda\bar{a} +
  \lambda\bar{b}`
- assocative: :math:`(\bar{a} + \bar{b}) + \bar{c} = \bar{a} + (\bar{b} +
  \bar{c})`
- common orientation: :math:`\bar{b} = k\bar{a}` where :math:`\bar{b}` and
  :math:`\bar{b}` have the same orientation

Unit vectors are vectors with a magnitude of :math:`1`. Any vector has an
assocated unit vector with the same orientation and sense, found by:

.. math::

   \hat{v} = \frac{\bar{v}}{||\bar{v}||}

where :math:`||\bar{v}||` is the `Euclidean norm`_ (2-norm) of the vector.

.. _Euclidean norm: https://en.wikipedia.org/wiki/Norm_(mathematics)#Euclidean_norm

Vector Functions
================

Vectors can be functions of scalar variables. :math:`\bar{v}` is a vector
function of scalar variable :math:`q` in reference frame :math:`A` if, when
:math:`q` changes, :math:`\bar{v}` changes when viewed in :math:`A`. This
implies that :math:`\bar{v}` may not be a function of scalar variable :math:`q`
in a different reference frame.

Staring with reference frame :math:`A` and vector :math:`\bar{v}` which is a
function of :math:`n` scalars :math:`q_1,q_2,\ldots,q_n` in :math:`A`, let
:math:`\hat{a}_x,\hat{a}_y,\hat{a}_z` be a set of mutually perpendicular unit
vectors fixed in :math:`A`, i.e. :math:`\hat{a}_x,\hat{a}_y,\hat{a}_z` are
constant when observed from :math:`A`. Then there are three unique scalar
functions :math:`v_x,v_y,v_z` of :math:`q_1,q_2,\ldots,q_n` such that:

.. math::

   \bar{v} = v_x \hat{a}_x + v_y \hat{a}_y + v_y \hat{a}_y

:math:`v_x \hat{a}_x` is called the :math:`\hat{a}_x` component of
:math:`\bar{v}` and :math:`v_x` is called measure number of :math:`\bar{v}`.
Since the components are mutually perpendicular the measure number is the dot
product
