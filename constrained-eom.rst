===============================
Constrained Equations of Motion
===============================

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`constrained-eom` or Jupyter Notebook:
   :jupyter-download:notebook:`constrained-eom`.

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   me.init_vprinting(use_latex='mathjax')

When there are holonomic constraints present the equations of motion are
comprised of the kinematical differential equations
