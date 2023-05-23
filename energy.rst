================
Energy and Power
================

.. note::

   You can download this example as a Python script:
   :jupyter-download-script:`loads` or Jupyter Notebook:
   :jupyter-download-notebook:`loads`.

.. jupyter-execute::

   import sympy as sm
   import sympy.physics.mechanics as me
   me.init_vprinting(use_latex='mathjax')

.. container:: invisible

   .. jupyter-execute::

      class ReferenceFrame(me.ReferenceFrame):

          def __init__(self, *args, **kwargs):

              kwargs.pop('latexs', None)

              lab = args[0].lower()
              tex = r'\hat{{{}}}_{}'

              super(ReferenceFrame, self).__init__(*args,
                                                   latexs=(tex.format(lab, 'x'),
                                                           tex.format(lab, 'y'),
                                                           tex.format(lab, 'z')),
                                                   **kwargs)
      me.ReferenceFrame = ReferenceFrame

Learning Objectives
===================

After completing this chapter readers will be able to:

Introduction
============

So far we have invesitaged multibody systems from the perspect of forces and
their relationship to motion. It is also useful to understand these systems
from a power and energy perspective. Power is the time rate of change in work
done. Work is the eenrgy gained, disspated, or exchanged in a system.

.. power:: https://en.wikipedia.org/wiki/Power_(physics)

.. math::

   P = \frac{dW}{dt}

Knowing that work is a force dotted with a change in position, power can be
written as a force dotted with a velocity.

.. math::

   P = \bar{F} \cdot \bar{v}

Power can enter into a system, exit a system, or be exhanged within a system.

The time integral of power is work or energy. Energy of a multibody system can
be classified as kinetic, potential (conservative), or non-conservative.

.. math::

   E_k = m \bar{v} \cdot \bar{v} / 2  + \bar{\omega} \cdot \breve{I} \cdot \bar{\omega} / 2

.. math::

   E_p = 

Jumping
=======

.. jupyter-execute::

   mu, mt, mc, mf = sm.symbols('m_u, m_t, m_c, m_f')
   It, Ic = sm.symbols('I_t, I_c')
   kc, cc = sm.symbols('k_c, c_c')
