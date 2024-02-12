============
Introduction
============

What You Will Learn
===================

- How to formulate the equations of motion for a set of interacting rigid
  bodies, i.e. a `multibody system`_.
- How to manage and incorporate kinematic constraints.
- How to simulate a multibody system.
- How to visualize the motion of a multibody system in 2D and 3D.
- How to interpret the behavior of multibody systems.

.. _multibody system: https://en.wikipedia.org/wiki/Multibody_system

Prerequisites
=============

- Linear algebra
- Vector calculus
- Calculus based physics
- Statics
- Dynamics
- Introductory numerical methods
- Introductory scientific computing

Purpose
=======

The goal of this text is to help you learn multibody dynamics via a mixture of
computation and traditional mathematical presentation. Most existing textbooks
on the subject are either purely mathematical and problems are solved by pencil
and paper or there are computational elements that are tacked on rather than
integrated. I hope to weave the two much more fluidly here so that you can
learn the principles of multibody dynamics through computing.

This text is less about teaching deep theory in multibody dynamics and more
about application and doing. After following the text and practicing, you
should be able to correctly model, simulate, and visualize multibody dynamic
systems so that you can use them as a tool to answer scientific questions and
solve engineering problems.

Choice of dynamics formalism
============================

To teach multibody dynamics, one must choose a formalism for notation and
deriving the equations of motion. There are numerous methods for doing so, from
`Newton and Euler`_'s to Lagrange_ and Hamilton_'s to Jain and Featherstone_'s.
Here I use an approach primarily derived from `Thomas R. Kane`_ and David
Levinson in their 1985 book "Dynamics, Theory and Application" [Kane1985]_. The
notation offers a precise way to track all of the nuances in multibody dynamics
bookkeeping and a realization of the equations of motion that obviates having
to introduce virtual motion concepts and that handles kinematic constraints
without the need of `Lagrange multipliers`_.

.. _Newton and Euler: https://en.wikipedia.org/wiki/Newton%E2%80%93Euler_equations
.. _Lagrange: https://en.wikipedia.org/wiki/Lagrangian_mechanics
.. _Hamilton: https://en.wikipedia.org/wiki/Hamiltonian_mechanics
.. _Featherstone: https://en.wikipedia.org/wiki/Featherstone%27s_algorithm
.. _Thomas R. Kane: https://en.wikipedia.org/wiki/Thomas_R._Kane
.. _Lagrange multipliers: https://en.wikipedia.org/wiki/Lagrange_multiplier

Choice of programming language
==============================

With the goal of teaching through computation, it means I need to also choose a
programming language. There are many programming languages well suited to
multibody dynamics computation, but I select Python_ for several reasons: 1)
Python is open source and freely available for use, 2) Python is currently one,
if not the, most popular programming language in the world, 3) the scientific
libraries available in Python are voluminous and widely used in academia and
industry, and 4) Python has SymPy_ which provides a foundation for computer
aided-algebra and calculus.

.. _Python: http://www.python.org
.. _SymPy: http://www.sympy.org

History
=======

The primary presentation of multibody dynamics in this text is based on the
presentation I and my fellow graduate students received in the graduate
Multibody Dynamics course taught by Mont Hubbard and the late Fidelis O. Eke at
the University of California, Davis in the early 2000's. Profs. Hubbard an Eke
taught the course from the late 80s or early 90s until they retired in 2013
(Prof. Hubbard) and 2016 (Prof. Eke). The 10-week course was based on Thomas R.
Kane's and David A. Levinson's 1985 book "Dynamics, Theory and Application" and
followed the book and companion computational materials closely. Prof. Eke was
a PhD student of Thomas R. Kane at Stanford and Prof. Hubbard adopted Prof.
Kane's approach to dynamics after moving to UC Davis from Stanford [#]_. I
helped with Prof. Eke's 2015 course and taught the course in `2017
<https://moorepants.github.io/mae223/2017/>`_ and `2019
<https://moorepants.github.io/mae223/>`_ at UC Davis and this text is a
continuation of the notes and materials I developed based on Profs. Hubbard and
Eke's notes and materials which now includes some elements of TU Delft's past
multibody dynamics course.

When I took the UC Davis course in 2006 as a graduate student, I naively
decided to derive and analyze the nonlinear and linear Carvallo-Whipple bicycle
model [Meijaard2007]_ as my course project [#]_. Fortunately, another student
visiting from Aachen University, Thomas Engelhardt, also choose the same model
and his success finally helped me squash the bugs in my formulation. Luke
Peterson, Gilbert Gede, and Angadh Nanjangud subsequently joined Hubbard and
Eke's labs and with Luke's lead we were sucked into the world of open source
scientific software. At that time, Python's use by scientists and engineers
began to gain traction and we fortunately jumped on the bandwagon. We had
become quite frustrated with the black box approach of the commercial software
tools most engineers used at that time, this included the tool Autolev that was
developed by Kane's collaborators for the automation of multibody dynamics
modeling. To remedy this frustration, Luke wrote the `first version of PyDy`_
as a `Google Summer of Code`_ participant in 2009. Gilbert followed him by
implementing a new version as `SymPy Mechanics`_ in 2011 also as a Google
Summer of Code participant. We use Gilbert's, now modified and extended,
implementation in this text. Combined with the power of SymPy and Jupyter
Notebooks (IPython Notebooks back then), SymPy Mechanics provides a
computational tool that is especially helpful for learning and teaching
multibody dynamics. It is also equally useful for advanced modeling in research
and industry.

.. _first version of PyDy: https://github.com/hazelnusse/pydy
.. _Google Summer of Code: https://en.wikipedia.org/wiki/Google_Summer_of_Code
.. _SymPy Mechanics: https://docs.sympy.org/latest/modules/physics/mechanics/index.html

I have stewarded and developed the software as well as taught and researched
with it over the last decade with the help of a long list of contributors. This
text is a presentation of the methods and lessons learned from over the years
of doing multibody dynamics with open source Python software tools.

.. [#] The project is shared at https://github.com/moorepants/MAE-223
.. [#] Mont was working on a skateboard dynamics model in the late 70s and
   presented his model to an audience that included Thomas Kane. As the story
   goes, Prof. Kane approached Mont after the lecture to privately tell him his
   dynamics model was incorrect. Mont then took it upon himself to learn Kane's
   approach to dynamics so that his future models would be less likely to have
   such errors.

Acknowledgements
================

Wouter Wolfslag contributed the "Equations of Motion with the Lagrange Method"
chapter, "Alternatives for Representing Orientation" section, and reviewed
updates for version 0.2. Peter Stahlecker and Jan Heinen provided page-by-page
review of the text while drafting version 0.1. Peter did the same for version
0.2. Arthur Ryman contributed edits to the first version. Their feedback has
helped improve the text in many ways. We also thank the students of TU Delft's
Multibody Dynamics course who test the materials while learning.

These are the primary contributors to the SymPy Mechanics software presented in
the text, in approximate order of first contribution:

- Dr. Luke Peterson, 2009
- Dr. Gilbert Gede, 2011
- Dr. Angadh Nanjangud, 2012
- Tarun Gaba, 2013
- Oliver Lee, 2013
- Dr. Chris Dembia, 2013
- Jim Crist, 2014
- Sahil Shekhawat, 2015
- James McMillan, 2016
- Nikhil Pappu, 2018
- Sudeep Sidhu, 2020
- Abhinav Kamath, 2020
- Timo Stienstra, 2022
- Dr. Sam Brockie, 2023

SymPy Mechanics is built on top of SymPy, whose `1000+ contributors`_ have also
greatly helped SymPy Mechanics be what it is. Furthermore, the software sits on
the top of a large ecosystem of open source software written by thousands and
thousands of contributors who we owe for the solid foundation.

.. _1000+ contributors: https://github.com/sympy/sympy/blob/master/AUTHORS

Tools Behind the Book
=====================

I write the contents in plain text using the reStructuredText_ markup language
for processing by Sphinx_. The mathematics are rendered with MathJax_ in the
HTML version. I use the `Jupyter Sphinx`_ extension which executes the code in
each chapter as if it were a Jupyter notebook and embeds the Jupyter generated
outputs into the resulting HTML page. The extension also converts each chapter
into a Python script and Jupyter notebook for download. I use the `Material
Sphinx Theme`_ and `sphinx-togglebutton`_ for the dropdown information boxes. I
host the source for the book on Github_, where I use Github Actions to build
the website and push it to a Github Pages host using `ghp-import`_. I use
Github's issue tracker and pull request tools to manage tasks and changes. The
figures are drawn with a Wacom One tablet and the `Xournal++`_ application.

.. _reStructuredText: https://en.wikipedia.org/wiki/ReStructuredText
.. _Sphinx: https://www.sphinx-doc.org
.. _MathJax: https://www.mathjax.org
.. _Jupyter Sphinx: https://github.com/jupyter/jupyter-sphinx
.. _Material Sphinx Theme: https://github.com/bashtage/sphinx-material
.. _sphinx-togglebutton: https://github.com/executablebooks/sphinx-togglebutton
.. _Github: https://github.com
.. _ghp-import: https://github.com/c-w/ghp-import
.. _Xournal++: https://xournalpp.github.io
