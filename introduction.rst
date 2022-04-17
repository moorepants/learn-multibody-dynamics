============
Introduction
============

What You Will Learn
===================

- Formulate the equations of motion for a set of interacting rigid bodies.
- Incorporate kinematic constraints.
- Simulate a multibody system.
- Visualize the motion of a multibody system.
- Interpret the behavior of multibody systems.

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
integrated. I hope to weave the two much more fluidly here. This means that we
need to choose a computing language to do so.  There are many programming
languages well suited to multibody dynamics computation, but we select Python_
for a few reasons: 1) Python is open source and freely available for use, 2)
Python is currently one, if not the, most popular programming language in the
world, 3) the scientific libraries available in Python are voluminous and
widely used in academia and industry, and 4) Python has SymPy_ which provides a
foundation for computer aided algebra and calculus.

.. _Python: http://www.python.org
.. _SymPy: http://www.sympy.org

This text is not about teaching deep theory in multibody dynamics, but is
focused on application and doing. After following the text you should be able
to correctly model, simulate, and visualize multibody dynamic systems so that
you can use them as a tool to answer the engineering and scientific questions.
There are plenty of other books to learn the theory and we will cite them so
you can study further. In particular, this text provides an alternative
approach to understanding and utilizing the dynamics approach given by Thomas
R. Kane and David Levinson in their 1985 book "Dynamics, Theory and
Application" [Kane1985]_.

History
=======

The primary presentation of multibody dynamics in this text is based on the
presentation I and my fellow graduate students received in the graduate
Multibody Dynamics course taught by Mont Hubbard and the late Fidelis O. Eke at
the University of California, Davis from sometime in the 80s until 2016. Prof.
Eke was a PhD student of Thomas Kane at Stanford and Prof. Hubbard adopted
Prof. Kane's approach to dynamics after moving to UC Davis from Stanford[*]_.
The 10 week course was based on Thomas Kane's and David Levinsion's 1985 book
"Dynamics, Theory and Application" and followed the book and companion
computational materials closely.

In 2006, I took the course and naively decided to derive and analyze the
non-linear and linear Carvallo-Whipple bicycle model (see [Meijaard2007]_ for
an explanation of this model) as my course project. Fortunately, another
student visiting from Aachen university also choose the same model and his
success finally helped me squash the bugs in my model. Luke Peterson, Gilbert
Gede, and Angadh Nanjangud subsequently joined Hubbard and Eke's labs and with
Luke's lead we began being sucked into the world of open source scientific
software. At that time, Python's use by scientists and engineers began to
rapidly climb and we fortunately jumped on the bandwagon. We became quite
frustrated with the black box approach of the commercial software tools most
engineers used at that time, this included the tool Autolev that was developed
by Kane's collaborators for the automation of multibody dynamics modeling. Luke
wrote the first version of PyDy as a Google Summer of Code participant in 2009.
Gilbert followed him by implementing a new version as SymPy Mechanics in 2011
also as a Google SUmmer of Code participant. We Gilbert's implementation in
this this text. Combined with the power of SymPy and Jupyter Notebooks (IPython
Notebooks back then), SymPy Mechanics provides a computational tool that is
especially helpful for learning and teaching multibody dynamics.  It is also
equally useful for advanced modeling in research and industry.

I've stewarded and developed the software as well as taught and researched with
it over the last decade with the help of a long list of contributors. This text
is a presentation of the methods and lessons learned from over the last 10
years of doing multibody dynamics with open source Python software tools.

Acknowledgements
================

These are the primary contributors to the software presented in the text, in
approximate order of first contribution:

.. todo:: Find all the info and dates.

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

Peter Stahlecker and Jan Heinen provided page by page review of the while
drafting the first version. Their feedback has helped improve the text in many
ways.

.. [*] Mont was working on a skateboard dynamics model in the the late 70s and
   presented his model to an audience that included Thomas Kane. As the story
   goes, Prof. Kane approached Mont after the lecture to privately tell him his
   dynamics model was incorrect. Mont then took it upon himself to learn Kane's
   approach to dynamics so that his future models would be less likely to have
   such errors.
