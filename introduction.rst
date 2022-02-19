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

The goal of this book is to help you learn multibody dynamics via a mixture of
computation and traditional mathematic presenation. Most existing textbooks on
the subject are either purely mathematical and problems are solved by pencil
and paper or there are computational elements that are tacked on either in
specific problems or sections of the book. We hope to weave the two much more
fluidly here. This means that we need to choose a computing language to do so.
There are many programming lanugages well suited to multibody dynamics
computation, but we select Python for a few reasons: 1) Python is open source
and freely available for use, 2) Python is currently one, if not the, most
popular programming language in the world, and 3) the scientific libraries
available in Python are volumnious and widely used in academia and industry.

This book is not about teaching deep theory in multibody dynamics, but is
focused on application and doing. After following the book you should be able
to correctly model, simulate, and visualize multibody dynamic systems so that
you can use them as a tool to answer the engineering and scientific questions
you need to. There are plently of other books to learn the theory and we will
cite them so you can study further.

History
=======

The primary presentation of multibody dynamics in this book is based on the
presentation in the graduate Multibody Dynamics course taught by Mont Hubbard
and the late Fidelis O. Eke at the University of California, Davis from
sometime in the 80s until 2016. Prof. Eke was a PhD student of Thomas Kane at
Stanford and Prof. Hubbard adopted Prof. Kane's approach to dynamics after
moving to UC Davis from Stanford[*]_. The 10 week course was based on Thomas
Kane's and David Levinsion's 1985 book "Dynamics, Theory and Application" and
followed the book and companion computational materials closely.

In 2006, I took the course and naively decided to derive and analyzie the
non-linear and linear Carvallo-Whipple bicycle model as my course project.
Fortunately, another student visiting from Aachen university also choose the
same model and his success finally helped me squash the bugs in my model. Luke
Peterson, Gilbert Gede, and Angadh Nanjangud subsequently joined Hubbard and
Eke's labs and with Luke's lead we began being sucked into the world of open
source scientific software. At that time Python's use by scientists and
engineers began to rapidly climb and we jumped on the bandwagon. We became
quite frustrated with the black box approach of the commercial software tools
most engineers used at that time, this included the tool Autolev that was
developed by Kane and his collaborators for the automation of multibody
dynamics modeling. Luke wrote the first version of PyDy as a Google Summer of
Code participant in 2009. Gilbert followed him by implementing a new version as
SymPy Mechanics in 2011 also as a Google SUmmer of Code participant. We use
this implementation in this this book. Combined with the power of SymPy and
Jupyter Notebooks (IPython Notebooks back then), it provides a computational
tool that is especially helpful for learning and teachign multibody dynamics
and can also be used for advanced modeling too.

I've stewarded and developed the software as well as taught and researched with
it over the last decade with the help of a long list of contributors. This book
is a presentation of the methods and lessons learned from over the last 10
years of doing multibody dynamics with open source Python software tools.

These are the primary contibutors to the software presented in the book, in
apparoximate order of first contribution:

.. todo:: Find all the info and dates.

- Dr. Luke Peterson, 2009
- Dr. Gilbert Gede, 2011
- Dr. Angadh Nanjangud, 2012
- Jim Crist, 2014
- Tarun Gaba
- Oliver Lee
- Sahil Shekhawat
- Dr. Chris Dembia
- James McMillan
- Sudeep Sidhu
- Nikhil Pappu

.. [*] Mont was working on a skateboard dynamics model in the the late 70s and
   presented his model to an audience that included Thomas Kane. As the story
   goes, Prof. Kane approached Mont after the lecture to privately tell him his
   dynamics model was incorrect. Mont then took it upon himself to learn Kane's
   approach to dynamics so that his future models would be less likely to have
   such errors.

https://github.com/hazelnusse/pydy
