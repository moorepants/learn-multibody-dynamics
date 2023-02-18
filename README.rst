======
README
======

This repository contains the source for the website `Learn Multibody Dynamics
<https://moorepants.github.io/learn-multibody-dynamics/>`_.

Authors and Contributors
========================

Authors: Jason K. Moore

Contributors: Peter Stahlecker, Jan Heinen, @tamoora, Christopher Dembia, Arthur Ryman, Wouter Wolfslag

License
=======

The contents of this repository are licensed under the CC-BY 4.0 license. See
``license.rst`` for more information.

Design
======

The design of the book has these general goals:

- use computational thinking to teach multibody dynamics
- the book is online and interactive (run/modify code, interact with figures)
- don't just teach the math and expect the students to figure out the
  computation on their own, the computations should be first class and could
  even go as far as replacing the math as the teach mechanism
- computational difficulties should not be hidden and should be taught
- after completing the book, you should be able to create a physics engine if
  you desired
- most lessons should be motivated by a real problem statement, instead of
  teaching the concepts with little connection to reality or strictly from a
  mathematical sense
- allow collaborative contributions

Building the Website
====================

Setup the code and environment
------------------------------

Clone the repository::

   git clone https://github.com/moorepants/learn-multibody-dynamics.git
   cd learn-multibody-dynamics

Install miniconda_ or a similar tool (e.g.  Anaconda_) and create a conda
environment for the book::

   conda env create -f multibody-book-env.yml

.. _miniconda: https://docs.conda.io/en/latest/miniconda.html
.. _Anaconda: https://www.anaconda.com/products/individual

Activate the conda environment::

   conda activate multibody-book

Build and view the website
--------------------------

To build the website run::

   make html

When complete, the website is then viewable in your browser::

   <yourbrowser> _build/html/index.html

You can also run sphinx-autobuild (updates while while you edit) with::

   make autobuild

Working with a remote git branch
--------------------------------

If you want to build one of the branches (for example a pull request), you'll
need to fetch and checkout the branch. First fetch down all the branches::

   git fetch origin

Then checkout the branch (this command is only need the first time you check it
out)::

   git checkout -b branch-name origin/branch-name

The branch name is listed on the pull request just under the title "...wants to
merge X commits into master from branch-name." Or you can find all branches
here: https://github.com/moorepants/learn-multibody-dynamics/branches

Now run::

   make clean
   make html

The ``make clean`` makes sure you don't keep any remnants from prior builds
around before building the new branch.

After you have a new branch setup you can switch between the master branch and
any branch name with just::

   git checkout master
   git checkout branch-name

If the master branch or any other branch has been updated on github you can
pull down the latest changes with::

   git checkout branch-name
   git pull origin branch-name

Editing Guide
=============

restructuredtext
----------------

The text is written in reStructuredText and processed with Sphinx. The `Sphinx
reStructuredText documentation page
<https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html>`_ is a
good starting point to learn the syntax.

reStructuredText doesn't enforce a specific heading order, but this should be
followed for this text:

.. code:: rst

   =======
   Chapter
   =======

   Section
   =======

   Subsection
   ----------

   Subsubsection
   ^^^^^^^^^^^^^

Autoreferencing is enabled so the above sections can be referenced with:

.. code:: rst

   :ref:`Chapter`
   :ref:`Section`
   :ref:`Subsection`
   :ref:`Subsubsection`

For equations and figures, they need to be manually labeled for numbered
referencing. Use these patterns:

.. code:: rst

   :label:`eq-my-equation-name`
   :math:numref:`eq-my-equation-name`

   .. _fig-my-figure-name:
   :numref:`fig-my-figure-name`

When citing Kane & Levinson 1985 or other books include the page number:

.. code:: rst

   ([Kane1985_], pg. 23)

Cross-referencing API documentation in various libraries:

.. code:: rst

   :external:py:meth:`~sympy.physics.vector.frame.ReferenceFrame.ang_acc_in`
   :external:py:class:`~sympy.physics.vector.frame.ReferenceFrame`
   :external:py:func:`~sympy.physics.vector.functions.cross`

Exercises look like this:

.. code:: rst

   .. admonition:: Exercise

      What is 1 + 1?

   .. admonition:: Solution
      :class: dropdown

      The answer is 2.

jupyer-sphinx
-------------

We use jupyter-sphinx to transform each page with code cells into a Jupyter
Notebook and Python script. Any page that includes ``.. jupyter-execute::``
directives will be processed in this way. The documentation for jupyter-sphinx
is here:

https://jupyter-sphinx.readthedocs.io

Each page that has executable code should include these download links at the
top of the page. If the filename is ``page.rst`` then include:

.. code:: rst

   .. note::

      You can download this example as a Python script:
      :jupyter-download-script:`page` or Jupyter Notebook:
      :jupyter-download-notebook:`page`.

Xournal++
---------

I draw the figures, one per page, in Xournal++. The I export as -> svg ->
choose None for background and "current page" to get a single exported svg.

The SVG figures should be cropped to the bounding box of the drawn elements.
One can do so using Inkscape with these button presses: File -> Document
Properties -> Resize Page to Content. With Inkscape > 1.0 this command will
crop the figure:

.. code:: bash

   inkscape --export-type=svg --export-area-drawing ./my-figure.svg

Live rebuilding with sphinx-autobuild
-------------------------------------

`Sphinx autobuild`_ is a pretty good way to get almost instaneous rendered HTML
versions of the reStructuredText file. You can open a window with your text
editor and a window with your broswer side-by-side for almost instant feedback
on the formatting and Jupyter code execution.

.. _Sphinx autobuild: https://github.com/executablebooks/sphinx-autobuild

.. code:: bash

   sphinx-autobuild -b html . _build/html/

This is also encoded in the Makefile:

.. code:: bash

   make autobuild

Execute code cells in IPython while writing
-------------------------------------------

tmux
^^^^

https://tmuxcheatsheet.com/

https://medium.com/hackernoon/a-gentle-introduction-to-tmux-8d784c404340

::

   tmux new
   <Ctrl>+b %  # side by side panes
   <Ctrl>+<arrow key>  # jump between panes

vim-slime
^^^^^^^^^

https://github.com/jpalardy/vim-slime

create a vim slime config file for rst

::

   <Ctrl>+cc  # execute line(s) in ipython pane

Content Resources
=================

Here are links to various resources that use SymPy for dynamics that could be
incorporated into this repository, as is or as inspiration:

- UC Davis MAE223 notebooks: https://moorepants.github.io/mae223/schedule.html
- PyDy tutorial: https://github.com/pydy/pydy-tutorial-human-standing
- PyDy documentation examples: https://pydy.readthedocs.io/en/latest/index.html#examples
- PyDy source examples: https://github.com/pydy/pydy/tree/master/examples
- SymPy Mechanics documentation: https://docs.sympy.org/dev/modules/physics/mechanics/index.html
- Resonance notebooks: https://moorepants.github.io/resonance/
- Yeadon example: https://nbviewer.jupyter.org/github/chrisdembia/yeadon/blob/v1.2.1/examples/bicyclerider/bicycle_example.ipynb
- Problems from EME 134: https://moorepants.github.io/eme134/labs.html
- TU Delft MAE41055 2021 homework notebooks
- Oliver's solutions to the TUD advanced dynamics course examples: https://github.com/pydy/pydy/pull/137
