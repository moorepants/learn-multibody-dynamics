======
README
======

The source for the website "Learning Multibody Dynamics". Viewable at:

https://moorepants.github.io/learn-multibody-dynamics/

Building the Website
====================

Clone the repository::

   git clone https://github.com/moorepants/learn-multibody-dynamics.git
   cd learn-multibody-dynamics

Install miniconda_ or Anaconda_.

.. _minconda: https://docs.conda.io/en/latest/miniconda.html
.. _Anaconda: https://www.anaconda.com/products/individual

Create a conda environment for the book::

   conda env create -f multibody-book-env.yml

Activate the environment::

   conda activate multibody-book

To build once run::

   make html

When complete, the website is then viewable in your browser::

   <yourbrowser> _build/html/index.html

You can also run sphinx-autobuild (updates while while you edit) with::

   make autobuild

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

License
=======

The contents of this repository are licensed under the CC-BY 4.0 license. See
``license.rst``.

Editing Guide
=============

restructuredtext
----------------

The text is written in reStructuredText and procesed with Sphinx.

Heading order:

.. code:: rst

   =====
   Title
   =====

   H1
   ==

   H2
   --

   H3
   ^^

jupyer-sphinx
-------------

Any page that includes ``.. jupyter-execute::`` directives will be processed
with Jupyter Sphinx.

https://jupyter-sphinx.readthedocs.io

Xournal++
---------

I draw the figures, one per page, in Xournal++.

Reize xournal++ svg exports to just the extents

seems to require gui to open (--without-gui doesn't work with verbs)

inkscape --verb=FitCanvasToDrawing --verb=FileSave --verb=FileQuit orientation-camera-gimbal.svg

Live rebuilding with sphinx-autobuild
-------------------------------------

https://github.com/executablebooks/sphinx-autobuild

::

   sphinx-autobuild -b html . _build/html/

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
