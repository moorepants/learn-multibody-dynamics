======
README
======

The source for the website "Learning Multibody Dynamics". Viewable at:

https://moorepants.github.io/learn-multibody-dynamics/

License
=======

The contents of this repository are licensed under the CC-BY 4.0 license.

Editing Guide
=============

restructuredtext
----------------

The text is written in reStructuredText and procesed with Sphinx.

Heading order::

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

tmux
----

https://medium.com/hackernoon/a-gentle-introduction-to-tmux-8d784c404340

::

   tmux new
   <Ctrl>+b %  # side by side panes
   <Ctrl>+<arrow key>  # jump between panes

vim-slime
---------

https://github.com/jpalardy/vim-slime

create a vim slime config file for rst

::

   <Ctrl>+cc  # execute line(s) in ipython pane

sphinx-autobuild
----------------

https://github.com/executablebooks/sphinx-autobuild

::

   sphinx-autobuild -b html . _build/html/

Xournal++
---------

Reize xournal++ svg exports to just the extents

seems to require gui to open (--without-gui doesn't work with verbs)

inkscape --verb=FitCanvasToDrawing --verb=FileSave --verb=FileQuit orientation-camera-gimbal.svg
