======
README
======

The source for the website "Learning Multibody Dynamics". Viewable at:

https://moorepants.github.io/learn-multibody-dynamics/

restructuredtext
================

The text is written in reStructuredText.

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

tmux
====

https://medium.com/hackernoon/a-gentle-introduction-to-tmux-8d784c404340

::

   tmux new
   <Ctrl>+b %  # side by side panes
   <Ctrl>+<arrow key>  # jump between panes

vim-slime
=========

create a vim slime config file for rst

::

   <Ctrl>+cc  # execute line(s) in ipython pane

sphinx-autobuild
================

::

   sphinx-autobuild -b html . _build/html/
