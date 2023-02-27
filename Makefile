# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build
CHAPTER       =

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

resizesvg:
	# this seems to work without opening the gui
	inkscape --export-plain-svg=test.svg --export-area-drawing ./generalized-forces-partial-velocities.svg
	# once i have inkscape 1.0 this should work
	#inkscape --export-type=svg --export-area-drawing ./generalized-forces-partial-velocities.svg
	# this requires opening the gui
	#inkscape --verb=FitCanvasToDrawing --verb=FileSave --verb=FileQuit *.svg

autobuild:
	sphinx-autobuild -j "auto" -D CHAPTER=$(CHAPTER) -b html --ignore "$(shell pwd)/_build/jupyter_execute/*" . _build/html/
