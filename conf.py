import os
import subprocess

# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

# -- Project information -----------------------------------------------------

project = 'Learn Multibody Dynamics'
html_title = project
copyright = '2022-2024, Jason K. Moore'
author = 'Jason K. Moore'
version = '0.2.dev0'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx_togglebutton',  # this has to be first so that the material css doesn't clobber its css
    'jupyter_sphinx',
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.graphviz',
    'sphinx.ext.imgconverter',  # converts the svg for latex builds
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.todo',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

ONGITHUB = "ONGITHUB" in os.environ

if ONGITHUB:
    commit_id = subprocess.check_output(['git', 'rev-parse', '--short',
                                         'HEAD']).strip().decode('ascii')
    version += '+' + commit_id

if not ONGITHUB:
    # Display TODO notes.
    todo_include_todos = True

# This configures sphinx to number figures and allow referencing them, if
# labeled, using :numref:`my_figure`.
numfig = True

# Sphinx >=4 default to MathJax v3, but v3 does not support wrapping lines. So
# force Sphinx to use v2 and config MathJax to wrap long lines.
mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@2/MathJax.js?config=TeX-AMS-MML_HTMLorMML"
mathjax2_config = {
    "HTML-CSS": {
        "linebreaks": {"automatic": True}
    }
}

# Setup intersphinx so that we can reference the SymPy documentation.
# :external:py:func:`~sympy.physics.vector.functions.dot` for example.
intersphinx_mapping = {
    'matplotlib': ('https://matplotlib.org/stable/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
    'sympy': ('https://docs.sympy.org/latest/', None),
    'py3js': ('https://pythreejs.readthedocs.io/en/stable', None),
}

# Options for Math
# These options influence Math notations.
# Set this option to True if you want all displayed math to be numbered. The
# default is False.

math_number_all = True

# Don't parse these files with Sphinx:
exclude_patterns = [
    'README.rst',
]

if "CHAPTER" in os.environ:
    CHAPTER = os.environ['CHAPTER']
else:
    CHAPTER = None

chapters = [
    'angular',
    'configuration',
    'differentiation',
    'energy',
    'eom',
    'generalized-forces',
    'holonomic-eom',
    'jupyter-python',
    'lagrange',
    'loads',
    'mass',
    'motion',
    'noncontributing',
    'nonholonomic-eom',
    'orientation',
    'simulation',
    'sympy',
    'tmt',
    'translational',
    'vectors',
    'visualization',
]

if (CHAPTER is not None) and (CHAPTER.lower() in chapters):
    chapters.remove(CHAPTER.lower())
    exclude_patterns += [c + '.rst' for c in chapters]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_material'
html_sidebars = {
    "**": ["logo-text.html", "globaltoc.html", "localtoc.html",
           "searchbox.html"]

}
html_theme_options = {
    'base_url': 'https://moorepants.github.io/learn-multibody-dynamics/',
    'color_primary': 'teal',
    'color_accent': 'deep-orange',  # hover color of hyperlinks
    'repo_name': 'Learn Multibody Dynamics',
    'repo_url': 'https://github.com/moorepants/learn-multibody-dynamics/',
    "logo_icon": "&#xe52f",
    'master_doc': False,  # Doesn't show duplicate title
    'nav_links': [{"href": "index", "internal": True, "title": "Home"}],
}
if ONGITHUB:
    # Takes too long to build locally for autobuild, so only do it in
    # production.
    html_theme_options['css_minify'] = False
    html_theme_options['html_minify'] = False
html_css_files = ['css/custom.css']  # seems to load after the material css

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Options for jupyter-sphinx

jupyter_execute_data_priority = [
    'application/vnd.jupyter.widget-view+json',
    'text/html',
    'text/latex',  # put latex before png so sympy displays mathjax
    'image/svg+xml',
    'image/png',
    'image/jpeg',
    'text/plain',
]

# latex build options
latex_engine = 'xelatex'  # handles unicode
latex_toplevel_sectioning = 'chapter'
