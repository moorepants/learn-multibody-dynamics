====================
Install the Software
====================

If you would like to run the software on your own computer, follow these
instructions which recommend the use of a Conda_-based installation process.
The process below describes setting up a new base environment with the needed
software installed.

.. _Conda: https://en.wikipedia.org/wiki/Conda_(package_manager)

1) Miniconda
============

Miniconda is a stripped down version of Anaconda so that you can install only
what you desire. If you already have Miniconda (or Anaconda) on your computer,
you can skip this step or delete your prior Miniconda (or Anaconda) folder on
your computer to uninstall it. Download Miniconda for your operating system:

https://docs.conda.io/en/latest/miniconda.html

Install as a user, not an administrator, when asked. This will install the
package manager conda and configure your computer to use the Python installed
with Miniconda when you open a terminal or command prompt.

2) Create and Activate an Environment
=====================================

Open either the terminal (Linux/Mac) or the (Anaconda) command prompt (Windows)
and type the following series of commands followed each by the <enter> key to
execute the commands.

Create the environment with:

.. code-block:: bash

   conda create -c conda-forge -n learn-multibody-dynamics python=3.10

The ``-c conda-forge`` flag installs the packages from `Conda Forge`_. Conda
Forge is a community maintained collection of compatible software packages and
offers a larger number of packages than the default configuration.

.. _Conda Forge: https://conda-forge.org/

Now activate the environment:

.. code-block:: bash

   conda activate learn-multibody-dynamics

3) Install Packages
===================

Now you can install the packages that are required for executing the code in
this book with this command:

.. code-block:: bash

   conda install -c conda-forge ipympl ipython jupyter notebook matplotlib numpy pythreejs "scikits.odes" scipy "sympy>=1.11"

4) Open Jupyter Notebook
========================

To check that everything works, type the command to open Jupyter:

.. code-block:: bash

   jupyter notebook

Jupyter should open in your web browser and you should be able to run the
scripts and notebooks found on the other pages.

Software Versions
=================

This website was built with the following software versions:

.. jupyter-execute::

   import IPython
   IPython.__version__

.. jupyter-execute::

   import jupyter_sphinx
   jupyter_sphinx.__version__

.. jupyter-execute::

   import matplotlib
   matplotlib.__version__

.. jupyter-execute::

   import notebook
   notebook.__version__

.. jupyter-execute::

   import numpy
   numpy.__version__

.. jupyter-execute::

   import platform
   platform.python_version()

.. jupyter-execute::

   import pythreejs._version
   pythreejs._version.__version__

.. jupyter-execute::

   import pkg_resources
   pkg_resources.get_distribution("scikits.odes").version

.. jupyter-execute::

   import scipy
   scipy.__version__

.. jupyter-execute::

   import sphinx
   sphinx.__version__

.. jupyter-execute::

   import sphinx_material
   sphinx_material.__version__

.. jupyter-execute::

   import sphinx_togglebutton
   sphinx_togglebutton.__version__

.. jupyter-execute::

   import sympy
   sympy.__version__
