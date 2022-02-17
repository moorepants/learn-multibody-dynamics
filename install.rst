====================
Install the Software
====================

If you would like to run the software on your own computer, follow these
instructions.

1) Miniconda
============

If you already have Anaconda or Miniconda on your computer from another class,
you can skip this step or delete your prior Anaconda/Miniconda folder on your
computer to uninstall it. Anaconda or Miniconda are functionally the same, but
Anaconda comes with much more pre-installed software. Miniconda is a stripped
down version so that you can install only what you desire.

Download miniconda for your operating system:

https://docs.conda.io/en/latest/miniconda.html

Install as a user, not an administrator, when asked. This will install the
package manager conda and configure your computer to use the Python installed
with Miniconda when you open a terminal or command prompt.

2) Conda Forge
==============

Configure conda to download packages from `Conda Forge`_ as the primary
download source. Conda Forge is a community maintained collection of compatible
software packages and offers a larger number of packages than the default
configuration.

.. _Conda Forge: https://conda-forge.org/

Open either the terminal (Linux/Mac) or the Anaconda command prompt (Windows)
and type the following commands followed by the <enter> key to execute the
commands.

This first command adds Conda Forge as a download source for software
packages:

.. code-block:: bash

   conda config --add channels conda-forge

This ensures conda selects packages from Conda Forge first:

.. code-block:: bash

   conda config --set channel_priority strict

Now that Conda Forge packages are available, update everything that is already
installed with:

.. code-block:: bash

   conda update --all

This command could take some minutes, especially if you are using Anaconda
because it has many more package compatibilities to work out.

3) Install Packages
===================

Now you can install the packages that are required for these materials with
this command:

.. code-block:: bash

   conda install ipympl ipython ipywidgets jupyter matplotlib notebook numpy pythreejs scipy "sympy>=1.9"

4) Open Jupyter Notebook
========================

To check that everything works, type the command to open Jupyter:

.. code-block:: bash

   jupyter notebook

Jupyter should open in your web browser and you should be able to run the
scripts and notebooks found on the other pages.
