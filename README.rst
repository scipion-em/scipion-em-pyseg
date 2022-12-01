=======================
Scipion - PySeg plugin
=======================

This plugin allows to use PySeg_ - De novo analysis for cryo-electron tomography - within the Scipion framework.

=====
Setup
=====

- **System pre-requisites:**

    1. Cmake 2.6.3+. The intallation command in Ubuntu is:

    .. code-block::

        sudo apt-get install cmake

    2. GSL (GNU Scientific Library). In Ubuntu, the installation command is:

    .. code-block::

        sudo apt-get install libgsl-dev

    3. gcc/g++ 5, 6, or 7 (for DisPerSE_ compilation). In Ubuntu,
       the installation command is:

    .. code-block::

        sudo apt -y install gcc-5 g++-5

- **Install this plugin in devel mode:**

Using the command line:

.. code-block::

    scipion3 installp -p local/path/to/scipion-em-pyseg --devel

Installation can be checked out running some tests (Important: TestPosRec requires the plugin scipion-em-dynamo_
to be installed and plugin scipion-em-reliontomo_ to be updated:

.. code-block::

     scipion3 tests pyseg.tests.test_preseg_graphs_fils_picking.TestFromPresegToPicking

.. code-block::

    scipion3 tests pyseg.tests.test_pos_rec.TestPostRec

- **Contact information:**

If you experiment any problem, please contact us here: scipion-users@lists.sourceforge.net or open an issue_.

We'll be pleased to help.

*Scipion Team*


.. _PySeg: https://github.com/anmartinezs/pyseg_system
.. _DisPerSE: http://www2.iap.fr/users/sousbie/web/html/indexd41d.html
.. _scipion-em-dynamo: https://github.com/scipion-em/scipion-em-dynamo
.. _scipion-em-reliontomo: https://github.com/scipion-em/scipion-em-reliontomo
.. _issue: https://github.com/scipion-em/scipion-em-pyseg/issues
