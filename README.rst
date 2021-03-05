=======================
Scipion - PySeg plugin
=======================

This plugin allows to use **PySeg - De novo analysis for cryo-electron tomography -**
(https://github.com/anmartinezs/pyseg_system) within the Scipion framework.

=====
Setup
=====

- **System pre-requisites:**

    1. Cmake 2.6.3+
    2. GSL (GNU Scientific Library)

- **Install this plugin in devel mode:**

Using the command line:

.. code-block::

    scipion3 installp -p local/path/to/scipion-em-pyseg --devel

Installation can be checked out running some tests (Important: TestPosRec requires Scipion plugin for Dynamo
to be installed --> https://github.com/scipion-em/scipion-em-dynamo):

.. code-block::

     scipion3 tests pyseg.tests.test_preseg_graphs_fils_picking.TestFromPresegToPicking

.. code-block::

    scipion3 tests pyseg.tests.test_pos_rec.TestPostRec

- **Contact information:**

If you experiment any problem, please contact us here: scipion-users@lists.sourceforge.net or open an issue
--> https://github.com/scipion-em/scipion-em-pyseg/issues

We'll be pleased to help.

*Scipion Team*


