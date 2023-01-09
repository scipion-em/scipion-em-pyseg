=======================
Scipion - PySeg plugin
=======================

This plugin allows to use PySeg_ - De novo analysis for cryo-electron tomography - within the Scipion framework.

=====
Setup
=====

**System pre-requisites:**

    1. Cmake 2.6.3+. The intallation command in Ubuntu is:

    .. code-block::

        sudo apt-get install cmake

    2. GSL (GNU Scientific Library). In Ubuntu, the installation command is:

    .. code-block::

        sudo apt-get install libgsl-dev

    3. gcc/g++ version greater or equal to 5 and lower or equal to 7 (for DisPerSE_ compilation). In Ubuntu,
       the installation command is:

    .. code-block::

        sudo apt -y install gcc-5 g++-5

============
Installation
============
The plugin can be installed in user (stable) or developer (latest, may be unstable) mode:

**1. User (stable) version:**:

.. code-block::

    scipion3 installp -p scipion-em-pyseg

**2. Developer (latest, may be unstable) version:**:

* Clone the source code repository:

.. code-block::

    git clone https://github.com/scipion-em/scipion-em-pyseg.git
    
* Install:

.. code-block::

    scipion3 installp -p local/path/to/scipion-em-pyseg --devel
    
=========
Protocols
=========
The integrated protocols are:

1. pyseg - fils: Filters a MbGraphMCF (Mean Cumulative Function) object by extracting a filament network

2. pyseg - graphs: Analyzes a GraphMCF (Mean Cumulative Function) from a segmented membrane

3. pyseg - picking: Extracts particles from a filament network of a oriented single membrane graph

4. pyseg - 2D classification: Unsupervised and deterministic classification of membrane-bound particles

5. pyseg - posrec: post-process already reconstructed particles; rot angle randomization and membrane suppression

6. pyseg - preseg membranes: Segment membranes into membranes, inner surroundings and outer surroundings
    
=====
Tests
=====

The installation can be checked out running some tests (Important: TestPosRec requires the plugins scipion-em-xmipp_
and scipion-em-reliontomo_ to be installed:

.. code-block::

     scipion3 tests pyseg.tests.test_preseg_graphs_fils_picking.TestFromPresegToPicking

.. code-block::

    scipion3 tests pyseg.tests.test_pos_rec.TestPostRec
    
========
Tutorial
========
A tutorial about how to use PySeg within Scipion can be found here_.

==========
References
==========

* `Template-free detection and classification of heterogeneous membrane-bound complexes in cryo-electron tomograms. <http://doi.org/10.1038/s41592-019-0675-5>`_
  A. Martinez-Sanchez et al., Nature Methods, 2020.

===================
Contact information
===================

If you experiment any problem, please contact us here: scipion-users@lists.sourceforge.net or open an issue_.

We'll be pleased to help.

*Scipion Team*


.. _PySeg: https://github.com/anmartinezs/pyseg_system
.. _DisPerSE: http://www2.iap.fr/users/sousbie/web/html/indexd41d.html
.. _scipion-em-xmipp: https://github.com/I2PC/scipion-em-xmipp
.. _scipion-em-reliontomo: https://github.com/scipion-em/scipion-em-reliontomo
.. _issue: https://github.com/scipion-em/scipion-em-pyseg/issues
.. _here: https://scipion-em.github.io/docs/release-3.0.0/docs/user/denoising_mbSegmentation_pysegDirPicking/tomosegmemTV-pySeg-workflow.html#tomosegmemtv-pyseg-workflow
