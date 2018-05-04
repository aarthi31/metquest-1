
************
Installation
************

Use pip3 to install metquest from `PyPI: <https://pypi.python.org/pypi/metquest>`__

Via Python Package
==================

.. code:: bash

    pip3 install metquest


Direct installation
===================

1. Install `Python 3.4 or higher <https://www.python.org/downloads/>`__
2. Clone this repository to your computer using ``git`` or `download the
   repository <https://github.com/aarthi31/MetQuest/>`__ and decompress
   it. 
3. Navigate to the folder where metquest is downloaded and type

.. code:: bash

    python setup.py install

(Elevated ``sudo`` rights may be required depending on the platform. Replace python with python3 when multiple python distributions are found.)


Required python packages
========================
1. cobra >= 0.11.3
2. numpy >= 1.14.3
3. scipy
4. python-libsbml
5. networkx >= 2.1

Input
=====

Folder whose structure is as shown:

.. code-block:: text

    mainfolder/
    |-Example1/
    |   |-- SBML model(s) of metabolic networks # XML files of the metabolic networks(COBRA-compatible)
    |   |-- seed_mets.txt       # Text file containing the seed metabolites separated by a newline
    |   │-- source_mets.txt     # Text file containing the source metabolites separated by a newline
    |   |-- target_mets.txt     # Text file containing the target metabolites separated by a newline
    |   |-- cutoff.txt          # Text file containing the size cut-offs separated by a newline  
    |-Example2/
    |   ...

Kindly ensure that the SBML model has the field <model id> and the metabolites
are prefixed with the model identifiers, for instance, if the model identifier is 
'ecoli_core_model', and the seed metabolite is 'fum_c', the input text file
should contain ecoli_core_model fum_c

Running MetQuest
================

From command line
*****************
MetQuest can be directly run from the terminal as

.. code:: bash

    metquest.sh <path containing the input folder>


Navigate to the folder where metquest is installed and type

.. code:: bash

    python execute_metquest.py <path containing the input folder>

(Replace python with python3 when multiple python distributions are found.)


From python console
********************

.. code:: python

>>> import metquest
>>> metquest.execute_all_codes()


When prompted, enter the path containing the folder with all the data files

Running examples
****************

In the python console, type the following

.. code:: python

>>> import metquest
>>> metquest.example.run_this_example()


This will run the example files.

