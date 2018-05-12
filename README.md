Welcome to MetQuest\'s documentation!
=====================================

MetQuest is a dynamic programming based algorithm for identifying all
possible pathways from metabolic networks between the source and the
target metabolites. MetQuest requires the genome-scale metabolic
reconstructions, set of seed, source and target metabolites and the
pathway length cut-off. MetQuest is compatible with Python 3 and is
OS-independent. An extensive documentation can be found [here](http://metquestdoc.readthedocs.io/).

Authors
=======

-   [Aarthi Ravikrishnan](https://github.com/aarthi31)
-   Meghana Nasre
-   [Karthik Raman](https://github.com/ramanlab)


Acknowledgments
===============

This work was supported by the Indian Institute of Technology Madras
grant ERP/1314/004/RESF/KARH to KR and the INSPIRE fellowship,
Department of Science and Technology, Government of India to AR.

Installation
------------

Use pip3 to install metquest from
[PyPI:](https://pypi.python.org/pypi/metquest)

Via Python Package
==================

``` {.sourceCode .bash}
pip3 install metquest
```

Direct installation
===================

1.  Install [Python 3.4 or higher](https://www.python.org/downloads/)
2.  Clone this repository to your computer using `git` or [download the
    repository](https://github.com/aarthi31/MetQuest/) and decompress
    it.
3.  Navigate to the folder where metquest is downloaded and type

``` {.sourceCode .bash}
python3 setup.py install
```

(Elevated `sudo` rights may be required depending on the platform)

Required python packages
========================

1.  cobra \>= 0.11.3
2.  numpy \>= 1.14.3
3.  scipy
4.  python-libsbml
5.  networkx \>= 2.1

Input
=====

Folder whose structure is as shown:

``` {.sourceCode .text}
mainfolder/
|-Example1/
|   |-- SBML model(s) of metabolic networks # XML files of the metabolic networks(COBRA-compatible)
|   |-- seed_mets.txt       # Text file containing the seed metabolites separated by a newline
|   |-- source_mets.txt     # Text file containing the source metabolites separated by a newline
|   |-- target_mets.txt     # Text file containing the target metabolites separated by a newline
|   |-- cutoff.txt          # Text file containing the size cut-offs separated by a newline  
|-Example2/
|   ...
```

Kindly ensure that the SBML model has the field \<model id\> and the
metabolites are prefixed with the model identifiers, for instance, if
the model identifier is \'ecoli\_core\_model\', and the seed metabolite
is \'fum\_c\', the input text file should contain ecoli\_core\_model
fum\_c

Running MetQuest
================

### From command line

MetQuest can be directly run from the terminal as

``` {.sourceCode .bash}
metquest.sh <path containing the input folder>
```

Navigate to the folder where metquest is installed and type

``` {.sourceCode .bash}
python3 execute_metquest.py <path containing the input folder>
```

### From python console

``` 
>>> import metquest
>>> metquest.execute_all_codes()
```

When prompted, enter the path containing the folder with all the data files

### Running examples

In the python console, type the following

```
>>> import metquest
>>> metquest.example.run_this_example()
```

This will run the example files.
