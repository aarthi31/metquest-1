from setuptools import setup, find_packages
import os

requirements = [
    'cobra',
    'numpy',
    'scipy',
    'networkx>=2.1',
    'python-libsbml',
]

try:
    with open('README.md') as handle:
        description = handle.read()
except:
    description = ''
datadir = os.path.join('metquest','example','data')
datafiles = [(d, [os.path.join(d,f) for f in files])
    for d, folders, files in os.walk(datadir)]

setup(
    name='metquest',
    version='0.1.29',
    packages=find_packages(),
    install_requires=requirements,
    setup_requires=[],
    scripts=['bin/metquest.sh'],
    author='Aarthi Ravikrishnan',
    author_email='aarthiravikrishnan@gmail.com',
    description='MetQuest: Enumerating all possible biosynthetic pathways in metabolic networks ',
    long_description=description,
    license="LGPL/GPL v2+",
    keywords='metabolism biology graph-theory pathways',
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v2'
            ' or later (LGPLv2+)',
        'License :: OSI Approved :: GNU General Public License v2'
            ' or later (GPLv2+)',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    platforms='GNU/Linux, Mac OS X >= 10.7, Microsoft Windows >= 7',
    data_files = datafiles,
    include_package_data = True

)

