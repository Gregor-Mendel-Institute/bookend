import setuptools
import numpy
from distutils.core import setup
from distutils.extension import Extension

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="bookend", # Replace with your own username
    version="0.0.1",
    author="Michael A. Schon",
    author_email="michael.schon@gmi.oeaw.ac.at",
    description="End-guided transcript assembler for short and long RNA-seq reads.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Gregor-Mendel-Institute/bookend",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    keywords='transcriptome assembler bioinformatics rna sequencing',
    install_requires=['cython', 'numpy', 'pysam'],
    # scripts=['bin/sort-elr', 'bin/collapse-elr', 'bin/gff3-to-bed'],
    entry_points={
        'console_scripts': ['bookend = bookend.bookend:main'],
    },
    ext_modules = [
        Extension("bookend.core.cython_utils._assembly_utils", ["bookend/core/cython_utils/_assembly_utils.c"], include_dirs=[numpy.get_include()]),
        Extension("bookend.core.cython_utils._element_graph", ["bookend/core/cython_utils/_element_graph.c"], include_dirs=[numpy.get_include()]),
        Extension("bookend.core.cython_utils._fasta_utils", ["bookend/core/cython_utils/_fasta_utils.c"]),
        Extension("bookend.core.cython_utils._pq", ["bookend/core/cython_utils/_pq.c"]),
        Extension("bookend.core.cython_utils._rnaseq_utils", ["bookend/core/cython_utils/_rnaseq_utils.c"], include_dirs=[numpy.get_include()]),
    ],
    python_requires='>=3.6',
)
