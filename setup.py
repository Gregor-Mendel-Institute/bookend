from setuptools import setup
from setuptools.extension import Extension
try:
    from Cython.setuptools import build_ext
except:
    # If we couldn't import Cython, use the normal setuptools
    # and look for a pre-compiled .c file instead of a .pyx file
    from setuptools.command.build_ext import build_ext
    ext_modules = [
        Extension("bookend.core.cython_utils._assembly_utils", ["bookend/core/cython_utils/_assembly_utils.c"]),
        Extension("bookend.core.cython_utils._element_graph", ["bookend/core/cython_utils/_element_graph.c"]),
        Extension("bookend.core.cython_utils._fasta_utils", ["bookend/core/cython_utils/_fasta_utils.c"]),
        Extension("bookend.core.cython_utils._pq", ["bookend/core/cython_utils/_pq.c"]),
        Extension("bookend.core.cython_utils._rnaseq_utils", ["bookend/core/cython_utils/_rnaseq_utils.c"]),
    ]
else:
    # If we successfully imported Cython, look for a .pyx file
    ext_modules = [
        Extension("bookend.core.cython_utils._assembly_utils", ["bookend/core/cython_utils/_assembly_utils.pyx"]),
        Extension("bookend.core.cython_utils._element_graph", ["bookend/core/cython_utils/_element_graph.pyx"]),
        Extension("bookend.core.cython_utils._fasta_utils", ["bookend/core/cython_utils/_fasta_utils.pyx"]),
        Extension("bookend.core.cython_utils._pq", ["bookend/core/cython_utils/_pq.pyx"]),
        Extension("bookend.core.cython_utils._rnaseq_utils", ["bookend/core/cython_utils/_rnaseq_utils.pyx"]),
    ]

class CustomBuildExtCommand(build_ext):
    """build_ext command for use when numpy headers are needed."""
    def run(self):
        # Import numpy here, only when headers are needed
        import numpy
        # Add numpy headers to include_dirs
        self.include_dirs.append(numpy.get_include())
        # Call original build_ext command
        build_ext.run(self)

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="bookend_rna", # Replace with your own username
    version="1.1.0",
    author="Michael A. Schon",
    author_email="michael.schon@wur.nl",
    description="End-guided transcript assembler for short and long RNA-seq reads.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Gregor-Mendel-Institute/bookend",
    install_requires=['cython', 'numpy', 'pysam'],
    packages=['bookend', 'bookend.core'],
    cmdclass = {'build_ext': CustomBuildExtCommand},
    entry_points={'console_scripts': ['bookend = bookend.__main__:main']},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    keywords='transcriptome assembler bioinformatics rna sequencing',
    ext_modules = ext_modules,
    python_requires='>=3.6',
)
