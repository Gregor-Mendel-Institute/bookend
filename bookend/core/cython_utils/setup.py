#cython: language_level=3
from distutils.core import setup
from Cython.Build import cythonize
import numpy

cythonize('_assembly_utils.pyx')
cythonize('_element_graph.pyx')
cythonize('_fasta_utils.pyx')
cythonize('_pq.pyx')
cythonize('_rnaseq_utils.pyx')

# setup(name="_assembly_utils", ext_modules=cythonize('_assembly_utils.pyx', ), include_dirs=[numpy.get_include()])
# setup(name="_element_graph", ext_modules=cythonize('_element_graph.pyx', ), include_dirs=[numpy.get_include()])
# setup(name="_fasta_utils", ext_modules=cythonize('_fasta_utils.pyx'),)
# setup(name="_pq", ext_modules=cythonize('_pq.pyx'),)
# setup(name="_rnaseq_utils", ext_modules=cythonize('_rnaseq_utils.pyx'), include_dirs=[numpy.get_include()])

