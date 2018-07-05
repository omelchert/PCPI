import setuptools
from Cython.Build import cythonize
import numpy as np


PyPCPI_ext = cythonize(
    (setuptools.Extension('PyPCPI.layeredMedia.poissonIntegralSolver.poissonIntegral_cython.customPolarSolverMcml',
               sources=['PyPCPI/layeredMedia/poissonIntegralSolver/poissonIntegral_cython/customPolarSolverMcml.pyx'],
               include_dirs=[np.get_include()]),
     setuptools.Extension('PyPCPI.voxelizedMedia.poissonIntegralSolver.poissonIntegral_cython.customCartesianSolverMcxyz',
               sources=['PyPCPI/voxelizedMedia/poissonIntegralSolver/poissonIntegral_cython/customCartesianSolverMcxyz.pyx'],
               include_dirs=[np.get_include()]),
))


setuptools.setup(
    name='PyPCPI',
    version='1.0.0.dev1',
    author='Oliver Melchert',
    author_email='',
    description='A python module for optoacoustic signal generation via'
                '[P]olar [C]onvolution and [P]oisson [I]ntegral evaluation',
    url='https://github.com/omelchert/PCPI',
    packages=setuptools.find_packages(),
    ext_modules=PyPCPI_ext,
)
