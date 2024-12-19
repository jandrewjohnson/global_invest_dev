from distutils.core import setup
from distutils.extension import Extension
from setuptools import setup, find_packages
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import os

# MAKE THIS AN EDITABLE INSTALL by running the following command in the terminal while in the directory of this setup.py file:
# pip install -e .

packages=find_packages()
include_package_data=True

# extensions = [
#         Extension(
#         "seals.seals_cython_functions",  # This corresponds to the Python import path
#         ["seals/seals_cython_functions.pyx"],  # Path to the .pyx file
#     )
# ]

setup(
    name = 'global_invest',
    packages = packages,
    version = '0.0.1',
    description = 'Global ecosystem service models',
    author = 'Justin Andrew Johnson',
    url = 'https://github.com/jandrewjohnson/global_invest_dev',
    keywords = ['geospatial', 'raster', 'shapefile', 'ecoystem services'],
    classifiers = [],
    # ext_modules=cythonize(extensions),
    # include_dirs=[numpy.get_include()],
    cmdclass={'build_ext': build_ext},
)
