from setuptools import setup
from Cython.Build import cythonize

setup(
    name='Amplitude',
    ext_modules=cythonize("amplitude.pyx"),
    zip_safe=False,
)
