"""
I need my modules to be a "package" for FireWorks to be able to detect it.
The code below is inspired from the Internet and is the minimal amount of work
which does the trick. This "package" is not meant for distribution, is still under
heavy development and is not guaranteed to work at all!
"""
from setuptools import setup

setup(name='VaspDrive',
        version='0.1',
        description='Collection of utilities to run my vasp calculations with Fireworks',
        url='https://github.com/rousseab/VaspDrive',
        author='Bruno Rousseau',
        author_email='rousseau.bruno@gmail.com',
        license='NONE',
        packages=['VaspDrive'],
        zip_safe=False)

