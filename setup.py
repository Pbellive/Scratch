import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "Scratch",
    version = "0.0.1",
    author = "Dave Marchant",
    author_email = "dave@compgeoinc.com",
    description = ("Code that doesn't yet have a home"),
    license = "MIT",
    keywords = "Geophysics,miscellaneous",
    url = "https://github.com/Pbellive/Scratch",
    packages=['Scratch'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT License",
    ],
) 
