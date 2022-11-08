from setuptools import setup, find_packages

from codecs import open
from os import path

# The directory containing this file
HERE = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(HERE, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="PolCurveFit",
    version="1.2.1",
    description="library for the analysis of polarization curves",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://polcurvefit.readthedocs.io/",
    author="Meeke van Ede",
    author_email="meekevanede@gmail.com",
    license="MIT",
    classifiers=[
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Operating System :: OS Independent"
    ],
    packages=["PolCurveFit"],
    include_package_data=True,
    package_data={'': ['data/*.txt']},
    install_requires=["numpy","scipy","matplotlib","pandas"] 
)
