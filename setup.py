from setuptools import setup, find_packages
import codecs
import os

VERSION = '0.0.10'
DESCRIPTION = 'Dbar Algorithm for EIT'
LONG_DESCRIPTION = 'This package complements the pyEIT package, providing a framework of the Dbar Algorithm for Electrical Impedance Tomography initially implemented by S. Siltanen in Matlab based on the theoretical of A. Nachman'

setup(
    name="py_dbar",
    version=VERSION,
    author="NablaIp",
    author_email = "<ivanpombo.eigen@gmail.com>",
    description=DESCRIPTION,
    long_description_content_type = "text/markdown",
    long_description=LONG_DESCRIPTION,
    packages = find_packages(),
    install_requires = ["numpy", "scipy", "pyamg" ],
    keywords=['python', 'EIT', 'DBar Algorithm'],
    classifiers=[
       'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
    ]

)
