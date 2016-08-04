#!/usr/bin/env python3

from setuptools import setup

DESCRIPTION = (
    'Interpolate between stationary points on the potential energy surface'
)

CLASSIFIERS = [
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Topic :: Scientific/Engineering :: Chemistry',
    'Topic :: Scientific/Engineering :: Physics'
]

setup(
    name='interpPES',
    version='0.1.0',
    description=DESCRIPTION,
    long_description=open('README.md', 'r').read(),
    author='Tschijnmo TSCHAU',
    author_email='tschijnmotschau@gmail.com',
    url='https://github.com/tschijnmo/interpPES',
    license='MIT',
    packages=['interpPES'],
    install_requires=open('requirements.txt', 'r').read().split(),
    classifiers=CLASSIFIERS
)
